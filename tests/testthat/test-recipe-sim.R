options(stringsAsFactors = FALSE)

r_dir <- if (dir.exists("R")) "R" else file.path("..", "..", "R")
invisible(lapply(sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE)), source))

test_that("validate_recipe normalizes PH exponential and sets defaults", {
  rec <- list(
    n = 10,
    covariates = list(defs = list()),
    treatment = list(assignment = "randomization", allocation = "1:1"),
    event_time = list(
      model = "ph_exp",
      baseline = list(rate = 0.1),
      effects = list(intercept = 0, treatment = 0)
    ),
    censoring = list(mode = "target_overall", target = 0.2, admin_time = 5)
  )

  v <- validate_recipe(rec)
  expect_equal(v$event_time$model, "cox_pwexp")
  expect_true(is.list(v$event_time$baseline))
  expect_true(is.numeric(v$event_time$baseline$rates))
})

test_that("simulate_from_recipe handles explicit censoring and logistic PS", {
  defs <- list(
    list(name = "x1", type = "continuous", dist = "normal",
         params = list(mean = 0, sd = 1), transform = list("center(0)", "scale(2)")),
    list(name = "grp", type = "categorical", dist = "categorical",
         params = list(prob = c(0.5, 0.5), labels = c("g1", "g2")))
  )

  rec <- list(
    n = 30,
    covariates = list(defs = defs),
    treatment = list(
      assignment = "logistic_ps",
      ps_model = list(formula = "~ x1")
    ),
    event_time = list(
      model = "aft_weibull",
      baseline = list(shape = 1.2, scale = 2),
      effects = list(intercept = 0, treatment = -0.1, formula = "~ x1", beta = c(0, 0.2)),
      frailty = list(type = "lognormal", var = 0.2, group = "grp")
    ),
    censoring = list(
      mode = "explicit",
      administrative = list(time = 3),
      random = list(dist = "exponential", params = list(rate = 0.1)),
      dependent = list(formula = "~ x1", beta = c(0, 0.1), base = 0.05)
    ),
    seed = 1
  )

  dat <- simulate_from_recipe(rec)
  expect_true(is.data.frame(dat))
  expect_true(all(c("time", "status", "arm", "x1", "grp") %in% names(dat)))
  expect_true(is.numeric(attr(dat, "achieved_censoring")))
})

test_that("simulate_from_recipe supports gamma frailty with PH model", {
  defs <- list(
    list(name = "grp", type = "categorical", dist = "categorical",
         params = list(prob = c(0.4, 0.6), labels = c("a", "b")))
  )

  rec <- list(
    n = 20,
    covariates = list(defs = defs),
    treatment = list(assignment = "randomization", allocation = "1:1"),
    event_time = list(
      model = "ph_weibull",
      baseline = list(shape = 1.2, scale = 1.5),
      effects = list(intercept = 0, treatment = 0.1),
      frailty = list(type = "gamma", var = 0.3, group = "grp")
    ),
    censoring = list(
      mode = "explicit",
      administrative = list(time = 2),
      random = list(dist = "exponential", params = list(rate = 0.05))
    )
  )

  dat <- simulate_from_recipe(rec)
  expect_equal(nrow(dat), 20)
  expect_true(is.numeric(attr(dat, "achieved_censoring")))
})

test_that("recipe_grid expands paths and gen_covariates handles factors", {
  defs <- list(
    list(name = "cat", type = "categorical", dist = "categorical",
         params = list(prob = c(0.5, 0.5), labels = c("A", "B"))),
    list(name = "ord", type = "ordinal", dist = "ordinal",
         params = list(prob = c(0.2, 0.8), labels = c("L", "H")))
  )

  rec <- list(
    n = 12,
    covariates = list(defs = defs),
    treatment = list(assignment = "randomization", allocation = "1:1"),
    event_time = list(
      model = "aft_lognormal",
      baseline = list(mu = 2.2, sigma = 0.5),
      effects = list(intercept = 0, treatment = -0.2)
    ),
    censoring = list(mode = "target_overall", target = 0.2, admin_time = 2)
  )

  grid <- recipe_grid(rec, vary = list(n = c(10, 12), "event_time.effects.treatment" = c(-0.1, -0.2)))
  expect_equal(length(grid), 4)

  X <- gen_covariates(5, list(defs = defs))
  expect_true(is.factor(X$cat))
  expect_true(is.ordered(X$ord))
})

test_that("validate_recipe catches malformed effects and frailty settings", {
  base <- list(
    n = 8,
    covariates = list(defs = list(list(
      name = "grp", type = "categorical", dist = "categorical",
      params = list(prob = c(0.5, 0.5), labels = c("a", "b"))
    ))),
    event_time = list(
      model = "aft_lognormal",
      baseline = list(mu = 1.5, sigma = 0.8),
      effects = list(intercept = 0, treatment = 0)
    ),
    censoring = list(mode = "target_overall", target = 0.2, admin_time = 3)
  )

  rec_bad_beta <- base
  rec_bad_beta$event_time$effects$formula <- "~ grp"
  expect_error(validate_recipe(rec_bad_beta), "effects\\$beta missing")

  rec_bad_formula <- base
  rec_bad_formula$event_time$effects$beta <- c(0, 0.1)
  expect_error(validate_recipe(rec_bad_formula), "effects\\$formula missing")

  rec_bad_frailty <- base
  rec_bad_frailty$event_time$frailty <- list(type = "bad", var = 0.1, group = "grp")
  expect_error(validate_recipe(rec_bad_frailty), "frailty\\$type")
})

test_that("allocation parsing and treatment assignment error branches are covered", {
  expect_equal(.parse_allocation("1:3")[["p1"]], 0.25)
  expect_error(.parse_allocation("1-3"), "allocation must be")

  X <- data.frame(x1 = rnorm(5))
  tr <- list(assignment = "logistic_ps", ps_model = list(formula = "~ x1", beta = c(0.1)))
  expect_error(.assign_treatment(5, X, tr), "beta length mismatch")
})

test_that("linear predictor and frailty edge/error paths are exercised", {
  X <- data.frame(x1 = c(0, 1), grp = c("a", "b"))
  expect_error(
    .build_lp(list(intercept = 0, treatment = 0, covariates = list(missing = 1)), X, arm = c(0, 1)),
    "unknown covariate"
  )

  eta <- c(0, 0)
  expect_error(
    .apply_frailty(eta, model = "aft_weibull",
                   frailty = list(type = "gamma", var = 0.2, group = "grp"),
                   groups = list(grp = c("a", "b"))),
    "PH/Cox models only"
  )
  expect_error(
    .apply_frailty(eta, model = "ph_weibull",
                   frailty = list(type = "lognormal", var = 0.2, group = "grp"),
                   groups = list()),
    "not found"
  )
})

test_that("piecewise simulator/time model and censoring solver branches are covered", {
  set.seed(99)
  expect_error(.sim_pwexp(rates = c(0.1, 0.2), cuts = numeric(0), lp = c(0, 0)), "must be m-1")
  t1 <- .sim_pwexp(rates = c(0.2), cuts = numeric(0), lp = c(0, 0, 0))
  expect_true(all(is.finite(t1) & t1 > 0))

  t2 <- .sim_time("cox_pwexp", baseline = list(rates = c(0.1, 0.2), cuts = 1), eta = c(0, 0, 0), n = 3)
  expect_length(t2, 3)
  expect_true(all(is.finite(t2) & t2 > 0))

  # If admin censoring alone already exceeds target, solver should return 0.
  rate <- .solve_rate_for_target(T_event = c(10, 12, 14), target = 0.2, admin_time = 1, tol = 0.01)
  expect_equal(rate, 0)
})

test_that("simulate_from_recipe explicit censoring validation errors are covered", {
  rec_random_bad <- list(
    n = 20,
    covariates = list(defs = list(list(name = "x1", type = "continuous", dist = "normal",
                                       params = list(mean = 0, sd = 1)))),
    treatment = list(assignment = "randomization", allocation = "1:1"),
    event_time = list(
      model = "aft_lognormal",
      baseline = list(mu = 2, sigma = 0.5),
      effects = list(intercept = 0, treatment = 0)
    ),
    censoring = list(
      mode = "explicit",
      random = list(dist = "weibull", params = list(shape = 1, scale = 2))
    )
  )
  expect_error(simulate_from_recipe(rec_random_bad, seed = 3), "exponential")

  rec_dep_bad <- rec_random_bad
  rec_dep_bad$censoring <- list(
    mode = "explicit",
    dependent = list(formula = "~ x1", beta = c(0, 0.1, 0.2), base = 0.05)
  )
  expect_error(simulate_from_recipe(rec_dep_bad, seed = 3), "beta length mismatch")
})
