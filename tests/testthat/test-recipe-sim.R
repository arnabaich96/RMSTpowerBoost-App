options(stringsAsFactors = FALSE)

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
