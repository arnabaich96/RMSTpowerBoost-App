options(stringsAsFactors = FALSE)

repo_root <- function() {
  if (file.exists("app.R")) return(".")
  candidate <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(candidate, "app.R"))) return(candidate)
  ""
}

test_that("report input snapshot includes reproducibility-related fields", {
  root <- repo_root()
  if (!nzchar(root)) skip("Repo app.R not available in this test context.")
  app_env <- new.env(parent = globalenv())
  old_wd <- getwd()
  setwd(root)
  on.exit(setwd(old_wd), add = TRUE)
  app_path <- "app.R"
  expect_true(file.exists(app_path))
  sys.source(app_path, envir = app_env)
  fn <- app_env$report_inputs_builder
  expect_true(is.function(fn))
  x <- list(
    model_selection = "Linear IPCW Model",
    analysis_type = "Power",
    time_var = "time",
    status_var = "status",
    arm_var = "arm",
    strata_var = "strata",
    dc_linear_terms = c("x1"),
    calc_method = "Repeated",
    L = 365,
    alpha = 0.05,
    sample_sizes = "100,150",
    target_power = 0.8,
    sim_seed = 11,
    seed_reps = 22,
    R_reps = 500,
    data_mode = "Upload"
  )
  out <- fn(x)
  expect_named(out, c(
    "model_selection", "analysis_type", "time_var", "status_var", "arm_var",
    "strata_var", "dc_linear_terms", "calc_method", "L", "alpha",
    "sample_sizes", "target_power", "sim_seed", "seed_reps", "replications",
    "data_mode"
  ))
})

test_that("modernized app features are present in app source", {
  app_root <- repo_root()
  if (!nzchar(app_root)) skip("Repo app.R not available in this test context.")
  app_path <- file.path(app_root, "app.R")
  expect_true(file.exists(app_path))
  txt <- paste(readLines(app_path, warn = FALSE), collapse = "\n")
  
  expect_match(txt, "accordion\\(", perl = TRUE)
  expect_match(txt, "download_csv_template", perl = TRUE)
  expect_match(txt, "download_toy_pilot", perl = TRUE)
  expect_match(txt, "AFT \\(Lognormal\\)", perl = TRUE)
  expect_match(txt, "PH \\(Piecewise Exponential\\)", perl = TRUE)
  expect_match(txt, "run_checklist_ui", perl = TRUE)
  expect_match(txt, "sim_validation_ui", perl = TRUE)
  expect_match(txt, "auto_run_pending", perl = TRUE)
  expect_match(txt, "Reproducibility", perl = TRUE)
  expect_match(txt, "tabPanel\\(\\s*\"Pipeline\"", perl = TRUE)
  expect_match(txt, "tabPanel\\(\\s*\"About\"", perl = TRUE)
  expect_match(txt, "selected = \"Pipeline\"", perl = TRUE)
  expect_match(txt, "Bug a? Report|Report a Bug", perl = TRUE)
  expect_match(txt, "Coverage", perl = TRUE)
  expect_match(txt, "Citation / Contact", perl = TRUE)
})

test_that("pipeline html page exists with expected sections", {
  app_root <- repo_root()
  if (!nzchar(app_root)) skip("Repo pipeline.html not available in this test context.")
  html_path <- file.path(app_root, "www", "pipeline.html")
  expect_true(file.exists(html_path))
  html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  expect_match(html, "RMSTpowerBoost Workflow", perl = TRUE)
  expect_match(html, "Step 1b\\. Upload-Only Data Cleaning", perl = TRUE)
  expect_match(html, "Quick Troubleshooting", perl = TRUE)
})
