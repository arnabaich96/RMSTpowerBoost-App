# Extracted from test-shiny-e2e-linear.R:45

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RMSTpowerBoostApp", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_no_shinytest2()
app <- safe_app_driver(
    app_dir = app_dir_path(),
    name = "e2e-linear",
    seed = 1,
    load_timeout = 60000,
    timeout = 60000
  )
if (inherits(app, "shinytest2_unavailable")) {
    testthat::succeed()
    return()
  }
on.exit(app$stop(), add = TRUE)
app$wait_for_idle(timeout = 60000)
app$set_inputs(data_mode = "Upload")
app$click("confirm_step1")
app$upload_file(pilot_data_upload = fixture_linear_csv())
app$wait_for_idle(timeout = 60000)
app$click("confirm_step2_data_upload")
app$wait_for_idle(timeout = 60000)
app$set_inputs(
    time_var = "time",
    status_var = "status",
    arm_var = "arm",
    sample_sizes = "50,100"
  )
app$wait_for_idle(timeout = 60000)
app$click("confirm_step2")
app$wait_for_idle(timeout = 60000)
app$click("run_analysis")
app$wait_for_idle(timeout = 180000)
