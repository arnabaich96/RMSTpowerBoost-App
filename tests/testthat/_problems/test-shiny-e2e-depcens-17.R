# Extracted from test-shiny-e2e-depcens.R:17

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RMSTpowerBoostApp", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_no_shinytest2()
app <- shinytest2::AppDriver$new(
    app_dir = app_dir_path(),
    name = "e2e-depcens",
    seed = 2
  )
on.exit(app$stop(), add = TRUE)
app$wait_for_idle(timeout = 60000)
app$set_inputs(
    data_mode = "Upload",
    pilot_data_upload = shinytest2::upload_file(fixture_linear_csv())
  )
