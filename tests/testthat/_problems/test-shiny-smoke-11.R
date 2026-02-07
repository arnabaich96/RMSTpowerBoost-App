# Extracted from test-shiny-smoke.R:11

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RMSTpowerBoostApp", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_no_shinytest2()
app <- safe_app_driver(
    app_dir = app_dir_path(),
    name = "smoke",
    seed = 1,
    load_timeout = 60000,
    timeout = 60000
  )
