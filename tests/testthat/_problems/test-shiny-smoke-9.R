# Extracted from test-shiny-smoke.R:9

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RMSTpowerBoostApp", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_no_shinytest2()
app <- shinytest2::AppDriver$new(
    app_dir = app_dir_path(),
    name = "smoke",
    seed = 1
  )
