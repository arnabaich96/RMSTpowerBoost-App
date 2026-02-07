# Opt-in end-to-end testing for the Shiny app (browser-driven).
#
# Why opt-in?
# - Requires shinytest2 + a headless browser (Chromium via chromote)
# - Slower than unit tests
#
# Enable locally / in CI with:
#   Sys.setenv(RUN_SHINYTEST2 = "true")

skip_if_no_shinytest2 <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("shinytest2")

  flag <- Sys.getenv("RUN_SHINYTEST2", unset = "")
  if (!nzchar(flag) || tolower(flag) %in% c("0", "false", "no")) {
    testthat::skip("Set RUN_SHINYTEST2=true to run shinytest2 E2E tests.")
  }
}

app_dir_path <- function() {
  # Directory containing inst/app.R and inst/www/*
  testthat::test_path("..", "..", "inst")
}

fixture_linear_csv <- function() {
  testthat::test_path("..", "..", "sim.data", "linear_ipcw_data.csv")
}
