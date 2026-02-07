
test_that("shiny e2e: upload pilot data and run default analysis", {
  skip_if_no_shinytest2()

  app <- shinytest2::AppDriver$new(
    app_dir = app_dir_path(),
    name = "e2e-linear",
    seed = 1
  )
  on.exit(app$stop(), add = TRUE)

  app$wait_for_idle(timeout = 60000)

  # Upload the built-in fixture (no missing values)
  app$set_inputs(
    data_mode = "Upload",
    pilot_data_upload = shinytest2::upload_file(fixture_linear_csv())
  )
  app$wait_for_idle(timeout = 60000)

  # Step 2 confirm should unlock Step 3 controls
  app$click("confirm_step2")
  app$wait_for_idle(timeout = 60000)

  # Run analysis
  app$click("run_analysis")
  app$wait_for_idle(timeout = 180000)

  # Check primary outputs exist
  testthat::expect_false(is.null(app$get_value(output = "results_plot")))
  testthat::expect_false(is.null(app$get_value(output = "survival_plotly_output")))
  testthat::expect_false(is.null(app$get_value(output = "key_results_ui")))
})
