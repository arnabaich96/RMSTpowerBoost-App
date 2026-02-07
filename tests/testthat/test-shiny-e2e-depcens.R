
test_that("shiny e2e: dependent censoring model runs", {
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
  app$wait_for_idle(timeout = 60000)

  # Switch model and ensure required censoring covariates are set
  app$set_inputs(model_selection = "Dependent Censoring Model")
  app$wait_for_idle(timeout = 60000)

  # Choose a couple of available covariates for the censoring Cox model
  app$set_inputs(dc_linear_terms = c("age", "sex"))
  app$wait_for_idle(timeout = 60000)

  app$click("confirm_step2")
  app$wait_for_idle(timeout = 60000)

  app$click("run_analysis")
  app$wait_for_idle(timeout = 180000)

  testthat::expect_false(is.null(app$get_value(output = "results_plot")))
  testthat::expect_false(is.null(app$get_value(output = "run_log_summary_ui")))
})
