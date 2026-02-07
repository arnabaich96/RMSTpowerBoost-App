
test_that("shiny e2e: dependent censoring model runs", {
  skip_if_no_shinytest2()

  app <- safe_app_driver(
    app_dir = app_dir_path(),
    name = "e2e-depcens",
    seed = 2,
    load_timeout = 60000,
    timeout = 60000
  )
  on.exit(app$stop(), add = TRUE)

  app$wait_for_idle(timeout = 60000)

  app$set_inputs(data_mode = "Upload")
  app$click("confirm_step1")
  app$upload_file(pilot_data_upload = fixture_linear_csv())
  app$wait_for_idle(timeout = 60000)

  # Confirm upload cleaning step
  app$click("confirm_step2_data_upload")
  app$wait_for_idle(timeout = 60000)

  # Switch model and ensure required censoring covariates are set
  app$set_inputs(model_selection = "Dependent Censoring Model")
  app$set_inputs(
    time_var = "time",
    status_var = "status",
    arm_var = "arm"
  )
  app$wait_for_idle(timeout = 60000)

  # Choose a couple of available covariates for the censoring Cox model
  app$set_inputs(dc_linear_terms = c("age", "sex"))
  app$wait_for_idle(timeout = 60000)

  app$click("confirm_step2")
  app$wait_for_idle(timeout = 60000)

  app$click("run_analysis")
  app$wait_for_idle(timeout = 180000)

  testthat::expect_false(is.null(app$get_value(output = "results_plot")))
  app$set_inputs(main_tabs = "Run Log")
  app$wait_for_idle(timeout = 60000)
  app$wait_for_value(output = "run_log_summary_ui", timeout = 180000)
})
