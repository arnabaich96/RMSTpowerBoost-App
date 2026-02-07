
test_that("shiny e2e: upload pilot data and run default analysis", {
  skip_if_no_shinytest2()

  app <- safe_app_driver(
    app_dir = app_dir_path(),
    name = "e2e-linear",
    seed = 1,
    load_timeout = 60000,
    timeout = 60000
  )
  on.exit(app$stop(), add = TRUE)

  app$wait_for_idle(timeout = 60000)

  # Upload the built-in fixture (no missing values)
  app$set_inputs(data_mode = "Upload")
  app$click("confirm_step1")
  app$upload_file(pilot_data_upload = fixture_linear_csv())
  app$wait_for_idle(timeout = 60000)

  # Confirm upload cleaning step
  app$click("confirm_step2_data_upload")
  app$wait_for_idle(timeout = 60000)

  # Map required fields (auto-map fallback if needed)
  app$set_inputs(
    time_var = "time",
    status_var = "status",
    arm_var = "arm"
  )
  app$wait_for_idle(timeout = 60000)

  # Confirm Step 3 (model & mapping)
  app$click("confirm_step2")
  app$wait_for_idle(timeout = 60000)

  # Run analysis
  app$click("run_analysis")
  app$wait_for_idle(timeout = 180000)

  # Check primary outputs exist
  testthat::expect_false(is.null(app$get_value(output = "results_plot")))
  app$set_inputs(main_tabs = "KM Plot")
  app$wait_for_idle(timeout = 60000)
  app$wait_for_value(output = "survival_plotly_output", timeout = 180000)
  app$set_inputs(main_tabs = "Analysis")
  app$wait_for_idle(timeout = 60000)
  app$wait_for_value(output = "key_results_ui", timeout = 180000)
})
