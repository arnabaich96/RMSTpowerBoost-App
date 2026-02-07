
test_that("shiny smoke: app launches", {
  skip_if_no_shinytest2()

  app <- shinytest2::AppDriver$new(
    app_dir = app_dir_path(),
    name = "smoke",
    seed = 1
  )
  on.exit(app$stop(), add = TRUE)

  app$wait_for_idle(timeout = 60000)

  # Minimal assertion: key UI output exists (pipeline panel)
  val <- app$get_value(output = "pipeline_page_ui")
  testthat::expect_false(is.null(val))
})
