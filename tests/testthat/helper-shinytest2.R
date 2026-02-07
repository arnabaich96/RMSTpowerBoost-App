# Opt-in end-to-end testing for the Shiny app (browser-driven).
#
# Why opt-in?
# - Requires shinytest2 + a headless browser (Chromium via chromote)
# - Slower than unit tests
#
# Enable locally / in CI with:
#   Sys.setenv(RUN_SHINYTEST2 = "true")

skip_if_no_shinytest2 <- function() {
  # Force non-CRAN behavior for local/CI runs
  Sys.setenv(
    NOT_CRAN = "true",
    SHINYTEST2_APP_DRIVER_TEST_ON_CRAN = "1",
    TESTTHAT_PARALLEL = "false"
  )
  if (!requireNamespace("shinytest2", quietly = TRUE)) {
    stop("shinytest2 is required for E2E tests but is not installed.", call. = FALSE)
  }
  init_chromote_for_tests()
}

app_dir_path <- function() {
  # Prefer source repo if available; otherwise fall back to package root
  root <- testthat::test_path("..", "..")
  src <- find_source_app_dir()
  if (nzchar(src)) return(src)
  if (file.exists(file.path(root, "app.R"))) return(root)
  root
}

find_source_app_dir <- function() {
  env_path <- Sys.getenv("RMSTPOWERBOOSTAPP_SOURCE", unset = "")
  if (nzchar(env_path) && file.exists(file.path(env_path, "app.R"))) {
    return(normalizePath(env_path, winslash = "/", mustWork = FALSE))
  }
  # Common local repo locations (fast checks only)
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE", unset = ""), "Git-Repos", "RMSTSpowerBoost", "RMSTpowerBoost-App"),
    file.path(Sys.getenv("USERPROFILE", unset = ""), "Repos", "RMSTSpowerBoost", "RMSTpowerBoost-App"),
    "d:/Git-Repos/RMSTSpowerBoost/RMSTpowerBoost-App"
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(file.path(p, "app.R"))) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  ""
}

fixture_linear_csv <- function() {
  candidates <- c(
    testthat::test_path("..", "fixtures", "e2e", "linear_ipcw_data.csv"),
    testthat::test_path("..", "..", "sim.data", "linear_ipcw_data.csv")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(found[[1]])
  stop(
    "Missing test fixture 'linear_ipcw_data.csv'. Checked: ",
    paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

safe_app_driver <- function(...) {
  tryCatch(
    shinytest2::AppDriver$new(...),
    error = function(e) {
      # Retry once with a fresh chromote config
      init_chromote_for_tests()
      shinytest2::AppDriver$new(...)
    }
  )
}

init_chromote_for_tests <- function() {
  if (!requireNamespace("chromote", quietly = TRUE)) return(invisible(FALSE))
  chrome <- chromote::find_chrome()
  if (nzchar(chrome)) Sys.setenv(CHROMOTE_CHROME = chrome)
  # Prefer modern headless for recent Chrome versions and isolate a temp profile
  prof <- file.path(tempdir(), "chromote-profile")
  if (!dir.exists(prof)) dir.create(prof, recursive = TRUE, showWarnings = FALSE)
  chromote::set_chrome_args(c(
    "--headless=new",
    "--disable-gpu",
    "--force-color-profile=srgb",
    "--disable-extensions",
    "--mute-audio",
    "--no-first-run",
    "--no-default-browser-check",
    "--remote-allow-origins=*",
    paste0("--user-data-dir=", prof)
  ))
  invisible(TRUE)
}
