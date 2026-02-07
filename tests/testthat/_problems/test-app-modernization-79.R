# Extracted from test-app-modernization.R:79

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RMSTpowerBoostApp", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
repo_root <- function() {
  if (file.exists("app.R")) return(".")
  candidate <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(candidate, "app.R"))) return(candidate)
  ""
}

# test -------------------------------------------------------------------------
app_root <- repo_root()
if (!nzchar(app_root)) skip("Repo pipeline.html not available in this test context.")
html_path <- file.path(app_root, "www", "pipeline.html")
expect_true(file.exists(html_path))
html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
expect_match(html, "RMSTpowerBoost Workflow", perl = TRUE)
