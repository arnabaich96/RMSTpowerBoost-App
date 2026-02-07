source("renv/activate.R")
renv::autoload()
# This app sources R/*.R explicitly inside app.R; disable Shiny's automatic
# R/ autoload to avoid package-dir warning and double sourcing.
options(shiny.autoload.r = FALSE)
if (dir.exists(".git")) {
  try({
    message("Running git pull...")
    result <- if (.Platform$OS.type == "windows") {
      shell("git pull", intern = TRUE)
    } else {
      system("git pull", intern = TRUE)
    }
    cat(result, sep = "\n")
  }, silent = TRUE)
}
