options(stringsAsFactors = FALSE)

get_app_root <- function() {
  if (file.exists("app.R")) {
    "."
  } else {
    candidate <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(candidate, "app.R"))) candidate else ""
  }
}

read_uploaded_dataset_contract <- function(path, original_name) {
  ext <- tolower(tools::file_ext(original_name))
  if (identical(ext, "csv")) {
    return(utils::read.csv(path, check.names = FALSE))
  }
  if (ext %in% c("txt", "tsv")) {
    df <- tryCatch(utils::read.delim(path, check.names = FALSE), error = function(e) NULL)
    if (!is.null(df) && ncol(df) > 1) return(df)
    df2 <- tryCatch(utils::read.csv(path, check.names = FALSE), error = function(e) NULL)
    if (!is.null(df2)) return(df2)
    stop("Could not parse text upload as TSV or CSV.", call. = FALSE)
  }
  if (identical(ext, "rds")) {
    obj <- readRDS(path)
    if (is.data.frame(obj)) return(obj)
    if (is.matrix(obj)) return(as.data.frame(obj, check.names = FALSE))
    stop("RDS file must contain a data.frame or matrix.", call. = FALSE)
  }
  if (ext %in% c("rdata", "rda")) {
    e <- new.env(parent = emptyenv())
    loaded <- load(path, envir = e)
    if (!length(loaded)) stop("RData file does not contain any objects.", call. = FALSE)
    objs <- mget(loaded, envir = e, inherits = FALSE)
    pick_name <- names(objs)[which(vapply(objs, is.data.frame, logical(1)))[1]]
    if (is.na(pick_name) || !nzchar(pick_name)) {
      pick_name <- names(objs)[which(vapply(objs, is.matrix, logical(1)))[1]]
      if (is.na(pick_name) || !nzchar(pick_name)) {
        stop("RData file must contain a data.frame or matrix object.", call. = FALSE)
      }
    }
    obj <- objs[[pick_name]]
    if (is.matrix(obj)) obj <- as.data.frame(obj, check.names = FALSE)
    return(obj)
  }
  stop("Unsupported upload type. Use CSV, TXT/TSV, RDS, or RData.", call. = FALSE)
}

test_that("missing-data fixtures exist for all formats and manifest is complete", {
  root <- get_app_root()
  if (!nzchar(root)) skip("Repo fixtures are unavailable in this test context.")
  fixtures_dir <- file.path(root, "tests", "fixtures", "missing_pipeline")
  manifest_path <- file.path(fixtures_dir, "manifest.csv")
  expect_true(file.exists(manifest_path), info = "Run sim.data/generate_missing_pipeline_fixtures.R first.")

  manifest <- utils::read.csv(manifest_path, check.names = FALSE)
  expect_equal(nrow(manifest), 20L)
  expect_true(all(c("fixture_name", "format", "path", "nrow", "ncol", "missing_total", "missing_pct") %in% names(manifest)))

  fixtures <- c("fixture_nomiss", "fixture_mildmiss", "fixture_heavymiss", "fixture_mixed_edge")
  formats <- c("csv", "txt", "tsv", "rds", "rdata")
  for (fx in fixtures) {
    for (fmt in formats) {
      ext <- if (fmt == "rdata") "RData" else fmt
      expect_true(file.exists(file.path(fixtures_dir, paste0(fx, ".", ext))), info = paste("Missing fixture:", fx, fmt))
    }
  }
})

test_that("fixture parsers load all supported formats with stable shape", {
  root <- get_app_root()
  if (!nzchar(root)) skip("Repo fixtures are unavailable in this test context.")
  fixtures_dir <- file.path(root, "tests", "fixtures", "missing_pipeline")
  fixtures <- c("fixture_nomiss", "fixture_mildmiss", "fixture_heavymiss", "fixture_mixed_edge")
  formats <- c("csv", "txt", "tsv", "rds", "RData")

  for (fx in fixtures) {
    for (fmt in formats) {
      fpath <- file.path(fixtures_dir, paste0(fx, ".", fmt))
      df <- read_uploaded_dataset_contract(fpath, basename(fpath))
      expect_true(is.data.frame(df), info = paste("Not data.frame:", basename(fpath)))
      expect_equal(nrow(df), 100L, info = basename(fpath))
      expect_equal(ncol(df), 8L, info = basename(fpath))
      expect_true(all(c("time", "status", "arm", "strata", "x_num1", "x_num2", "x_cat1", "x_chr1") %in% names(df)))
    }
  }
})

test_that("missingness integrity is deterministic and heavy fixture has threshold-drop candidates", {
  root <- get_app_root()
  if (!nzchar(root)) skip("Repo fixtures are unavailable in this test context.")
  fixtures_dir <- file.path(root, "tests", "fixtures", "missing_pipeline")
  read_csv <- function(name) utils::read.csv(file.path(fixtures_dir, paste0(name, ".csv")), check.names = FALSE)

  nomiss <- read_csv("fixture_nomiss")
  mild <- read_csv("fixture_mildmiss")
  heavy <- read_csv("fixture_heavymiss")
  mixed <- read_csv("fixture_mixed_edge")

  expect_equal(sum(is.na(nomiss)), 0L)
  expect_equal(sum(is.na(mild)), 23L)
  expect_equal(sum(is.na(heavy)), 183L)
  expect_equal(sum(is.na(mixed)), 24L)
  expect_gt(sum(is.na(mild$strata)), 0L)
  expect_gt(sum(is.na(heavy$strata)), 0L)
  expect_gt(sum(is.na(mixed$strata)), 0L)

  missing_totals <- c(
    nomiss = sum(is.na(nomiss)),
    mild = sum(is.na(mild)),
    heavy = sum(is.na(heavy)),
    mixed = sum(is.na(mixed))
  )
  expect_identical(names(which.max(missing_totals)), "heavy")

  missing_pct_cols <- colMeans(is.na(heavy))
  drop_candidates <- names(missing_pct_cols)[missing_pct_cols > 0.35]
  expect_true(all(c("x_num1", "x_num2") %in% drop_candidates))
})

test_that("forward-compat: cleaning mode gating and final dataset pointer tests are pending", {
  skip("Pending implementation: missing-data cleaning UI/server pipeline.")
  expect_true(FALSE)
})
