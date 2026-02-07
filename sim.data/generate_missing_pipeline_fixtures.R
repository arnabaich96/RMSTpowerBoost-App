#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

get_app_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, cmd_args, value = TRUE)
  if (length(hit)) {
    script_path <- normalizePath(sub(file_arg, "", hit[[1]]), winslash = "/", mustWork = FALSE)
    return(dirname(dirname(script_path)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

build_base_df <- function(n = 100L) {
  data.frame(
    time = round(seq(0.5, by = 0.4, length.out = n), 3),
    status = rep(c(1L, 0L), length.out = n),
    arm = rep(c(0L, 1L), each = n / 2L),
    strata = rep(c("A", "B", "C"), length.out = n),
    x_num1 = round(rnorm(n, mean = 0.5, sd = 1.1), 3),
    x_num2 = round(rnorm(n, mean = 10, sd = 2.5), 3),
    x_cat1 = factor(rep(c("low", "mid", "high"), length.out = n), levels = c("low", "mid", "high")),
    x_chr1 = rep(c("alpha", "beta", "gamma", "delta", "omega"), length.out = n)
  )
}

make_fixtures <- function() {
  set.seed(20260207)
  n <- 100L

  fixture_nomiss <- build_base_df(n)

  fixture_mildmiss <- build_base_df(n)
  fixture_mildmiss$x_num1[c(3, 17, 41, 66, 92)] <- NA
  fixture_mildmiss$x_cat1[c(5, 8, 48, 53, 79)] <- NA
  fixture_mildmiss$time[c(9, 58, 91)] <- NA
  fixture_mildmiss$status[c(11, 63, 99)] <- NA
  fixture_mildmiss$x_chr1[c(14, 20, 35, 76)] <- NA
  fixture_mildmiss$strata[c(12, 44, 90)] <- NA

  fixture_heavymiss <- build_base_df(n)
  fixture_heavymiss$x_num1[1:50] <- NA
  fixture_heavymiss$x_num2[21:70] <- NA
  fixture_heavymiss$x_cat1[seq(2, 40, by = 2)] <- NA
  fixture_heavymiss$x_chr1[seq(1, 50, by = 2)] <- NA
  fixture_heavymiss$status[71:80] <- NA
  fixture_heavymiss$time[81:95] <- NA
  fixture_heavymiss$arm[96:100] <- NA
  fixture_heavymiss$strata[c(11, 22, 33, 44, 55, 66, 77, 88)] <- NA

  fixture_mixed_edge <- build_base_df(n)
  fixture_mixed_edge$status <- ifelse(fixture_mixed_edge$status == 1, "event", "censored")
  fixture_mixed_edge$arm <- ifelse(fixture_mixed_edge$arm == 1, "treatment", "control")
  fixture_mixed_edge$status[c(4, 18, 46, 75)] <- NA
  fixture_mixed_edge$arm[c(7, 54)] <- NA
  fixture_mixed_edge$strata[c(25, 50, 75, 100)] <- NA
  fixture_mixed_edge$x_cat1[c(10, 11, 12, 70, 71)] <- NA
  fixture_mixed_edge$x_chr1[c(6, 19, 64, 89)] <- NA
  fixture_mixed_edge$x_num2[c(15, 20, 21, 67, 88)] <- NA

  list(
    fixture_nomiss = fixture_nomiss,
    fixture_mildmiss = fixture_mildmiss,
    fixture_heavymiss = fixture_heavymiss,
    fixture_mixed_edge = fixture_mixed_edge
  )
}

write_fixture_formats <- function(df, fixture_name, out_dir) {
  files <- list(
    csv = file.path(out_dir, paste0(fixture_name, ".csv")),
    txt = file.path(out_dir, paste0(fixture_name, ".txt")),
    tsv = file.path(out_dir, paste0(fixture_name, ".tsv")),
    rds = file.path(out_dir, paste0(fixture_name, ".rds")),
    rdata = file.path(out_dir, paste0(fixture_name, ".RData"))
  )

  utils::write.csv(df, files$csv, row.names = FALSE)
  utils::write.table(df, files$txt, sep = "\t", row.names = FALSE, quote = TRUE)
  utils::write.table(df, files$tsv, sep = "\t", row.names = FALSE, quote = TRUE)
  saveRDS(df, files$rds)
  obj <- df
  assign(fixture_name, obj)
  save(list = fixture_name, file = files$rdata)
  rm(list = fixture_name, inherits = TRUE)

  files
}

main <- function() {
  app_root <- get_app_root()
  out_dir <- file.path(app_root, "tests", "fixtures", "missing_pipeline")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fixtures <- make_fixtures()
  manifest_rows <- list()

  for (nm in names(fixtures)) {
    df <- fixtures[[nm]]
    files <- write_fixture_formats(df, nm, out_dir)
    miss_total <- sum(is.na(df))
    miss_pct <- round(100 * miss_total / (nrow(df) * ncol(df)), 3)

    for (fmt in names(files)) {
      rel_path <- file.path("tests", "fixtures", "missing_pipeline", basename(files[[fmt]]))
      manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
        fixture_name = nm,
        format = fmt,
        path = rel_path,
        nrow = nrow(df),
        ncol = ncol(df),
        missing_total = miss_total,
        missing_pct = miss_pct,
        stringsAsFactors = FALSE
      )
    }
  }

  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(out_dir, "manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE)

  cat("Missing-data fixtures generated:\n")
  print(manifest)
  cat("\nManifest written to:", manifest_path, "\n")
}

main()
