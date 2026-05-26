#!/usr/bin/env Rscript
# Combine single-view comparator HPC outputs.

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      i <- i + 1L
      next
    }
    nm <- sub("^--", "", key)
    val <- if (i + 1L <= length(args) && !grepl("^--", args[[i + 1L]])) args[[i + 1L]] else TRUE
    out[[nm]] <- val
    i <- i + if (isTRUE(val)) 1L else 2L
  }
  out
}

bind_fill <- function(dfs) {
  if (length(dfs) == 0) return(data.frame())
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(x) {
    miss <- setdiff(all_cols, names(x))
    for (m in miss) x[[m]] <- NA
    x[, all_cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

read_with_file <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  x$source_file <- basename(path)
  x
}

parse_singleview_summary_name <- function(path) {
  b <- basename(path)
  data.frame(
    job_id = as.integer(sub("^job([0-9]+)_.*$", "\\1", b)),
    sample_name = sub("^.*_(n[0-9]+)_.*$", "\\1", b),
    nsample = as.integer(sub("^.*_n([0-9]+)_.*$", "\\1", b)),
    correction = sub("^.*_q([A-Za-z]+)grid.*$", "\\1", b),
    snr = suppressWarnings(as.numeric(sub("^.*_snr([0-9p.]+)_.*$", "\\1", b))),
    n_pathways = as.integer(sub("^.*_path([0-9]+)_.*$", "\\1", b)),
    tie = suppressWarnings(as.numeric(gsub("p", ".", sub("^.*_tie([0-9p.]+)_summary[.]csv$", "\\1", b)))),
    stringsAsFactors = FALSE
  )
}

add_summary_metrics <- function(x, path) {
  meta <- parse_singleview_summary_name(path)
  x <- cbind(meta[rep(1L, nrow(x)), , drop = FALSE], x)
  if (!("correction" %in% names(x)) && "fdr_method" %in% names(x)) x$correction <- x$fdr_method
  if ("correction" %in% names(x) && "fdr_method" %in% names(x)) {
    x$correction[is.na(x$correction) | !nzchar(x$correction)] <- x$fdr_method[is.na(x$correction) | !nzchar(x$correction)]
  }
  if (all(c("precision", "recall") %in% names(x))) {
    denom <- x$precision + x$recall
    x$f1_score <- ifelse(is.finite(denom) & denom > 0, 2 * x$precision * x$recall / denom, NA_real_)
  }
  x
}

combine_type <- function(results_dir, pattern, out_file, label, transform = NULL) {
  files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
  cat(label, "files:", length(files), "\n")
  if (length(files) == 0) return(invisible(data.frame()))
  out <- bind_fill(lapply(files, function(f) {
    tryCatch({
      x <- read_with_file(f)
      if (!is.null(transform)) x <- transform(x, f)
      x
    }, error = function(e) {
      data.frame(source_file = basename(f), read_error = conditionMessage(e), stringsAsFactors = FALSE)
    })
  }))
  write.csv(out, out_file, row.names = FALSE)
  cat("Wrote", label, ":", normalizePath(out_file, mustWork = FALSE), "\n")
  invisible(out)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
results_dir <- args$`results-dir` %||% args$results_dir %||% "singleview_comparator_results"
out_dir <- args$`out-dir` %||% args$out_dir %||% results_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

combine_type(
  results_dir,
  "_summary[.]csv$",
  file.path(out_dir, "combined_singleview_summary.csv"),
  "summary",
  transform = add_summary_metrics
)
combine_type(results_dir, "_detail[.]csv$", file.path(out_dir, "combined_singleview_detail.csv"), "detail")
combine_type(results_dir, "_mediators[.]csv$", file.path(out_dir, "combined_singleview_mediators.csv"), "mediators")
combine_type(results_dir, "_settings[.]csv$", file.path(out_dir, "combined_singleview_settings.csv"), "settings")

err_files <- list.files(results_dir, pattern = "_errors[.]txt$", full.names = TRUE)
if (length(err_files) > 0) {
  errors <- do.call(rbind, lapply(err_files, function(f) {
    txt <- readLines(f, warn = FALSE)
    data.frame(
      source_file = basename(f),
      has_error = any(nzchar(txt)),
      message = paste(txt[nzchar(txt)], collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(errors, file.path(out_dir, "combined_singleview_errors.csv"), row.names = FALSE)
}

cat("Done.\n")
