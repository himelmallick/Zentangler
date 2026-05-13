#!/usr/bin/env Rscript
# Generate one paired SIMMBA/InterSIM simulation cache object for one grid row.
#
# The saved RDS is intentionally method-agnostic. Zentangler, HIMA, HIMA2, and
# multimedia runners should all read the same jobXXXX_sim.rds so benchmarking is
# paired on the exact same simulated data and truth labels.

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

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.numeric(x)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

row_value <- function(row, name, default = NULL) {
  if (!(name %in% colnames(row))) return(default)
  val <- row[[name]][1]
  if (length(val) == 0 || is.na(val) || !nzchar(as.character(val))) return(default)
  val
}

load_zentangler_sources <- function(package_dir = NULL) {
  if (!is.null(package_dir) && nzchar(package_dir) && file.exists(file.path(package_dir, "DESCRIPTION"))) {
    r_files <- sort(list.files(file.path(package_dir, "R"), pattern = "[.]R$", full.names = TRUE))
    for (f in r_files) source(f)
    return(invisible(TRUE))
  }
  suppressPackageStartupMessages(library(Zentangler))
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$grid)) stop("Missing --grid")
if (is.null(args$`task-id`)) stop("Missing --task-id")

grid_file <- args$grid
out_dir <- args$`out-dir` %||% args$out_dir %||% "simmba_cache"
package_dir <- args$`package-dir` %||% args$package_dir %||% NULL
task_offset <- as.integer(args$`task-offset` %||% args$task_offset %||% 0L)
array_task_id <- as.integer(args$`task-id`)
task_id <- array_task_id + task_offset
overwrite <- as_bool(args$overwrite %||% FALSE, default = FALSE)

grid <- utils::read.csv(grid_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop("effective task-id must be between 1 and ", nrow(grid), ". Got ", task_id)
}

row <- grid[task_id, , drop = FALSE]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sim_id <- as.integer(row_value(row, "sim_id", row_value(row, "job_id", task_id)))
job_id <- as.integer(row_value(row, "job_id", sim_id))
out_rds <- file.path(out_dir, sprintf("sim%04d_sim.rds", sim_id))
out_meta <- file.path(out_dir, sprintf("sim%04d_sim_metadata.csv", sim_id))
out_grid_row <- file.path(out_dir, sprintf("sim%04d_grid_row.csv", sim_id))
out_session <- file.path(out_dir, sprintf("sim%04d_sessionInfo.txt", sim_id))
out_log <- file.path(out_dir, sprintf("sim%04d_sim_log.txt", sim_id))

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., sep = ""))
  cat(line, "\n")
  cat(line, "\n", file = out_log, append = TRUE)
}

if (file.exists(out_rds) && !overwrite) {
  log_msg("Simulation cache already exists and overwrite=FALSE: ", out_rds)
  quit(save = "no", status = 0)
}

log_msg("Generating paired SIMMBA cache")
log_msg("Array task ID: ", array_task_id)
log_msg("Task offset: ", task_offset)
log_msg("Effective grid row: ", task_id, " of ", nrow(grid))
log_msg("Simulation ID: ", sim_id)
log_msg("Job ID: ", job_id)
log_msg("Grid file: ", normalizePath(grid_file, mustWork = FALSE))
log_msg("Output RDS: ", normalizePath(out_rds, mustWork = FALSE))
print(row)

load_zentangler_sources(package_dir)
utils::write.csv(row, out_grid_row, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)

TIE <- as_num_or_null(row_value(row, "TIE", 1))
if (is.null(TIE)) TIE <- 1

sim_args <- list(
  nsample = as.integer(row_value(row, "nsample", 600L)),
  nrep = as.integer(row_value(row, "nrep", 10L)),
  seed = as.integer(row_value(row, "sim_seed", row_value(row, "seed", 1L))),
  p.train = 1,
  ygen.mode = as.character(row_value(row, "ygen_mode", "LM")),
  outcome.type = as.character(row_value(row, "outcome_type", "continuous")),
  snr = as.numeric(row_value(row, "snr", 1)),
  n.pathways = as.integer(row_value(row, "n_pathways", 10L)),
  TIE = rep(TIE, 3)
)

log_msg("gen_simmba args: ", paste(names(sim_args), unlist(sim_args), sep = "=", collapse = ", "))
sim <- do.call(gen_simmba, sim_args)
attr(sim, "zentangler_grid_row") <- row
attr(sim, "zentangler_sim_args") <- sim_args
attr(sim, "zentangler_job_id") <- job_id
attr(sim, "zentangler_sim_id") <- sim_id
attr(sim, "zentangler_task_id") <- task_id

saveRDS(sim, out_rds)
meta <- data.frame(
  job_id = job_id,
  sim_id = sim_id,
  task_id = task_id,
  array_task_id = array_task_id,
  task_offset = task_offset,
  sim_path = out_rds,
  nsample = sim_args$nsample,
  nrep = sim_args$nrep,
  seed = sim_args$seed,
  p_train = sim_args$p.train,
  ygen_mode = sim_args$ygen.mode,
  outcome_type = sim_args$outcome.type,
  snr = sim_args$snr,
  n_pathways = sim_args$n.pathways,
  TIE = TIE,
  stringsAsFactors = FALSE
)
utils::write.csv(meta, out_meta, row.names = FALSE)

log_msg("Saved simulation cache: ", out_rds)
log_msg("Saved metadata: ", out_meta)
