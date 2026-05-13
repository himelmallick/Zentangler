#!/usr/bin/env Rscript
# Run one row of a Zentangler InterSIM/SIMMBA HPC grid.
#
# Usage locally:
#   Rscript inst/scripts/run_intersim_hpc_task.R \
#     --grid intersim_zentangler_grid.csv \
#     --task-id 1 \
#     --out-dir hpc_results \
#     --package-dir /path/to/Zentangler
#
# Usage with SLURM:
#   Rscript inst/scripts/run_intersim_hpc_task.R \
#     --grid intersim_zentangler_grid.csv \
#     --task-id ${SLURM_ARRAY_TASK_ID} \
#     --task-offset 0 \
#     --out-dir hpc_results \
#     --package-dir /home/fs01/naa4050/Zentangler
#
# If your cluster has MaxArraySize <= 1000, submit repeated arrays with
# --array=1-1000 and increase --task-offset by 1000 for each chunk.

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

as_chr_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.character(x)
}

as_int_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.integer(x)
}

as_num_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.numeric(x)
}

as_num_vec <- function(x, default = 0.25) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  vals <- unlist(strsplit(as.character(x), "[,;|[:space:]]+"))
  vals <- vals[nzchar(vals)]
  out <- suppressWarnings(as.numeric(vals))
  out <- out[is.finite(out)]
  if (length(out) == 0) default else out
}

row_value <- function(row, name, default = NULL) {
  if (!(name %in% colnames(row))) return(default)
  val <- row[[name]][1]
  if (length(val) == 0 || is.na(val) || !nzchar(as.character(val))) return(default)
  val
}

safe_tag <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(as.character(x)), "NA", as.character(x))
  gsub("[^A-Za-z0-9]+", "", x)
}

load_zentangler <- function(package_dir = NULL) {
  if (!is.null(package_dir) && nzchar(package_dir) && file.exists(file.path(package_dir, "DESCRIPTION"))) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(package_dir, quiet = TRUE)
      return(invisible(TRUE))
    }
    r_files <- sort(list.files(file.path(package_dir, "R"), pattern = "[.]R$", full.names = TRUE))
    for (f in r_files) source(f)
    return(invisible(TRUE))
  }

  suppressPackageStartupMessages(library(Zentangler))
  invisible(TRUE)
}

write_if_nonempty <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.data.frame(x) && nrow(x) == 0L) return(invisible(FALSE))
  utils::write.csv(x, file, row.names = FALSE)
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$grid)) stop("Missing --grid")
if (is.null(args$`task-id`)) stop("Missing --task-id")

args$out_dir <- if (!is.null(args$`out-dir`)) args$`out-dir` else if (!is.null(args$out_dir)) args$out_dir else "zentangler_hpc_results"
args$package_dir <- if (!is.null(args$`package-dir`)) args$`package-dir` else if (!is.null(args$package_dir)) args$package_dir else NULL
args$sim_dir <- if (!is.null(args$`sim-dir`)) args$`sim-dir` else if (!is.null(args$sim_dir)) args$sim_dir else NULL

grid <- utils::read.csv(args$grid, stringsAsFactors = FALSE, check.names = FALSE)
task_offset <- if (!is.null(args$`task-offset`)) {
  as.integer(args$`task-offset`)
} else if (!is.null(args$task_offset)) {
  as.integer(args$task_offset)
} else {
  0L
}
if (!is.finite(task_offset) || task_offset < 0L) stop("task-offset must be a non-negative integer.")
array_task_id <- as.integer(args$`task-id`)
task_id <- array_task_id + task_offset
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop(
    "effective task-id must be between 1 and ", nrow(grid),
    ". Got array task-id ", args$`task-id`, " + offset ", task_offset,
    " = ", task_id
  )
}

row <- grid[task_id, , drop = FALSE]
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

method_preset <- as.character(row_value(row, "method_preset", "custom"))
a_stage_model <- as.character(row_value(row, "a_stage_model", "lm"))
screen_method <- as.character(row_value(row, "screen_method", "sis"))
sis_n <- if (identical(screen_method, "none")) NULL else as_int_or_null(row_value(row, "sis_n", NULL))
glmnet_alpha <- as_num_or_null(row_value(row, "glmnet_alpha", 1))
if (is.null(glmnet_alpha)) glmnet_alpha <- 1
fdr_method <- as.character(row_value(row, "fdr_method", "BH"))
q_thresholds <- as_num_vec(row_value(row, "q_thresholds", row_value(row, "q_threshold", 0.25)), default = 0.25)
fdr_scope <- as.character(row_value(row, "fdr_scope", "both"))

maaslin2_random_effect <- as_chr_or_null(row_value(row, "maaslin2_random_effect", NULL))
maaslin2_normalization <- as.character(row_value(row, "maaslin2_normalization", "NONE"))
maaslin2_transform <- as.character(row_value(row, "maaslin2_transform", "NONE"))
maaslin2_analysis_method <- as.character(row_value(row, "maaslin2_analysis_method", "LM"))
maaslin2_standardize <- as_bool(row_value(row, "maaslin2_standardize", FALSE), default = FALSE)

TIE <- as_num_or_null(row_value(row, "TIE", 1))
if (is.null(TIE)) TIE <- 1

tag <- sprintf(
  "job%04d_%s_%s_%s_%s_%s_%s_%s_n%s_sis%s_lam%s_a%s_q%s_scope%s_snr%s_path%s_tie%s",
  as.integer(row_value(row, "job_id", task_id)),
  safe_tag(method_preset),
  safe_tag(a_stage_model),
  safe_tag(row_value(row, "fusion_mode", "early")),
  safe_tag(screen_method),
  safe_tag(row_value(row, "outcome_type", "continuous")),
  safe_tag(row_value(row, "ygen_mode", "LM")),
  safe_tag(row_value(row, "b_inference", "debiased_lasso")),
  safe_tag(row_value(row, "nsample", NA)),
  safe_tag(row_value(row, "sis_n", NA)),
  safe_tag(row_value(row, "lambda_choice", "lambda.1se")),
  gsub("\\.", "p", as.character(glmnet_alpha)),
  paste0(safe_tag(fdr_method), "grid", gsub("\\.", "p", as.character(min(q_thresholds))), "to", gsub("\\.", "p", as.character(max(q_thresholds)))),
  safe_tag(fdr_scope),
  safe_tag(row_value(row, "snr", NA)),
  safe_tag(row_value(row, "n_pathways", NA)),
  gsub("\\.", "p", as.character(row_value(row, "TIE", 1)))
)

out_rds <- file.path(args$out_dir, paste0(tag, ".rds"))
out_summary <- file.path(args$out_dir, paste0(tag, "_summary.csv"))
out_detail <- file.path(args$out_dir, paste0(tag, "_detail.csv"))
out_errors <- file.path(args$out_dir, paste0(tag, "_errors.txt"))
out_settings <- file.path(args$out_dir, paste0(tag, "_settings.csv"))
out_grid_row <- file.path(args$out_dir, paste0(tag, "_grid_row.csv"))
out_model_summary <- file.path(args$out_dir, paste0(tag, "_fit_model_summaries.csv"))
out_view_summary <- file.path(args$out_dir, paste0(tag, "_fit_view_summaries.csv"))
out_session <- file.path(args$out_dir, paste0(tag, "_sessionInfo.txt"))
out_sim_source <- file.path(args$out_dir, paste0(tag, "_simulation_source.csv"))

sim_cache_file <- function(sim_dir, row, task_id) {
  if (is.null(sim_dir) || !nzchar(sim_dir)) return(NULL)
  sim_id <- as.integer(row_value(row, "sim_id", row_value(row, "job_id", task_id)))
  file.path(sim_dir, sprintf("sim%04d_sim.rds", sim_id))
}

maaslin2_output_dir <- as_chr_or_null(row_value(row, "maaslin2_output_dir", NULL))
if (identical(a_stage_model, "maaslin2") && is.null(maaslin2_output_dir)) {
  maaslin2_output_dir <- file.path(args$out_dir, "maaslin2_a_stage", tag)
}

cat("Running Zentangler HPC task\n")
cat("Task:", task_id, "of", nrow(grid), "\n")
print(row)
cat("SIM_DIR:", ifelse(is.null(args$sim_dir), "<none; generate internally>", args$sim_dir), "\n")

load_zentangler(args$package_dir)
utils::write.csv(row, out_grid_row, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)

sim_obj <- NULL
sim_path <- sim_cache_file(args$sim_dir, row, task_id)
if (!is.null(sim_path)) {
  if (!file.exists(sim_path)) {
    stop("Requested simulation cache does not exist: ", sim_path)
  }
  cat("Loading paired simulation cache:", sim_path, "\n")
  sim_obj <- readRDS(sim_path)
} else {
  cat("No simulation cache supplied; Zentangler will generate data internally.\n")
}
utils::write.csv(
  data.frame(
    task_id = task_id,
    job_id = as.integer(row_value(row, "job_id", task_id)),
    sim_id = as.integer(row_value(row, "sim_id", row_value(row, "job_id", task_id))),
    sim_path = ifelse(is.null(sim_path), NA_character_, sim_path),
    used_cache = !is.null(sim_obj),
    stringsAsFactors = FALSE
  ),
  out_sim_source,
  row.names = FALSE
)

start_time <- Sys.time()
res <- try(
  run_intersim_zentangler(
    nrep = as.integer(row_value(row, "nrep", 10L)),
    nsample = as.integer(row_value(row, "nsample", 600L)),
    fusion_modes = as.character(row_value(row, "fusion_mode", "early")),
    method_preset = method_preset,
    sis_n = sis_n,
    sis_rank = as.character(row_value(row, "sis_rank", "abs_a")),
    screen_method = screen_method,
    a_stage_model = a_stage_model,
    maaslin2_random_effect = maaslin2_random_effect,
    maaslin2_normalization = maaslin2_normalization,
    maaslin2_transform = maaslin2_transform,
    maaslin2_analysis_method = maaslin2_analysis_method,
    maaslin2_standardize = maaslin2_standardize,
    maaslin2_output_dir = maaslin2_output_dir,
    q_threshold = q_thresholds,
    fdr_method = fdr_method,
    fdr_scope = fdr_scope,
    top_n = as.integer(row_value(row, "top_n", 50L)),
    seed = as.integer(row_value(row, "sim_seed", row_value(row, "seed", 1L))),
    p.train = 1,
    ygen.mode = as.character(row_value(row, "ygen_mode", "LM")),
    outcome.type = as.character(row_value(row, "outcome_type", "continuous")),
    lambda_choice = as.character(row_value(row, "lambda_choice", "lambda.1se")),
    glmnet_alpha = glmnet_alpha,
    b_inference = as.character(row_value(row, "b_inference", "debiased_lasso")),
    debias_max_targets = as.integer(row_value(row, "debias_max_targets", 200L)),
    coop_rho = as.numeric(row_value(row, "coop_rho", 0.20)),
    residualize = as_bool(row_value(row, "residualize", FALSE), default = FALSE),
    bootstrap_repeats = as.integer(row_value(row, "bootstrap_repeats", 0L)),
    snr = as.numeric(row_value(row, "snr", 1)),
    n.pathways = as.integer(row_value(row, "n_pathways", 10L)),
    TIE = rep(TIE, 3),
    sim = sim_obj,
    verbose = TRUE
  ),
  silent = TRUE
)
end_time <- Sys.time()

if (inherits(res, "try-error")) {
  msg <- conditionMessage(attr(res, "condition"))
  writeLines(msg, out_errors)
  stop(msg)
}

res$hpc <- list(
  task_id = task_id,
  array_task_id = array_task_id,
  task_offset = task_offset,
  grid_row = row,
  started = start_time,
  finished = end_time,
  runtime_minutes = as.numeric(difftime(end_time, start_time, units = "mins")),
  package_dir = args$package_dir,
  maaslin2_output_dir = maaslin2_output_dir
)

fit_model_summaries <- if (length(res$fits) > 0L) {
  do.call(rbind, lapply(names(res$fits), function(nm) {
    x <- res$fits[[nm]]$model_summary
    if (is.null(x) || nrow(x) == 0L) return(NULL)
    data.frame(fit_id = nm, x, stringsAsFactors = FALSE)
  }))
} else {
  NULL
}

fit_view_summaries <- if (length(res$fits) > 0L) {
  do.call(rbind, lapply(names(res$fits), function(nm) {
    x <- res$fits[[nm]]$view_summary
    if (is.null(x) || nrow(x) == 0L) return(NULL)
    data.frame(fit_id = nm, x, stringsAsFactors = FALSE)
  }))
} else {
  NULL
}

saveRDS(res, out_rds)
utils::write.csv(res$summary, out_summary, row.names = FALSE)
utils::write.csv(res$detail, out_detail, row.names = FALSE)
write_if_nonempty(fit_model_summaries, out_model_summary)
write_if_nonempty(fit_view_summaries, out_view_summary)

settings_df <- data.frame(
  setting = names(res$settings),
  value = vapply(res$settings, function(x) paste(as.character(x), collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(settings_df, out_settings, row.names = FALSE)
if (length(res$errors) > 0) writeLines(res$errors, out_errors) else writeLines(character(0), out_errors)

cat("Finished task", task_id, "\n")
cat("Runtime minutes:", res$hpc$runtime_minutes, "\n")
cat("Saved:", out_rds, "\n")
print(res$summary)
