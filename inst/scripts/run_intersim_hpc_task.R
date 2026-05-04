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
#     --out-dir hpc_results

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

as_bool <- function(x) {
  if (length(x) == 0 || is.na(x)) return(FALSE)
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

load_zentangler <- function(package_dir = NULL) {
  if (!is.null(package_dir) && nzchar(package_dir) && file.exists(file.path(package_dir, "DESCRIPTION"))) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(package_dir, quiet = TRUE)
      return(invisible(TRUE))
    }
    source(file.path(package_dir, "R", "zentangler.R"))
    source(file.path(package_dir, "R", "trigger_intersim.R"))
    source(file.path(package_dir, "R", "simulation.R"))
    source(file.path(package_dir, "R", "run_intersim_zentangler.R"))
    return(invisible(TRUE))
  }

  suppressPackageStartupMessages(library(Zentangler))
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$grid)) stop("Missing --grid")
if (is.null(args$`task-id`)) stop("Missing --task-id")

args$out_dir <- if (!is.null(args$`out-dir`)) args$`out-dir` else if (!is.null(args$out_dir)) args$out_dir else "zentangler_hpc_results"
args$package_dir <- if (!is.null(args$`package-dir`)) args$`package-dir` else if (!is.null(args$package_dir)) args$package_dir else NULL

grid <- utils::read.csv(args$grid, stringsAsFactors = FALSE, check.names = FALSE)
task_id <- as.integer(args$`task-id`)
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop("task-id must be between 1 and ", nrow(grid), ". Got: ", args$`task-id`)
}

row <- grid[task_id, , drop = FALSE]
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
glmnet_alpha <- if ("glmnet_alpha" %in% colnames(row)) as_num_or_null(row$glmnet_alpha) else 1
if (is.null(glmnet_alpha)) glmnet_alpha <- 1

tag <- sprintf(
  "job%04d_%s_%s_n%s_sis%s_lam%s_a%s_q%s_snr%s_path%s",
  as.integer(row$job_id),
  as.character(row$fusion_mode),
  as.character(row$screen_method),
  as.character(row$nsample),
  ifelse(is.na(row$sis_n) || !nzchar(as.character(row$sis_n)), "NA", as.character(row$sis_n)),
  gsub("[^A-Za-z0-9]+", "", as.character(row$lambda_choice)),
  gsub("\\.", "p", as.character(glmnet_alpha)),
  gsub("\\.", "p", as.character(row$q_threshold)),
  as.character(row$snr),
  as.character(row$n_pathways)
)

out_rds <- file.path(args$out_dir, paste0(tag, ".rds"))
out_summary <- file.path(args$out_dir, paste0(tag, "_summary.csv"))
out_detail <- file.path(args$out_dir, paste0(tag, "_detail.csv"))
out_errors <- file.path(args$out_dir, paste0(tag, "_errors.txt"))
out_settings <- file.path(args$out_dir, paste0(tag, "_settings.csv"))

cat("Running Zentangler HPC task\n")
cat("Task:", task_id, "of", nrow(grid), "\n")
print(row)

load_zentangler(args$package_dir)

screen_method <- as.character(row$screen_method)
sis_n <- if (identical(screen_method, "none")) NULL else as_int_or_null(row$sis_n)
a_stage_model <- if ("a_stage_model" %in% colnames(row) && nzchar(as.character(row$a_stage_model))) {
  as.character(row$a_stage_model)
} else {
  "lm"
}
maaslin2_random_effect <- if ("maaslin2_random_effect" %in% colnames(row)) as_chr_or_null(row$maaslin2_random_effect) else NULL
maaslin2_normalization <- if ("maaslin2_normalization" %in% colnames(row) && nzchar(as.character(row$maaslin2_normalization))) as.character(row$maaslin2_normalization) else "NONE"
maaslin2_transform <- if ("maaslin2_transform" %in% colnames(row) && nzchar(as.character(row$maaslin2_transform))) as.character(row$maaslin2_transform) else "NONE"
maaslin2_analysis_method <- if ("maaslin2_analysis_method" %in% colnames(row) && nzchar(as.character(row$maaslin2_analysis_method))) as.character(row$maaslin2_analysis_method) else "LM"
maaslin2_standardize <- if ("maaslin2_standardize" %in% colnames(row)) as_bool(row$maaslin2_standardize) else FALSE
TIE <- as_num_or_null(row$TIE)
if (is.null(TIE)) TIE <- 1

start_time <- Sys.time()
res <- try(
  run_intersim_zentangler(
    nrep = as.integer(row$nrep),
    nsample = as.integer(row$nsample),
    fusion_modes = as.character(row$fusion_mode),
    sis_n = sis_n,
    sis_rank = as.character(row$sis_rank),
    screen_method = screen_method,
    a_stage_model = a_stage_model,
    maaslin2_random_effect = maaslin2_random_effect,
    maaslin2_normalization = maaslin2_normalization,
    maaslin2_transform = maaslin2_transform,
    maaslin2_analysis_method = maaslin2_analysis_method,
    maaslin2_standardize = maaslin2_standardize,
    q_threshold = as.numeric(row$q_threshold),
    top_n = as.integer(row$top_n),
    seed = as.integer(row$seed),
    p.train = as.numeric(row$p_train),
    ygen.mode = as.character(row$ygen_mode),
    outcome.type = as.character(row$outcome_type),
    lambda_choice = as.character(row$lambda_choice),
    glmnet_alpha = glmnet_alpha,
    b_inference = as.character(row$b_inference),
    debias_max_targets = as.integer(row$debias_max_targets),
    coop_rho = as.numeric(row$coop_rho),
    residualize = as_bool(row$residualize),
    bootstrap_repeats = as.integer(row$bootstrap_repeats),
    snr = as.numeric(row$snr),
    n.pathways = as.integer(row$n_pathways),
    TIE = rep(TIE, 3),
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
  grid_row = row,
  started = start_time,
  finished = end_time,
  runtime_minutes = as.numeric(difftime(end_time, start_time, units = "mins"))
)

saveRDS(res, out_rds)
utils::write.csv(res$summary, out_summary, row.names = FALSE)
utils::write.csv(res$detail, out_detail, row.names = FALSE)
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
