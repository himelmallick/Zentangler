#!/usr/bin/env Rscript
# Create a Zentangler InterSIM/SIMMBA HPC parameter grid.
#
# This version is synced to the current Zentangler API and is designed for a
# realistic HPC simulation campaign rather than an enormous full Cartesian grid.
#
# Usage:
#   Rscript inst/scripts/make_intersim_hpc_grid.R \
#     --out intersim_zentangler_grid.csv
#
# Optional profiles:
#   --profile focused   small smoke/iteration grid, about 320 jobs
#   --profile hpc5000   default, continuous-outcome method grid
#   --profile expanded  alias for hpc5000 to avoid accidental huge grids
#   --profile full      larger follow-up grid, includes nonlinear continuous settings
#
# The companion script run_intersim_hpc_task.R reads one row and calls
# run_intersim_zentangler().

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

args <- parse_args(commandArgs(trailingOnly = TRUE))
out_file <- if (!is.null(args$out)) args$out else "intersim_zentangler_grid.csv"
profile <- if (!is.null(args$profile)) as.character(args$profile) else "hpc5000"
profile <- match.arg(profile, choices = c("focused", "hpc5000", "expanded", "full"))
if (identical(profile, "expanded")) profile <- "hpc5000"

# -------------------------------------------------------------------------
# Grid profiles
# -------------------------------------------------------------------------
# focused:
#   Fast smoke grid for checking package/HPC setup.
# hpc5000:
#   Recommended first serious HPC grid. It covers the important method axes
#   for continuous outcomes, where current de-biased lasso B-path inference is valid.
# full:
#   Follow-up grid after hpc5000 identifies promising regions.
# -------------------------------------------------------------------------

if (identical(profile, "focused")) {
  nsample_values <- c(300L, 600L)
  nrep_per_job <- 5L
  a_stage_models <- c("lm", "maaslin2")
  lambda_choices <- c("lambda.1se", "lambda.min")
  glmnet_alpha_values <- c(1, 0.5)
  fdr_methods <- c("BH", "BY")
  q_thresholds <- seq(0.05, 0.25, by = 0.05)
  sis_n_values <- c(100L)
  sis_rank_values <- c("abs_a")
  include_screen_none <- TRUE
  outcome_grid <- expand.grid(
    outcome_type = "continuous",
    ygen_mode = "LM",
    stringsAsFactors = FALSE
  )
  signal_grid <- expand.grid(
    snr = c(1, 2),
    n_pathways = c(10L, 30L),
    TIE = c(1),
    stringsAsFactors = FALSE
  )
  b_inference_values <- c("debiased_lasso")
} else if (identical(profile, "hpc5000")) {
  # Count logic with current default settings:
  # base model rows = 200. Five q cutoffs are evaluated inside each fitted run.
  # screening rows = 9: SIS 100/200/250/500 x abs_a/pvalue + none
  # outcome rows = 1: continuous LM only
  # signal rows = 4: snr 1/2 x pathways 10/30, TIE fixed at 1
  # total = 200 * 9 * 1 * 4 = 7200 jobs
  nsample_values <- c(500L, 1000L, 1500L, 2000L, 3000L)
  nrep_per_job <- 10L
  a_stage_models <- c("lm", "maaslin2")
  lambda_choices <- c("lambda.1se", "lambda.min")
  glmnet_alpha_values <- c(1, 0.5) # 1 = lasso, 0.5 = elastic net
  fdr_methods <- c("BH", "BY")
  q_thresholds <- seq(0.05, 0.25, by = 0.05)
  sis_n_values <- c(100L, 200L, 250L, 500L)
  sis_rank_values <- c("abs_a", "pvalue")
  include_screen_none <- TRUE
  outcome_grid <- expand.grid(
    outcome_type = "continuous",
    ygen_mode = "LM",
    stringsAsFactors = FALSE
  )
  signal_grid <- expand.grid(
    snr = c(1, 2),
    n_pathways = c(10L, 30L),
    TIE = c(1),
    stringsAsFactors = FALSE
  )
  b_inference_values <- c("debiased_lasso")
} else {
  # Larger follow-up grid, not recommended as the first run.
  # Continuous-only follow-up with nonlinear Y-generation modes.
  nsample_values <- c(300L, 600L, 1000L)
  nrep_per_job <- 10L
  a_stage_models <- c("lm", "maaslin2")
  lambda_choices <- c("lambda.1se", "lambda.min")
  glmnet_alpha_values <- c(1, 0.75, 0.5)
  fdr_methods <- c("BH", "BY")
  q_thresholds <- seq(0.05, 0.25, by = 0.05)
  sis_n_values <- c(50L, 100L, 200L)
  sis_rank_values <- c("abs_a", "pvalue")
  include_screen_none <- TRUE
  outcome_grid <- expand.grid(
    outcome_type = "continuous",
    ygen_mode = c("LM", "Friedman"),
    stringsAsFactors = FALSE
  )
  signal_grid <- expand.grid(
    snr = c(1, 2),
    n_pathways = c(10L, 30L),
    TIE = c(1),
    stringsAsFactors = FALSE
  )
  b_inference_values <- c("debiased_lasso")
}

# Screening grid.
screen_sis <- expand.grid(
  screen_method = "sis",
  sis_n = sis_n_values,
  sis_rank = sis_rank_values,
  stringsAsFactors = FALSE
)
if (isTRUE(include_screen_none)) {
  screen_none <- data.frame(
    screen_method = "none",
    sis_n = NA_integer_,
    sis_rank = "abs_a",
    stringsAsFactors = FALSE
  )
  screen_grid <- rbind(screen_sis, screen_none)
} else {
  screen_grid <- screen_sis
}

# Early/late fusion can use lasso and elastic net via glmnet_alpha.
early_late_grid <- expand.grid(
  nsample = nsample_values,
  method_preset = "custom",
  a_stage_model = a_stage_models,
  fusion_mode = c("early", "late"),
  lambda_choice = lambda_choices,
  glmnet_alpha = glmnet_alpha_values,
  fdr_method = fdr_methods,
  b_inference = b_inference_values,
  stringsAsFactors = FALSE
)

# Intermediate fusion is cooperative/lasso-style in the current package.
intermediate_grid <- expand.grid(
  nsample = nsample_values,
  method_preset = "custom",
  a_stage_model = a_stage_models,
  fusion_mode = "intermediate",
  lambda_choice = lambda_choices,
  glmnet_alpha = 1,
  fdr_method = fdr_methods,
  b_inference = b_inference_values,
  stringsAsFactors = FALSE
)

base_grid <- rbind(early_late_grid, intermediate_grid)
grid <- merge(base_grid, screen_grid, all = TRUE)
grid <- merge(grid, outcome_grid, all = TRUE)
grid <- merge(grid, signal_grid, all = TRUE)

# MaAsLin2 A-stage controls. These columns are passed directly to the task
# runner. InterSIM/SIMMBA does not create repeated-subject IDs by default, so
# maaslin2_random_effect is empty unless you edit it for a longitudinal design.
grid$maaslin2_random_effect <- NA_character_
grid$maaslin2_normalization <- "NONE"
grid$maaslin2_transform <- "NONE"
grid$maaslin2_analysis_method <- "LM"
grid$maaslin2_standardize <- FALSE

# Shared run controls.
grid$grid_profile <- profile
grid$nrep <- nrep_per_job
grid$top_n <- 50L
grid$q_thresholds <- paste(format(q_thresholds, trim = TRUE, scientific = FALSE), collapse = ",")
grid$fdr_scope <- "both"
grid$p_train <- 0.70
grid$debias_max_targets <- 200L
grid$coop_rho <- 0.20
grid$bootstrap_repeats <- 0L
grid$residualize <- FALSE

# Give every row a reproducible seed and stable ID.
grid <- grid[order(
  grid$nsample,
  grid$outcome_type,
  grid$ygen_mode,
  grid$a_stage_model,
  grid$screen_method,
  grid$sis_rank,
  grid$sis_n,
  grid$fusion_mode,
  grid$lambda_choice,
  grid$glmnet_alpha,
  grid$b_inference,
  grid$fdr_method,
  grid$snr,
  grid$n_pathways,
  grid$TIE
), , drop = FALSE]
grid$job_id <- seq_len(nrow(grid))
grid$seed <- 100000L + grid$job_id * 100L

grid <- grid[, c(
  "job_id", "seed", "grid_profile", "nrep", "nsample", "p_train",
  "method_preset", "a_stage_model", "fusion_mode",
  "screen_method", "sis_n", "sis_rank",
  "lambda_choice", "glmnet_alpha", "fdr_method", "q_thresholds", "fdr_scope", "top_n",
  "outcome_type", "ygen_mode", "b_inference", "debias_max_targets",
  "coop_rho", "bootstrap_repeats", "residualize",
  "maaslin2_random_effect", "maaslin2_normalization", "maaslin2_transform",
  "maaslin2_analysis_method", "maaslin2_standardize",
  "snr", "n_pathways", "TIE"
)]

utils::write.csv(grid, out_file, row.names = FALSE, na = "")
cat("Wrote grid:", normalizePath(out_file, mustWork = FALSE), "\n")
cat("Profile:", profile, "\n")
cat("Number of jobs:", nrow(grid), "\n")
cat("A-stage models:", paste(unique(grid$a_stage_model), collapse = ", "), "\n")
cat("Fusion modes:", paste(unique(grid$fusion_mode), collapse = ", "), "\n")
cat("Screening:", paste(unique(grid$screen_method), collapse = ", "), "\n")
cat("SIS ranks:", paste(unique(grid$sis_rank), collapse = ", "), "\n")
cat("glmnet alpha values:", paste(sort(unique(grid$glmnet_alpha)), collapse = ", "), "\n")
cat("FDR methods:", paste(unique(grid$fdr_method), collapse = ", "), "\n")
cat("q thresholds per fit:", paste(unique(grid$q_thresholds), collapse = "; "), "\n")
cat("FDR scopes per fit:", paste(unique(grid$fdr_scope), collapse = ", "), "\n")
cat("Outcome types:", paste(unique(grid$outcome_type), collapse = ", "), "\n")
cat("Y generation modes:", paste(unique(grid$ygen_mode), collapse = ", "), "\n")
cat("Signal-to-noise values:", paste(sort(unique(grid$snr)), collapse = ", "), "\n")
cat("Pathway counts:", paste(sort(unique(grid$n_pathways)), collapse = ", "), "\n")
