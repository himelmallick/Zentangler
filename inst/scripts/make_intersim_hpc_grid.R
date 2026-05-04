#!/usr/bin/env Rscript
# Create a Zentangler InterSIM/SIMMBA HPC parameter grid.
#
# Usage:
#   Rscript inst/scripts/make_intersim_hpc_grid.R --out intersim_grid.csv
#
# Edit the vectors below to expand or shrink the simulation study.

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

# Core simulation size grid.
nsample_values <- c(300L, 600L, 1000L)
nrep_per_job <- 10L

# Method grid.
fusion_modes <- c("early", "intermediate", "late")
lambda_choices <- c("lambda.1se", "lambda.min")
# glmnet_alpha controls the Y-stage penalty for early/late fusion:
# 1 = lasso, 0.5/0.75 = elastic net, 0 = ridge.
# Intermediate cooperative fusion currently uses alpha = 1, so it is not
# repeated over the elastic-net alpha grid.
glmnet_alpha_values <- c(1, 0.75, 0.5)
q_thresholds <- c(0.10, 0.20, 0.25, 0.30)

# Screening grid.
# For screen_method = "none", sis_n is ignored and stored as NA.
sis_n_values <- c(50L, 100L, 200L, 300L)
sis_rank_values <- c("abs_a")
screen_sis <- expand.grid(
  screen_method = "sis",
  sis_n = sis_n_values,
  sis_rank = sis_rank_values,
  stringsAsFactors = FALSE
)
screen_none <- data.frame(
  screen_method = "none",
  sis_n = NA_integer_,
  sis_rank = "abs_a",
  stringsAsFactors = FALSE
)
screen_grid <- rbind(screen_sis, screen_none)

# Data-generating grid.
# TIE is repeated across the three InterSIM views inside the task runner.
signal_grid <- expand.grid(
  snr = c(1, 2),
  n_pathways = c(10L, 30L),
  TIE = c(1),
  stringsAsFactors = FALSE
)

early_late_grid <- expand.grid(
  nsample = nsample_values,
  fusion_mode = c("early", "late"),
  lambda_choice = lambda_choices,
  glmnet_alpha = glmnet_alpha_values,
  q_threshold = q_thresholds,
  outcome_type = "continuous",
  ygen_mode = "LM",
  b_inference = "debiased_lasso",
  coop_rho = 0.20,
  stringsAsFactors = FALSE
)
intermediate_grid <- expand.grid(
  nsample = nsample_values,
  fusion_mode = "intermediate",
  lambda_choice = lambda_choices,
  glmnet_alpha = 1,
  q_threshold = q_thresholds,
  outcome_type = "continuous",
  ygen_mode = "LM",
  b_inference = "debiased_lasso",
  coop_rho = 0.20,
  stringsAsFactors = FALSE
)
base_grid <- rbind(early_late_grid, intermediate_grid)

grid <- merge(base_grid, screen_grid, all = TRUE)
grid <- merge(grid, signal_grid, all = TRUE)

grid$nrep <- nrep_per_job
grid$top_n <- 50L
grid$p_train <- 0.70
grid$debias_max_targets <- 200L
grid$bootstrap_repeats <- 0L
grid$residualize <- FALSE

# Give every row a reproducible seed and stable ID.
grid <- grid[order(
  grid$nsample,
  grid$screen_method,
  grid$sis_n,
  grid$fusion_mode,
  grid$lambda_choice,
  grid$glmnet_alpha,
  grid$q_threshold,
  grid$snr,
  grid$n_pathways
), , drop = FALSE]
grid$job_id <- seq_len(nrow(grid))
grid$seed <- 100000L + grid$job_id * 100L

grid <- grid[, c(
  "job_id", "seed", "nrep", "nsample", "p_train",
  "fusion_mode", "screen_method", "sis_n", "sis_rank",
  "lambda_choice", "glmnet_alpha", "q_threshold", "top_n",
  "outcome_type", "ygen_mode", "b_inference", "debias_max_targets",
  "coop_rho", "bootstrap_repeats", "residualize",
  "snr", "n_pathways", "TIE"
)]

utils::write.csv(grid, out_file, row.names = FALSE, na = "")
cat("Wrote grid:", normalizePath(out_file, mustWork = FALSE), "\n")
cat("Number of jobs:", nrow(grid), "\n")
