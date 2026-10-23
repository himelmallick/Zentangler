#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript make_parallel_nc_grid.R <best_csv> <out_csv>", call. = FALSE)
}

best_file <- args[[1L]]
out_file <- args[[2L]]
grid <- utils::read.csv(best_file, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c(
  "mae", "mae_file", "exposure_mode", "exposure_label",
  "fusion_mode", "sis_n", "glmnet_alpha", "outcome",
  "n_active_mediators", "x_var", "y_family", "covariates",
  "n_samples_after_exposure_subset", "exposure_levels",
  "screen_method", "lambda_choice", "b_inference",
  "fdr_method", "fdr_scope"
)
missing_needed <- setdiff(needed, colnames(grid))
if (length(missing_needed) > 0L) {
  stop("Best-results table is missing columns: ", paste(missing_needed, collapse = ", "), call. = FALSE)
}
# Carry the A-stage settings through when present (grids written after presets were removed).
grid <- grid[, c(needed, intersect(c("a_stage_model", "transform", "a_stage_random_effects"), colnames(grid))), drop = FALSE]
# Legacy B-path inference names all map to the single de-biased estimator.
grid$b_inference[grid$b_inference %in% c("debiased_lasso", "debiased_logistic_lasso", "debiased_cox_lasso", "refit")] <- "debiased"
grid$job_id <- seq_len(nrow(grid))
utils::write.csv(grid, out_file, row.names = FALSE)
message("Wrote parallel NC grid: ", out_file)
message("Rows: ", nrow(grid))
