#!/usr/bin/env Rscript

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

row_value <- function(row, name, default = NULL) {
  if (!(name %in% colnames(row))) return(default)
  val <- row[[name]][1]
  if (length(val) == 0L || is.na(val) || (!is.numeric(val) && !nzchar(as.character(val)))) return(default)
  val
}

safe_tag <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(as.character(x)), "NA", as.character(x))
  gsub("[^A-Za-z0-9]+", "", x)
}

load_mae <- function(mae_file) {
  if (!file.exists(mae_file)) stop("MAE file not found: ", mae_file, call. = FALSE)
  env <- new.env(parent = emptyenv())
  load(mae_file, envir = env)
  if (!exists("mae", envir = env, inherits = FALSE)) {
    stop("Expected object 'mae' in ", mae_file, call. = FALSE)
  }
  get("mae", envir = env, inherits = FALSE)
}

require_coldata_columns <- function(cd, needed, exposure_mode) {
  missing_needed <- setdiff(needed, colnames(cd))
  if (length(missing_needed) > 0L) {
    stop(
      "Exposure mode '", exposure_mode, "' requires colData columns: ",
      paste(missing_needed, collapse = ", "),
      call. = FALSE
    )
  }
}

prepare_mae_exposure <- function(mae, exposure_mode = "existing", x_var = "treatment_any") {
  exposure_mode <- as.character(exposure_mode[[1L]])
  cd <- MultiAssayExperiment::colData(mae)
  keep <- rep(TRUE, nrow(cd))
  derived_x_var <- x_var

  if (identical(exposure_mode, "existing")) {
    if (!(x_var %in% colnames(cd))) stop("x_var '", x_var, "' not found in colData.", call. = FALSE)
  } else if (identical(exposure_mode, "blueberry_any")) {
    require_coldata_columns(cd, "blueberry_dose_g", exposure_mode)
    derived_x_var <- ".__zentangler_exposure__"
    cd[[derived_x_var]] <- as.numeric(cd$blueberry_dose_g > 0)
  } else if (identical(exposure_mode, "blueberry_dose")) {
    require_coldata_columns(cd, "blueberry_dose_g", exposure_mode)
    derived_x_var <- ".__zentangler_exposure__"
    cd[[derived_x_var]] <- as.numeric(cd$blueberry_dose_g)
  } else if (identical(exposure_mode, "Y_vs_Z")) {
    require_coldata_columns(cd, "group", exposure_mode)
    keep <- as.character(cd$group) %in% c("Y", "Z")
    derived_x_var <- ".__zentangler_exposure__"
    cd[[derived_x_var]] <- ifelse(as.character(cd$group) == "Y", 1, ifelse(as.character(cd$group) == "Z", 0, NA_real_))
  } else if (identical(exposure_mode, "Y_vs_X")) {
    require_coldata_columns(cd, "group", exposure_mode)
    keep <- as.character(cd$group) %in% c("Y", "X")
    derived_x_var <- ".__zentangler_exposure__"
    cd[[derived_x_var]] <- ifelse(as.character(cd$group) == "Y", 1, ifelse(as.character(cd$group) == "X", 0, NA_real_))
  } else {
    stop(
      "Unsupported exposure mode: ", exposure_mode, ". ",
      "Use one of: existing, blueberry_any, blueberry_dose, Y_vs_Z, Y_vs_X",
      call. = FALSE
    )
  }

  MultiAssayExperiment::colData(mae) <- cd
  list(
    mae = mae[, keep],
    x_var = derived_x_var
  )
}

write_if_nonempty <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.data.frame(x) && nrow(x) == 0L) return(invisible(FALSE))
  utils::write.csv(x, file, row.names = FALSE)
  invisible(TRUE)
}

split_arg <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(paste(x, collapse = ""))) return(character(0))
  out <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

# Per-view specs in grid cells: "lm" (all views) or "species=tweedieverse;metabolites=lm".
parse_view_spec <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(trimws(as.character(x)))) return(default)
  parts <- trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 1L && !grepl("=", parts)) return(parts)
  kv <- strsplit(parts, "=", fixed = TRUE)
  stats::setNames(vapply(kv, function(p) trimws(p[[2]]), character(1)), vapply(kv, function(p) trimws(p[[1]]), character(1)))
}

# A-stage settings from the grid row: model per view, transform per view and
# random effects shared by the MaAsLin2 / Tweedieverse views.
a_stage_settings <- function(row) {
  model <- parse_view_spec(row_value(row, "a_stage_model", "lm"), "lm")
  re <- split_arg(row_value(row, "a_stage_random_effects", ""))
  method_args <- if (length(re)) list(random_effects = re) else list()
  list(
    a_stage_model = model,
    transform = parse_view_spec(row_value(row, "transform", "none"), "none"),
    maaslin2_args = if (any(model == "maaslin2")) method_args else list(),
    tweedieverse_args = if (any(model == "tweedieverse")) method_args else list(),
    label = as.character(row_value(row, "a_stage_model", "lm"))
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$grid)) stop("Missing --grid", call. = FALSE)
if (is.null(args$`task-id`)) stop("Missing --task-id", call. = FALSE)

out_dir <- if (!is.null(args$`out-dir`)) args$`out-dir` else "parallel_nc_results"
task_offset <- if (!is.null(args$`task-offset`)) as.integer(args$`task-offset`) else 0L
n_perm <- if (!is.null(args$`n-perm`)) as.integer(args$`n-perm`) else 100L
q_threshold <- if (!is.null(args$`q-threshold`)) as.numeric(args$`q-threshold`) else 0.25
seed <- if (!is.null(args$seed)) as.integer(args$seed) else 1L

grid <- utils::read.csv(args$grid, stringsAsFactors = FALSE, check.names = FALSE)
array_task_id <- as.integer(args$`task-id`)
task_id <- array_task_id + task_offset
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop("effective task-id must be between 1 and ", nrow(grid), ". Got ", task_id, call. = FALSE)
}

row <- grid[task_id, , drop = FALSE]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Map a fusion label to the zentangler outcome-model arguments:
#   "multiview" (default; cooperative learning, rho in [0, 1] tuned by CV)
#   "early" / "late" (supplementary IntegratedLearner-style fusion)
# The legacy grid value "intermediate" is read as "multiview".
fusion_args <- function(label) {
  label <- as.character(label)[1]
  if (label %in% c("multiview", "intermediate", "cooperative")) {
    return(list(fusion = "multiview", integrated_learner_fusion = "late"))
  }
  if (label %in% c("early", "late")) {
    return(list(fusion = "integrated_learner", integrated_learner_fusion = label))
  }
  stop("Unknown fusion label: ", label, call. = FALSE)
}

suppressPackageStartupMessages({
  library(Zentangler)
  library(MultiAssayExperiment)
})

mae_name <- as.character(row_value(row, "mae"))
mae_file <- as.character(row_value(row, "mae_file"))
exposure_mode <- as.character(row_value(row, "exposure_mode", "existing"))
# The grid keeps its fusion_mode column (early / late / multiview) so existing
# CSVs still work; it is translated to fusion / integrated_learner_fusion here.
fusion_mode <- as.character(row_value(row, "fusion_mode"))
fusion_settings <- fusion_args(fusion_mode)
sis_n <- as.integer(row_value(row, "sis_n"))
glmnet_alpha <- as.numeric(row_value(row, "glmnet_alpha"))
outcome <- as.character(row_value(row, "outcome"))
n_active_mediators <- as.integer(row_value(row, "n_active_mediators", NA))
x_var <- as.character(row_value(row, "x_var", "treatment_any"))
y_family <- as.character(row_value(row, "y_family", "gaussian"))
covars <- split_arg(row_value(row, "covariates", ""))
a_stage <- a_stage_settings(row)
screen_method <- as.character(row_value(row, "screen_method", "sis"))
lambda_choice <- as.character(row_value(row, "lambda_choice", "lambda.min"))
b_inference <- as.character(row_value(row, "b_inference", "debiased"))
fdr_method <- as.character(row_value(row, "fdr_method", "BH"))
fdr_scope <- as.character(row_value(row, "fdr_scope", "global"))

tag <- sprintf(
  "job%04d_%s_%s_%s_sis%s_alpha%s_%s_perm%s",
  as.integer(row_value(row, "job_id", task_id)),
  safe_tag(mae_name),
  safe_tag(exposure_mode),
  safe_tag(fusion_mode),
  safe_tag(sis_n),
  gsub("\\.", "p", as.character(glmnet_alpha)),
  safe_tag(outcome),
  safe_tag(n_perm)
)

out_rds <- file.path(out_dir, paste0(tag, ".rds"))
out_summary <- file.path(out_dir, paste0(tag, "_empirical_summary.csv"))
out_null <- file.path(out_dir, paste0(tag, "_null_summaries.csv"))
out_observed <- file.path(out_dir, paste0(tag, "_observed_summary.csv"))
out_mediators <- file.path(out_dir, paste0(tag, "_observed_mediators_active.csv"))
out_empirical_mediators <- file.path(out_dir, paste0(tag, "_empirical_mediator_summary.csv"))
out_grid_row <- file.path(out_dir, paste0(tag, "_grid_row.csv"))
out_session <- file.path(out_dir, paste0(tag, "_sessionInfo.txt"))

cat("Running generic parallel negative-control task\n")
cat("Task:", task_id, "of", nrow(grid), "\n")
print(row)

mae_obj0 <- load_mae(mae_file)
prepared <- prepare_mae_exposure(mae_obj0, exposure_mode = exposure_mode, x_var = x_var)
mae_obj <- prepared$mae
utils::write.csv(row, out_grid_row, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)

fit <- run_negative_controls(
  mae = mae_obj,
  x_var = prepared$x_var,
  y_var = outcome,
  n_perm = n_perm,
  q_threshold = q_threshold,
  seed = seed,
  y_family = y_family,
  covariates = covars,
  transform = a_stage$transform,
  a_stage_model = a_stage$a_stage_model,
  maaslin2_args = a_stage$maaslin2_args,
  tweedieverse_args = a_stage$tweedieverse_args,
  screen_method = screen_method,
  sis_n = sis_n,
  glmnet_alpha = glmnet_alpha,
  lambda_choice = lambda_choice,
  fusion = fusion_settings$fusion,
  integrated_learner_fusion = fusion_settings$integrated_learner_fusion,
  b_inference = b_inference,
  fdr_method = fdr_method,
  fdr_scope = fdr_scope
)

saveRDS(fit, out_rds)
write_if_nonempty(fit$empirical_summary, out_summary)
write_if_nonempty(fit$null_summaries, out_null)
write_if_nonempty(fit$observed_summary, out_observed)
write_if_nonempty(zentangler_mediators(fit$observed_fit, q_threshold = q_threshold), out_mediators)
write_if_nonempty(fit$empirical_units, out_empirical_mediators)
# Per-permutation null mediator tables are no longer returned; per-mediator null
# frequencies and p-values are in fit$empirical_units (written above).

cat("Saved RDS:", out_rds, "\n")
cat("Observed active mediators:", n_active_mediators, "\n")
cat("Done.\n")
