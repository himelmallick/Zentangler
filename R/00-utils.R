# Zentangler MAE-first K-View Parallel Mediation Core
#
# Active API:
#   fit_multiview_parallel_zentangler(mae, x_var, y_var, ...)
#
# Data contract:
# - mae is a MultiAssayExperiment.
# - each experiment is one modality/view, usually stored as features x samples.
# - colData(mae) contains phenotype variables, including x_var and y_var.
# - the internal engine still uses sample x feature matrices after extraction.
# - x_var and y_var are numeric columns in colData(mae).
#
# Method structure:
# - A-stage: per-view HIMA-like univariate screen, M_j^(k) ~ X + C,
#   or optional MaAsLin2 fixed/random-effect screen for longitudinal data.
# - SIS: keep top mediators per view.
# - B/Y-stage: fit selected mediators jointly using early, intermediate, or late
#   fusion. Early/late use glmnet and can be lasso or elastic net through
#   glmnet_alpha; intermediate currently uses a cooperative lasso-style update.
# - Score: a_j * b_j, reported for every mediator with view labels preserved.
#
# Scientific status:
# - This is a transparent method-development prototype.
# - p-values after sparse selection are approximate and should be treated as
#   screening/ranking evidence, not final publication-grade selective inference.

# -----------------------------------------------------------------------------
# Preprocessing helpers
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

validate_fdr_method <- function(fdr_method) {
  match.arg(fdr_method, choices = c("BH", "BY"))
}

validate_fdr_scope <- function(fdr_scope) {
  match.arg(fdr_scope, choices = c("global", "within_view"))
}

validate_primary_inference <- function(primary_inference) {
  match.arg(primary_inference, choices = c("model_based", "bootstrap_score"))
}

validate_study_design <- function(study_design) {
  match.arg(study_design, choices = c("standard", "case_control", "time", "case_control_time"))
}

validate_exposure_role <- function(exposure_role) {
  match.arg(exposure_role, choices = c("time", "case", "interaction"))
}

zentangler_nonmissing_unique <- function(x) {
  out <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  sort(out)
}

zentangler_binary_indicator <- function(x, ref = NULL, compare = NULL, var_label = "variable") {
  vals <- zentangler_nonmissing_unique(x)
  if (is.null(ref) || is.null(compare)) {
    if (length(vals) != 2L) {
      stop(
        var_label, " must have exactly two non-missing levels when ref/compare are not supplied. ",
        "Observed levels: ", paste(vals, collapse = ", "),
        call. = FALSE
      )
    }
    if (is.null(ref)) ref <- vals[[1L]]
    if (is.null(compare)) compare <- vals[[2L]]
  }
  ref <- as.character(ref)
  compare <- as.character(compare)
  if (identical(ref, compare)) stop(var_label, " ref and compare levels must be different.", call. = FALSE)

  z <- as.character(x)
  out <- rep(NA_real_, length(z))
  out[z == ref] <- 0
  out[z == compare] <- 1
  out
}

prepare_zentangler_study_design <- function(
  pheno_df,
  x_var = NULL,
  covariates = NULL,
  study_design = c("standard", "case_control", "time", "case_control_time"),
  case_var = NULL,
  case_level = NULL,
  control_level = NULL,
  time_var = NULL,
  time_ref = NULL,
  time_compare = NULL,
  exposure_role = c("time", "case", "interaction"),
  add_design_covariates = TRUE,
  effect_x0 = NULL,
  effect_x1 = NULL
) {
  study_design <- validate_study_design(study_design)
  exposure_role <- validate_exposure_role(exposure_role)
  covariates <- unique(covariates)
  design_info <- list(study_design = study_design)

  if (identical(study_design, "standard")) {
    if (is.null(x_var) || length(x_var) != 1L || !nzchar(as.character(x_var))) {
      stop("x_var is required when study_design = 'standard'.", call. = FALSE)
    }
    return(list(
      pheno_df = pheno_df,
      x_var = as.character(x_var),
      covariates = covariates,
      effect_x0 = effect_x0,
      effect_x1 = effect_x1,
      design_info = design_info
    ))
  }

  if (identical(study_design, "case_control") || identical(study_design, "case_control_time")) {
    if (is.null(case_var) || !(case_var %in% colnames(pheno_df))) {
      stop("case_var must name a column in phenotype data for study_design = '", study_design, "'.", call. = FALSE)
    }
    pheno_df$.zentangler_case <- zentangler_binary_indicator(
      pheno_df[[case_var]],
      ref = control_level,
      compare = case_level,
      var_label = case_var
    )
    design_info$case_var <- case_var
    design_info$case_level <- case_level %||% zentangler_nonmissing_unique(pheno_df[[case_var]])[[2L]]
    design_info$control_level <- control_level %||% zentangler_nonmissing_unique(pheno_df[[case_var]])[[1L]]
  }

  if (identical(study_design, "time") || identical(study_design, "case_control_time")) {
    if (is.null(time_var) || !(time_var %in% colnames(pheno_df))) {
      stop("time_var must name a column in phenotype data for study_design = '", study_design, "'.", call. = FALSE)
    }
    pheno_df$.zentangler_time <- zentangler_binary_indicator(
      pheno_df[[time_var]],
      ref = time_ref,
      compare = time_compare,
      var_label = time_var
    )
    design_info$time_var <- time_var
    design_info$time_ref <- time_ref %||% zentangler_nonmissing_unique(pheno_df[[time_var]])[[1L]]
    design_info$time_compare <- time_compare %||% zentangler_nonmissing_unique(pheno_df[[time_var]])[[2L]]
  }

  if (identical(study_design, "case_control")) {
    x_var <- ".zentangler_case"
  } else if (identical(study_design, "time")) {
    x_var <- ".zentangler_time"
  } else {
    design_info$exposure_role <- exposure_role
    if (identical(exposure_role, "time")) {
      x_var <- ".zentangler_time"
      if (isTRUE(add_design_covariates)) covariates <- unique(c(covariates, ".zentangler_case"))
    } else if (identical(exposure_role, "case")) {
      x_var <- ".zentangler_case"
      if (isTRUE(add_design_covariates)) covariates <- unique(c(covariates, ".zentangler_time"))
    } else {
      pheno_df$.zentangler_case_time <- pheno_df$.zentangler_case * pheno_df$.zentangler_time
      x_var <- ".zentangler_case_time"
      if (isTRUE(add_design_covariates)) covariates <- unique(c(covariates, ".zentangler_case", ".zentangler_time"))
    }
  }

  if (is.null(effect_x0)) effect_x0 <- 0
  if (is.null(effect_x1)) effect_x1 <- 1
  design_info$derived_x_var <- x_var
  design_info$add_design_covariates <- isTRUE(add_design_covariates)

  list(
    pheno_df = pheno_df,
    x_var = x_var,
    covariates = covariates,
    effect_x0 = effect_x0,
    effect_x1 = effect_x1,
    design_info = design_info
  )
}

apply_primary_fdr_scope <- function(tab, fdr_method = c("BH", "BY"), fdr_scope = c("global", "within_view")) {
  fdr_method <- validate_fdr_method(fdr_method)
  fdr_scope <- validate_fdr_scope(fdr_scope)

  if (is.null(tab) || nrow(tab) == 0) return(tab)
  if (!("p_primary" %in% colnames(tab))) stop("tab must contain p_primary.")

  tab$q_primary_global <- p.adjust(tab$p_primary, method = fdr_method)
  tab$q_primary_within_view <- NA_real_

  if ("omics" %in% colnames(tab)) {
    idx_by_view <- split(seq_len(nrow(tab)), as.character(tab$omics))
  } else {
    idx_by_view <- list(all = seq_len(nrow(tab)))
  }

  for (idx in idx_by_view) {
    tab$q_primary_within_view[idx] <- p.adjust(tab$p_primary[idx], method = fdr_method)
  }

  tab$q_primary <- if (identical(fdr_scope, "global")) {
    tab$q_primary_global
  } else {
    tab$q_primary_within_view
  }

  tab
}

as_numeric_matrix <- function(df, block_name = "block") {
  if (is.matrix(df)) {
    if (!is.numeric(df)) stop(block_name, " must be numeric.")
    M <- df
  } else if (is.data.frame(df)) {
    keep <- vapply(df, is.numeric, logical(1))
    if (!any(keep)) stop(block_name, ": no numeric columns found.")

    dropped <- setdiff(colnames(df), colnames(df)[keep])
    if (length(dropped) > 0) {
      message(block_name, ": dropping non-numeric columns: ", paste(dropped, collapse = ", "))
    }

    M <- as.matrix(df[, keep, drop = FALSE])
  } else {
    stop(block_name, " must be a matrix or data.frame.")
  }

  storage.mode(M) <- "double"
  M
}

drop_zero_variance <- function(M, block_name = "block", tol = 1e-12) {
  sds <- apply(M, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & (sds > tol)
  if (!any(keep)) stop(block_name, ": all columns removed after zero-variance filtering.")

  dropped <- colnames(M)[!keep]
  if (length(dropped) > 0) {
    message(block_name, ": dropping zero-variance columns: ", paste(dropped, collapse = ", "))
  }

  M[, keep, drop = FALSE]
}

drop_nonfinite_columns <- function(M, block_name = "block") {
  bad <- colSums(!is.finite(M)) > 0
  if (any(bad)) {
    message(block_name, ": dropping columns with non-finite values: ", paste(colnames(M)[bad], collapse = ", "))
  }

  out <- M[, !bad, drop = FALSE]
  if (ncol(out) == 0) stop(block_name, ": all columns removed after non-finite filtering.")
  out
}

make_covariate_matrix <- function(pheno_df, covariates = NULL) {
  if (is.null(covariates) || length(covariates) == 0) return(NULL)

  missing <- setdiff(covariates, colnames(pheno_df))
  if (length(missing) > 0) {
    stop("Missing covariates in pheno_df: ", paste(missing, collapse = ", "))
  }

  mm <- model.matrix(~ ., data = pheno_df[, covariates, drop = FALSE])
  if (ncol(mm) <= 1) return(NULL)
  mm[, -1, drop = FALSE]
}

residualize_matrix <- function(M, C) {
  if (is.null(C) || ncol(C) == 0) return(as.matrix(M))

  Xc <- cbind(`(Intercept)` = 1, as.matrix(C))
  fit <- lm.fit(x = Xc, y = as.matrix(M))
  as.matrix(M) - Xc %*% fit$coefficients
}

sanitize_view_names <- function(x) {
  out <- make.names(as.character(x), unique = TRUE)
  out <- gsub("\\.+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  out[!nzchar(out)] <- paste0("view", seq_len(sum(!nzchar(out))))
  make.unique(out, sep = "_")
}

align_samples_multiview_blocks <- function(blocks, pheno_df, needed_pheno_cols) {
  if (!is.list(blocks) || length(blocks) < 1) {
    stop("blocks must be a named list of sample-by-feature matrices/data.frames.")
  }
  if (is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
    names(blocks) <- paste0("view", seq_along(blocks))
  }
  names(blocks) <- sanitize_view_names(names(blocks))

  if (is.null(rownames(pheno_df))) stop("pheno_df must have sample IDs as rownames.")
  if (any(vapply(blocks, function(x) is.null(rownames(x)), logical(1)))) {
    stop("Every block must have sample IDs as rownames.")
  }

  common_ids <- Reduce(intersect, c(lapply(blocks, rownames), list(rownames(pheno_df))))
  if (length(common_ids) == 0) stop("No overlapping sample IDs across blocks and pheno_df.")

  P <- pheno_df[common_ids, , drop = FALSE]
  cc <- complete.cases(P[, needed_pheno_cols, drop = FALSE])
  if (!all(cc)) message("Dropping ", sum(!cc), " samples with missing X/Y/covariate values.")
  ids <- common_ids[cc]

  out_blocks <- lapply(names(blocks), function(nm) {
    M <- as_numeric_matrix(blocks[[nm]][ids, , drop = FALSE], paste0("blocks$", nm))
    M <- drop_zero_variance(M, paste0("blocks$", nm))
    M <- drop_nonfinite_columns(M, paste0("blocks$", nm))
    M
  })
  names(out_blocks) <- names(blocks)

  list(sample_ids = ids, blocks = out_blocks, pheno = P[ids, , drop = FALSE])
}

prefix_multiview_blocks <- function(blocks) {
  mats <- list()
  map <- list()

  for (view in names(blocks)) {
    M <- as.matrix(blocks[[view]])
    raw_names <- colnames(M)
    design_names <- paste0(view, "::", raw_names)
    colnames(M) <- design_names
    mats[[view]] <- M
    map[[view]] <- data.frame(
      view = view,
      mediator = raw_names,
      design_name = design_names,
      stringsAsFactors = FALSE
    )
  }

  list(X = do.call(cbind, mats), map = do.call(rbind, map))
}
