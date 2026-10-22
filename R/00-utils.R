# Shared utilities: validation, study design, sample alignment, FDR.
#
# Method overview (see zentangler.R for the file map):
# - A-stage: M_j ~ X + C for every feature (lm), or MaAsLin2 with any of its
#   options through maaslin2_args.
# - Screen: top sis_n features per view by standardized |a| (or p-value).
# - Outcome model: multiview cooperative learning, rho tuned in [0, 1]
#   (0 = early, 1 = late fusion); IntegratedLearner-style early/late fusion is
#   available as a supplementary option.
# - B-path: de-biased lasso started from the fused fit, or a full-pipeline
#   bootstrap. Untested features enter the FDR step with p = 1.
# - FDR: HDMT-style composite-null procedure for H0: a = 0 or b = 0.

# -----------------------------------------------------------------------------
# Preprocessing helpers
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# Capture the user's global RNG state so a fit can call set.seed() internally
# and restore the caller's stream afterwards:
#   restore_rng <- zentangler_save_rng(); on.exit(restore_rng(), add = TRUE)
zentangler_save_rng <- function() {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}

zentangler_muffle_cox_ties <- function(expr) {
  # glmnet >= 5.1 announces a future change of its default Cox tie handling;
  # Zentangler fixes Breslow ties internally, so the notice is not actionable.
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("tie-handling|cox.ties", conditionMessage(w))) invokeRestart("muffleWarning")
  })
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
  # Natural level order: factor levels, then numeric order, then character
  # order. Only factor and numeric orders are used to choose reference levels
  # automatically; see zentangler_binary_indicator().
  keep <- !is.na(x) & nzchar(as.character(x))
  if (is.factor(x)) {
    lv <- levels(droplevels(x[keep]))
    return(as.character(lv))
  }
  if (is.numeric(x)) {
    return(as.character(sort(unique(x[keep]))))
  }
  sort(unique(as.character(x[keep])))
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
    if (!is.factor(x) && !is.numeric(x)) {
      # Alphabetical order silently decides the sign of every effect for
      # character columns (for example "case" sorts before "control"), so the
      # reference level must be stated explicitly.
      stop(
        var_label, " is a character column; supply the reference and comparison levels explicitly ",
        "(control_level/case_level or time_ref/time_compare), or store it as a factor whose first ",
        "level is the reference. Observed levels: ", paste(vals, collapse = ", "),
        call. = FALSE
      )
    }
    if (is.null(ref)) ref <- setdiff(vals, as.character(compare))[[1L]]
    if (is.null(compare)) compare <- setdiff(vals, as.character(ref))[[1L]]
  }
  ref <- as.character(ref)
  compare <- as.character(compare)
  if (identical(ref, compare)) stop(var_label, " ref and compare levels must be different.", call. = FALSE)
  z <- as.character(x)
  if (!(ref %in% z)) stop(var_label, ": reference level '", ref, "' not found.", call. = FALSE)
  if (!(compare %in% z)) stop(var_label, ": comparison level '", compare, "' not found.", call. = FALSE)

  out <- rep(NA_real_, length(z))
  out[z == ref] <- 0
  out[z == compare] <- 1
  attr(out, "ref") <- ref
  attr(out, "compare") <- compare
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
    case_ind <- zentangler_binary_indicator(
      pheno_df[[case_var]],
      ref = control_level,
      compare = case_level,
      var_label = case_var
    )
    design_info$case_var <- case_var
    design_info$case_level <- attr(case_ind, "compare")
    design_info$control_level <- attr(case_ind, "ref")
    pheno_df$.zentangler_case <- as.numeric(case_ind)
  }

  if (identical(study_design, "time") || identical(study_design, "case_control_time")) {
    if (is.null(time_var) || !(time_var %in% colnames(pheno_df))) {
      stop("time_var must name a column in phenotype data for study_design = '", study_design, "'.", call. = FALSE)
    }
    time_ind <- zentangler_binary_indicator(
      pheno_df[[time_var]],
      ref = time_ref,
      compare = time_compare,
      var_label = time_var
    )
    design_info$time_var <- time_var
    design_info$time_ref <- attr(time_ind, "ref")
    design_info$time_compare <- attr(time_ind, "compare")
    pheno_df$.zentangler_time <- as.numeric(time_ind)
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

validate_fdr_joint <- function(fdr_joint) {
  match.arg(fdr_joint, choices = c("hdmt", "maxp"))
}

zentangler_storey_pi0 <- function(p, lambda = 0.5) {
  p <- p[is.finite(p)]
  if (length(p) == 0L) return(1)
  min(1, mean(p > lambda) / (1 - lambda))
}

zentangler_hdmt_qvalues <- function(p_a, p_b, p_b_marginal = NULL, lambda = 0.5, lambda_a = 0.2, gamma_conf = 0.9, min_null_features = 50L) {
  # Composite-null FDR for H0: a = 0 or b = 0 (joint significance), in the
  # spirit of HDMT (Dai, Stanford & LeBlanc 2020, JASA), adapted to a screened
  # design in which the conditional p_b exists only for the B-stage-tested
  # set T.
  #
  # A feature is discovered at threshold t when it is in T and
  # max(p_a, p_b) <= t. Expected false discoveries are bounded by
  #   (i)  b-null features in T:    t * #{j in T : p_a,j <= t}
  #        (given the screened set, p_b of a b-null feature is uniform and
  #        independent of p_a; pi0_b(T) is bounded by 1 because T is small)
  #   (ii) a-null, b-non-null:      m * pi0_a * t * gamma * F1_b(t)
  #        where gamma is the share of b-non-null features among a-null
  #        features and F1_b the CDF of non-null p_b (estimated within T).
  # gamma is bounded from above with marginal outcome p-values
  # (Y ~ X + C + M_j) of the features with p_a > lambda_a (mostly a-null):
  # an upper `gamma_conf` confidence bound of 1 - pi0 from the binomial count
  # of marginal p-values above lambda. Marginal association includes every
  # conditional effect plus correlation with other mediators, which inflates
  # the bound further (conservative). Without marginal p-values gamma = 1.
  # pi0_a uses Storey's estimator. Untested features get q = 1.
  m_all <- length(p_a)
  q <- rep(1, m_all)
  tested <- is.finite(p_a) & is.finite(p_b)
  if (!any(tested)) return(list(q = q, gamma = NA_real_, pi0a = NA_real_, pi0b = NA_real_))
  m_a <- sum(is.finite(p_a))
  pi0a <- zentangler_storey_pi0(p_a, lambda)
  pb_t <- p_b[tested]
  pa_t <- p_a[tested]
  pi0b <- zentangler_storey_pi0(pb_t, lambda)

  gamma <- 1
  if (!is.null(p_b_marginal)) {
    null_a <- is.finite(p_a) & p_a > lambda_a & is.finite(p_b_marginal)
    m0 <- sum(null_a)
    if (m0 >= min_null_features) {
      k_hi <- sum(p_b_marginal[null_a] > lambda)
      frac_lo <- stats::qbeta(1 - gamma_conf, k_hi, m0 - k_hi + 1)
      if (!is.finite(frac_lo)) frac_lo <- 0
      gamma <- min(1, max(1 - frac_lo / (1 - lambda), 1 / m0))
    }
  }

  maxp <- pmax(pa_t, pb_t)
  t_grid <- sort(unique(maxp))
  Fb <- stats::ecdf(pb_t)(t_grid)
  F1 <- if (pi0b < 0.999) (Fb - pi0b * t_grid) / (1 - pi0b) else rep(1, length(t_grid))
  F1 <- pmin(1, pmax(t_grid, cummax(pmax(F1, 0))))
  R <- vapply(t_grid, function(t) sum(maxp <= t), numeric(1))
  n_pa <- vapply(t_grid, function(t) sum(pa_t <= t), numeric(1))
  efd <- t_grid * n_pa + m_a * pi0a * gamma * t_grid * F1
  fdr <- pmin(1, efd / pmax(R, 1))
  fdr <- rev(cummin(rev(fdr)))
  q[which(tested)] <- fdr[match(maxp, t_grid)]
  list(q = q, gamma = gamma, pi0a = pi0a, pi0b = pi0b)
}

zentangler_marginal_b_pvalues <- function(Y, X, blocks, C = NULL, y_family = c("gaussian", "binomial", "survival")) {
  # Vectorized score tests of M_j in Y ~ X + C + M_j, one per feature, for the
  # HDMT-style FDR step. Returns per-view named vectors.
  y_family <- match.arg(y_family)
  W <- zentangler_y_design(X, C)
  W1 <- cbind(1, W)
  base <- zentangler_unpenalized_fit(Y, W, y_family)
  eta0 <- base$intercept + as.numeric(W %*% base$coef)
  out <- lapply(blocks, function(M) {
    M <- as.matrix(M)
    if (identical(y_family, "survival")) {
      time <- as.numeric(Y[, 1])
      status <- as.numeric(Y[, 2])
      ord <- order(time, decreasing = TRUE)
      r <- exp(eta0 - max(eta0))[ord]
      Mo <- M[ord, , drop = FALSE]
      to <- time[ord]
      so <- status[ord]
      S0 <- 0
      S1 <- rep(0, ncol(M))
      S2 <- rep(0, ncol(M))
      U <- rep(0, ncol(M))
      V <- rep(0, ncol(M))
      pos <- 1L
      n <- nrow(M)
      for (t in unique(to)) {
        idx <- integer(0)
        while (pos <= n && to[pos] == t) {
          idx <- c(idx, pos)
          pos <- pos + 1L
        }
        mi <- Mo[idx, , drop = FALSE]
        S0 <- S0 + sum(r[idx])
        S1 <- S1 + colSums(mi * r[idx])
        S2 <- S2 + colSums(mi^2 * r[idx])
        ev <- idx[so[idx] == 1]
        d <- length(ev)
        if (d == 0L) next
        mbar <- S1 / S0
        U <- U + colSums(Mo[ev, , drop = FALSE]) - d * mbar
        V <- V + d * (S2 / S0 - mbar^2)
      }
      z <- U / sqrt(pmax(V, 1e-12))
    } else {
      if (identical(y_family, "binomial")) {
        mu <- stats::plogis(eta0)
        w <- mu * (1 - mu)
        resid <- as.numeric(Y) - mu
      } else {
        w <- rep(1, length(eta0))
        resid <- as.numeric(Y) - eta0
      }
      sw <- sqrt(w)
      Mt <- M - W1 %*% qr.coef(qr(W1 * sw), M * sw)
      U <- colSums(Mt * resid)
      sigma2 <- if (identical(y_family, "gaussian")) sum(resid^2) / (length(resid) - ncol(W1)) else 1
      V <- sigma2 * colSums(Mt^2 * w)
      z <- U / sqrt(pmax(V, 1e-12))
    }
    setNames(2 * stats::pnorm(-abs(z)), colnames(M))
  })
  names(out) <- names(blocks)
  out
}

apply_primary_fdr_scope <- function(tab, fdr_method = c("BH", "BY"), fdr_scope = c("global", "within_view"), fdr_joint = c("hdmt", "maxp"), joint_ok = TRUE) {
  # q-values for p_primary. With fdr_joint = "hdmt" (default) and model-based
  # joint evidence (joint_ok), q-values use the composite-null estimator
  # above; otherwise BH/BY is applied to p_primary.
  fdr_method <- validate_fdr_method(fdr_method)
  fdr_scope <- validate_fdr_scope(fdr_scope)
  fdr_joint <- validate_fdr_joint(fdr_joint)

  if (is.null(tab) || nrow(tab) == 0) return(tab)
  if (!("p_primary" %in% colnames(tab))) stop("tab must contain p_primary.")
  use_hdmt <- identical(fdr_joint, "hdmt") && isTRUE(joint_ok) && all(c("p_a", "p_b") %in% colnames(tab))
  qfun <- function(idx) {
    if (use_hdmt) {
      pb <- tab$p_b[idx]
      if ("tested_b_stage" %in% colnames(tab)) pb[!tab$tested_b_stage[idx]] <- NA_real_
      pm <- if ("p_b_marginal" %in% colnames(tab)) tab$p_b_marginal[idx] else NULL
      zentangler_hdmt_qvalues(tab$p_a[idx], pb, p_b_marginal = pm)$q
    } else {
      p.adjust(tab$p_primary[idx], method = fdr_method)
    }
  }

  tab$q_primary_global <- qfun(seq_len(nrow(tab)))
  tab$q_primary_within_view <- NA_real_

  if ("omics" %in% colnames(tab)) {
    idx_by_view <- split(seq_len(nrow(tab)), as.character(tab$omics))
  } else {
    idx_by_view <- list(all = seq_len(nrow(tab)))
  }

  for (idx in idx_by_view) {
    tab$q_primary_within_view[idx] <- qfun(idx)
  }

  tab$q_primary <- if (identical(fdr_scope, "global")) {
    tab$q_primary_global
  } else {
    tab$q_primary_within_view
  }
  attr(tab, "fdr_joint_used") <- if (use_hdmt) "hdmt" else paste0("maxp_", fdr_method)
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
  mm <- mm[, -1, drop = FALSE]
  # Prefix covariate columns so they can never collide with the exposure column
  # ("X") or with mediator design names inside the model matrices.
  colnames(mm) <- make.names(paste0("cov_", colnames(mm)), unique = TRUE)
  mm
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
    M <- drop_nonfinite_columns(M, paste0("blocks$", nm))
    M <- drop_zero_variance(M, paste0("blocks$", nm))
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

# -----------------------------------------------------------------------------
# Logging (verbose = TRUE)
# -----------------------------------------------------------------------------

zentangler_log <- function(verbose, ..., indent = 0L) {
  # Timestamped progress messages through message(), so they can be silenced
  # with verbose = FALSE or suppressMessages().
  if (!isTRUE(verbose)) return(invisible(NULL))
  message(sprintf("[Zentangler %s] ", format(Sys.time(), "%H:%M:%S")), strrep("  ", indent), ...)
  invisible(NULL)
}

zentangler_fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  trimws(if (p < 1e-4) formatC(p, format = "e", digits = 2) else formatC(p, format = "g", digits = 3))
}

zentangler_progress_every <- function(n) max(1L, as.integer(ceiling(n / 10)))

zentangler_exposure_outcome_test <- function(Y, X, C, y_family, x0, x1) {
  # Total-effect check before mediation: Y ~ X + covariates (lm / logistic /
  # Cox). Mediation can exist without a significant total effect (e.g. paths
  # of opposite sign), so this is reported, not used as a gate.
  r <- estimate_reduced_x_summary(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
  p <- r$p_reduced_total_effect[[1L]]
  data.frame(
    model = switch(y_family, gaussian = "linear", binomial = "logistic", survival = "Cox"),
    n = length(as.numeric(X)),
    n_covariates = if (is.null(C)) 0L else ncol(C),
    coefficient = r$reduced_total_effect_coefficient[[1L]],
    p_value = p,
    effect = r$reduced_total_effect[[1L]],
    effect_scale = switch(y_family, gaussian = "mean_difference", binomial = "risk_difference", survival = "log_hazard_ratio"),
    x0 = x0,
    x1 = x1,
    associated_p05 = is.finite(p) && p < 0.05,
    stringsAsFactors = FALSE
  )
}

zentangler_log_exposure_outcome <- function(verbose, eo) {
  if (!isTRUE(verbose) || is.null(eo) || nrow(eo) == 0L) return(invisible(NULL))
  zentangler_log(verbose, sprintf(
    "Exposure-outcome association (%s model, Y ~ X%s): coefficient %.4g, p = %s; effect of X %.3g -> %.3g = %.4g (%s).",
    eo$model, if (eo$n_covariates > 0L) " + covariates" else "", eo$coefficient, zentangler_fmt_p(eo$p_value),
    eo$x0, eo$x1, eo$effect, gsub("_", " ", eo$effect_scale)
  ))
  if (isTRUE(eo$associated_p05)) {
    zentangler_log(verbose, "Exposure and outcome are associated (p < 0.05).", indent = 1L)
  } else {
    zentangler_log(verbose, "Exposure and outcome are NOT significantly associated (p >= 0.05). Indirect effects can still exist (e.g. mediated paths of opposite sign), but interpret the total effect and proportion mediated with care.", indent = 1L)
  }
  invisible(NULL)
}
