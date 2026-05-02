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
# - A-stage: per-view HIMA-like univariate screen, M_j^(k) ~ X + C.
# - SIS: keep top mediators per view.
# - B/Y-stage: fit selected mediators jointly using early, intermediate, or late
#   lasso-based fusion.
# - Score: a_j * b_j, reported for every mediator with view labels preserved.
#
# Scientific status:
# - This is a transparent method-development prototype.
# - p-values after sparse selection are approximate and should be treated as
#   screening/ranking evidence, not final publication-grade selective inference.

# -----------------------------------------------------------------------------
# Preprocessing helpers
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# MultiAssayExperiment helpers
# -----------------------------------------------------------------------------

zentangler_require_mae <- function() {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required for Zentangler MAE input.", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required for Zentangler MAE input.", call. = FALSE)
  }
  invisible(TRUE)
}

mae_fill_nonfinite_zero <- function(df_like) {
  M <- as.matrix(df_like)
  M[!is.finite(M)] <- 0
  out <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- rownames(df_like)
  out
}

mae_rows_mean_by_id <- function(df_like) {
  if (nrow(df_like) == 0) return(df_like)
  M <- as.matrix(df_like)
  grp <- rownames(df_like)
  rs <- rowsum(M, group = grp, reorder = FALSE)
  n_grp <- as.numeric(table(grp)[rownames(rs)])
  out <- rs / n_grp
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

zentangler_mae_coldata <- function(mae) {
  zentangler_require_mae()
  cd <- as.data.frame(SummarizedExperiment::colData(mae), stringsAsFactors = FALSE)
  if (is.null(rownames(cd)) || any(!nzchar(rownames(cd)))) {
    stop("MAE colData must have rownames as primary sample IDs.", call. = FALSE)
  }
  cd
}

zentangler_mae_extract_block <- function(
  mae,
  experiment_name,
  assay_name = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  # Return one MAE experiment as samples x features.
  #
  # SummarizedExperiment assays are Bioconductor-standard features x samples, so
  # we transpose them. Plain matrices/data.frames are accepted too; if their
  # columns match MAE samples/sampleMap, they are treated as features x samples,
  # otherwise rows are treated as samples.
  duplicate_primary <- match.arg(duplicate_primary)
  zentangler_require_mae()

  exps <- MultiAssayExperiment::experiments(mae)
  if (!(experiment_name %in% names(exps))) {
    stop("experiment_name not found in MAE experiments: ", experiment_name, call. = FALSE)
  }

  exp_obj <- exps[[experiment_name]]
  if (inherits(exp_obj, "SummarizedExperiment")) {
    a_names <- SummarizedExperiment::assayNames(exp_obj)
    if (length(a_names) == 0) stop("No assays found in experiment: ", experiment_name, call. = FALSE)
    if (is.null(assay_name)) assay_name <- a_names[1]
    if (!(assay_name %in% a_names)) {
      stop("assay_name '", assay_name, "' not found in experiment '", experiment_name, "'.", call. = FALSE)
    }
    mat <- as.matrix(SummarizedExperiment::assay(exp_obj, assay_name))
    if (is.null(colnames(mat))) {
      stop("Experiment '", experiment_name, "' assay has no sample colnames.", call. = FALSE)
    }
    if (is.null(rownames(mat))) rownames(mat) <- paste0(experiment_name, "_feature_", seq_len(nrow(mat)))
    block <- as.data.frame(t(mat), check.names = FALSE, stringsAsFactors = FALSE)
    rownames(block) <- colnames(mat)
  } else if (is.matrix(exp_obj) || is.data.frame(exp_obj)) {
    mat <- as.matrix(exp_obj)
    if (is.null(rownames(mat))) rownames(mat) <- paste0(experiment_name, "_row_", seq_len(nrow(mat)))

    primary_ids <- rownames(zentangler_mae_coldata(mae))
    sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    mapped_cols <- character(0)
    if (all(c("assay", "colname") %in% colnames(sm))) {
      mapped_cols <- as.character(sm$colname[sm$assay == experiment_name])
    }
    looks_feature_by_sample <- !is.null(colnames(mat)) && (
      any(colnames(mat) %in% primary_ids) || any(colnames(mat) %in% mapped_cols)
    )

    if (looks_feature_by_sample) {
      block <- as.data.frame(t(mat), check.names = FALSE, stringsAsFactors = FALSE)
      rownames(block) <- colnames(mat)
    } else {
      if (is.null(rownames(mat))) {
        stop("Plain matrix experiment '", experiment_name, "' needs rownames as samples or colnames as samples.", call. = FALSE)
      }
      block <- as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE)
    }
  } else {
    stop(
      "Unsupported experiment class for '", experiment_name, "': ",
      paste(class(exp_obj), collapse = ", "),
      ". Use SummarizedExperiment-like or matrix/data.frame.",
      call. = FALSE
    )
  }

  sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  if (all(c("assay", "colname", "primary") %in% colnames(sm))) {
    smi <- sm[sm$assay == experiment_name, c("colname", "primary"), drop = FALSE]
    if (nrow(smi) > 0) {
      idx <- match(rownames(block), as.character(smi$colname))
      pri <- as.character(smi$primary[idx])
      use_map <- !is.na(pri) & nzchar(trimws(pri))
      rownames(block)[use_map] <- pri[use_map]
    }
  }

  keep_rows <- !is.na(rownames(block)) & nzchar(trimws(rownames(block)))
  block <- block[keep_rows, , drop = FALSE]
  if (nrow(block) == 0) stop("No valid sample rows after extracting experiment: ", experiment_name, call. = FALSE)

  if (anyDuplicated(rownames(block)) > 0) {
    if (identical(duplicate_primary, "first")) {
      block <- block[!duplicated(rownames(block)), , drop = FALSE]
    } else {
      block <- mae_rows_mean_by_id(block)
    }
  }

  raw <- as.matrix(block)
  num <- suppressWarnings(matrix(
    as.numeric(raw),
    nrow = nrow(raw),
    ncol = ncol(raw),
    dimnames = dimnames(raw)
  ))
  block_num <- as.data.frame(num, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(block_num) <- rownames(block)
  if (fill_nonfinite_zero) block_num <- mae_fill_nonfinite_zero(block_num)
  block_num
}

zentangler_mae_to_blocks <- function(
  mae,
  view_names = NULL,
  assay_names = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  duplicate_primary <- match.arg(duplicate_primary)
  zentangler_require_mae()

  exps <- MultiAssayExperiment::experiments(mae)
  if (is.null(view_names)) view_names <- names(exps)
  if (length(view_names) == 0) stop("MAE contains no experiments/views.", call. = FALSE)
  missing_views <- setdiff(view_names, names(exps))
  if (length(missing_views) > 0) {
    stop("Requested views not found in MAE: ", paste(missing_views, collapse = ", "), call. = FALSE)
  }

  if (is.null(assay_names)) {
    assay_names <- setNames(rep(NA_character_, length(view_names)), view_names)
  } else {
    if (is.null(names(assay_names))) {
      if (length(assay_names) != length(view_names)) {
        stop("Unnamed assay_names must have the same length as view_names.", call. = FALSE)
      }
      names(assay_names) <- view_names
    }
    assay_names <- assay_names[view_names]
  }

  safe_names <- sanitize_view_names(view_names)
  blocks <- vector("list", length(view_names))
  names(blocks) <- safe_names
  view_map <- data.frame(
    experiment_name = view_names,
    view = safe_names,
    assay_name = unname(assay_names),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(view_names)) {
    an <- assay_names[[view_names[i]]]
    if (is.na(an) || !nzchar(an)) an <- NULL
    blocks[[safe_names[i]]] <- zentangler_mae_extract_block(
      mae = mae,
      experiment_name = view_names[i],
      assay_name = an,
      duplicate_primary = duplicate_primary,
      fill_nonfinite_zero = fill_nonfinite_zero
    )
  }

  list(blocks = blocks, pheno_df = zentangler_mae_coldata(mae), view_map = view_map)
}

# -----------------------------------------------------------------------------
# A-stage screening
# -----------------------------------------------------------------------------

hima_sis_default_n <- function(n, p) {
  max(1L, min(p, floor(n / max(log(n), 1.0))))
}

hima_screen_mediators <- function(
  M,
  X,
  C = NULL,
  sis_n = NULL,
  rank_by = c("abs_a", "pvalue"),
  screen_method = c("sis", "none")
) {
  rank_by <- match.arg(rank_by)
  screen_method <- match.arg(screen_method)

  n <- nrow(M)
  p <- ncol(M)

  tab <- data.frame(
    mediator = colnames(M),
    a_hat = NA_real_,
    p_a = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    dat <- data.frame(y = M[, j], X = X)
    if (!is.null(C) && ncol(C) > 0) dat <- cbind(dat, as.data.frame(C))

    fit <- try(lm(y ~ ., data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) next

    co <- coef(summary(fit))
    if ("X" %in% rownames(co)) {
      tab$a_hat[j] <- co["X", "Estimate"]
      tab$p_a[j] <- co["X", "Pr(>|t|)"]
    }
  }

  tab$q_a <- p.adjust(tab$p_a, method = "BH")
  ord <- if (rank_by == "abs_a") {
    order(abs(tab$a_hat), decreasing = TRUE, na.last = TRUE)
  } else {
    order(tab$p_a, na.last = TRUE)
  }

  if (screen_method == "none") {
    # Keep every mediator that survived sample alignment and numeric filtering.
    # A-path estimates are still computed because they are part of the mediation score.
    sis_n_used <- p
    selected <- tab$mediator
  } else {
    if (is.null(sis_n)) sis_n <- hima_sis_default_n(n = n, p = p)
    sis_n_used <- max(1L, min(p, as.integer(sis_n)))
    selected <- tab$mediator[ord][seq_len(sis_n_used)]
  }

  tab$selected <- tab$mediator %in% selected
  tab <- tab[order(tab$p_a, na.last = TRUE), ]
  rownames(tab) <- NULL

  list(table = tab, selected = selected, sis_n = sis_n_used, screen_method = screen_method)
}

# -----------------------------------------------------------------------------
# Fusion / Y-stage helpers
# -----------------------------------------------------------------------------

choose_block_lambda <- function(Y, W, M, y_family, lambda_choice) {
  Z <- cbind(W, M)
  pf <- c(rep(0, ncol(W)), rep(1, ncol(M)))

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = y_family,
    alpha = 1,
    penalty.factor = pf,
    standardize = TRUE,
    intercept = FALSE
  )

  lam <- as.numeric(cvfit[[lambda_choice]])
  if (!is.finite(lam)) lam <- as.numeric(cvfit$lambda[min(10, length(cvfit$lambda))])
  lam <- max(lam, 1e-8)

  list(lambda = lam, cvfit = cvfit)
}

update_gamma_gaussian <- function(Y, W, offset) {
  fit <- lm.fit(x = W, y = Y - offset)
  gamma <- as.numeric(fit$coefficients)
  gamma[!is.finite(gamma)] <- 0
  names(gamma) <- colnames(W)
  gamma
}

fit_block_lasso_for_late <- function(
  Y,
  X,
  M,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  prefix = "B::"
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)

  M2 <- as.matrix(M)
  colnames(M2) <- paste0(prefix, colnames(M))

  Z <- cbind(X = X)
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, C)
  Z <- cbind(Z, M2)

  pf <- rep(1, ncol(Z))
  pf[colnames(Z) == "X"] <- 0
  if (!is.null(C) && ncol(C) > 0) pf[colnames(Z) %in% colnames(C)] <- 0

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = y_family,
    alpha = 1,
    penalty.factor = pf,
    standardize = TRUE
  )

  beta <- as.matrix(coef(cvfit, s = lambda_choice))[, 1]
  names(beta) <- gsub("`", "", names(beta), fixed = TRUE)

  pred <- as.numeric(predict(cvfit, newx = as.matrix(Z), s = lambda_choice, type = "response"))

  x_coef <- beta["X"]
  if (!is.finite(x_coef)) x_coef <- 0

  med_coef <- setNames(rep(0, ncol(M)), colnames(M))
  med_tag <- paste0(prefix, colnames(M))
  hit <- intersect(med_tag, names(beta))
  if (length(hit) > 0) med_coef[sub(paste0("^", prefix), "", hit)] <- beta[hit]

  list(cvfit = cvfit, pred = pred, x_coef = as.numeric(x_coef), med_coef = med_coef)
}

fit_y_multiview_glmnet_lasso <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min")
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)

  pref <- prefix_multiview_blocks(blocks)
  Z <- cbind(X = X)
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, C)
  Z <- cbind(Z, pref$X)

  pf <- rep(1, ncol(Z))
  pf[colnames(Z) == "X"] <- 0
  if (!is.null(C) && ncol(C) > 0) pf[colnames(Z) %in% colnames(C)] <- 0

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = y_family,
    alpha = 1,
    penalty.factor = pf,
    standardize = TRUE
  )

  beta <- as.matrix(coef(cvfit, s = lambda_choice))[, 1]
  names(beta) <- gsub("`", "", names(beta), fixed = TRUE)

  coef_by_view <- lapply(names(blocks), function(view) {
    dn <- paste0(view, "::", colnames(blocks[[view]]))
    out <- setNames(rep(0, length(dn)), colnames(blocks[[view]]))
    hit <- intersect(dn, names(beta))
    if (length(hit) > 0) out[sub(paste0("^", view, "::"), "", hit)] <- beta[hit]
    out
  })
  names(coef_by_view) <- names(blocks)

  x_coef <- beta["X"]
  if (!is.finite(x_coef)) x_coef <- 0

  list(
    method = "K-view early fusion lasso",
    coef_by_view = coef_by_view,
    x_coef = as.numeric(x_coef),
    fit = cvfit,
    feature_map = pref$map
  )
}

fit_y_multiview_cooperative_intermediate <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  rho = 0.2,
  maxit = 100,
  tol = 1e-04
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  if (!identical(y_family, "gaussian")) {
    stop("K-view cooperative intermediate fusion currently supports y_family='gaussian'.")
  }
  if (rho < 0) stop("rho must be non-negative.")

  W <- cbind(`(Intercept)` = 1, X = X)
  if (!is.null(C) && ncol(C) > 0) {
    C2 <- as.matrix(C)
    colnames(C2) <- paste0("C::", colnames(C2))
    W <- cbind(W, C2)
  }
  W <- as.matrix(W)
  blocks <- lapply(blocks, as.matrix)

  n <- length(Y)
  if (any(vapply(blocks, nrow, integer(1)) != n) || nrow(W) != n) {
    stop("Dimension mismatch in K-view cooperative fit.")
  }

  lambda_info <- lapply(blocks, function(M) {
    choose_block_lambda(Y = Y, W = W, M = M, y_family = y_family, lambda_choice = lambda_choice)
  })
  lambdas <- vapply(lambda_info, function(x) x$lambda, numeric(1))

  betas <- lapply(blocks, function(M) setNames(rep(0, ncol(M)), colnames(M)))
  gamma <- update_gamma_gaussian(Y = Y, W = W, offset = rep(0, n))
  converged <- FALSE
  iter <- 0L
  sqrt_rho <- sqrt(rho)

  for (it in seq_len(maxit)) {
    iter <- it
    betas_old <- betas
    gamma_old <- gamma

    eta_list <- Map(function(M, b) as.numeric(M %*% b), blocks, betas)
    offset_old <- Reduce(`+`, eta_list)
    gamma <- update_gamma_gaussian(Y = Y, W = W, offset = offset_old)
    r <- as.numeric(Y - W %*% gamma)

    for (view in names(blocks)) {
      M <- blocks[[view]]
      other_views <- setdiff(names(blocks), view)
      other_eta <- if (length(other_views) > 0) {
        Reduce(`+`, Map(function(M2, b2) as.numeric(M2 %*% b2), blocks[other_views], betas[other_views]))
      } else {
        rep(0, n)
      }
      target_main <- r - other_eta
      target_agree <- if (length(other_views) > 0) other_eta / length(other_views) else rep(0, n)

      if (rho > 0) {
        x_aug <- rbind(M, sqrt_rho * M)
        y_aug <- c(target_main, sqrt_rho * target_agree)
      } else {
        x_aug <- M
        y_aug <- target_main
      }

      fit_k <- glmnet::glmnet(
        x = as.matrix(x_aug),
        y = y_aug,
        family = "gaussian",
        alpha = 1,
        lambda = lambdas[view],
        intercept = FALSE,
        standardize = TRUE
      )
      b_new <- as.matrix(coef(fit_k, s = lambdas[view]))[-1, 1]
      names(b_new) <- colnames(M)
      betas[[view]] <- b_new
    }

    offset_new <- Reduce(`+`, Map(function(M, b) as.numeric(M %*% b), blocks, betas))
    gamma <- update_gamma_gaussian(Y = Y, W = W, offset = offset_new)

    max_db <- max(unlist(Map(function(a, b) max(abs(a - b)), betas, betas_old)))
    max_dg <- if (length(gamma) > 0) max(abs(gamma - gamma_old)) else 0
    delta <- max(max_db, max_dg)
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      break
    }
  }

  x_coef <- if ("X" %in% names(gamma)) as.numeric(gamma["X"]) else 0
  if (!is.finite(x_coef)) x_coef <- 0

  list(
    method = "K-view intermediate fusion cooperative lasso",
    coef_by_view = betas,
    x_coef = x_coef,
    fit = list(
      gamma = gamma,
      lambda = lambdas,
      rho = rho,
      iterations = iter,
      converged = converged,
      lambda_cv = lapply(lambda_info, `[[`, "cvfit")
    )
  )
}

fit_y_multiview_late_fusion_lasso <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min")
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)

  block_fits <- lapply(names(blocks), function(view) {
    fit_block_lasso_for_late(
      Y = Y,
      X = X,
      M = blocks[[view]],
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      prefix = paste0(view, "::")
    )
  })
  names(block_fits) <- names(blocks)

  meta_dat <- as.data.frame(lapply(block_fits, `[[`, "pred"), check.names = FALSE)
  colnames(meta_dat) <- paste0("pred_", names(blocks))
  meta_dat$Y <- Y
  meta_formula <- as.formula(paste("Y ~", paste(setdiff(colnames(meta_dat), "Y"), collapse = " + ")))

  meta_fit <- if (identical(y_family, "gaussian")) {
    lm(meta_formula, data = meta_dat)
  } else {
    glm(meta_formula, data = meta_dat, family = binomial())
  }

  meta_co <- coef(meta_fit)
  weights <- setNames(rep(0, length(blocks)), names(blocks))
  for (view in names(blocks)) {
    nm <- paste0("pred_", view)
    if (nm %in% names(meta_co) && is.finite(meta_co[nm])) weights[view] <- as.numeric(meta_co[nm])
  }

  coef_by_view <- lapply(names(blocks), function(view) block_fits[[view]]$med_coef * weights[view])
  names(coef_by_view) <- names(blocks)
  x_coef <- sum(vapply(names(blocks), function(view) block_fits[[view]]$x_coef * weights[view], numeric(1)))

  list(
    method = "K-view late fusion two-stage lasso",
    coef_by_view = coef_by_view,
    x_coef = as.numeric(x_coef),
    fit = list(blocks = lapply(block_fits, `[[`, "cvfit"), meta = meta_fit),
    meta_weights = weights
  )
}

fit_y_multiview_stage <- function(
  Y,
  X,
  blocks,
  C = NULL,
  fusion_mode = c("early", "intermediate", "late"),
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04
) {
  fusion_mode <- match.arg(fusion_mode)
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)

  if (identical(fusion_mode, "intermediate")) {
    return(fit_y_multiview_cooperative_intermediate(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      rho = coop_rho,
      maxit = coop_maxit,
      tol = coop_tol
    ))
  }

  if (identical(fusion_mode, "late")) {
    return(fit_y_multiview_late_fusion_lasso(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice
    ))
  }

  fit_y_multiview_glmnet_lasso(
    Y = Y,
    X = X,
    blocks = blocks,
    C = C,
    y_family = y_family,
    lambda_choice = lambda_choice
  )
}

# -----------------------------------------------------------------------------
# B-stage p-values and effect summaries
# -----------------------------------------------------------------------------

infer_p_b_multiview_refit <- function(Y, X, blocks, C = NULL, coef_by_view, y_family = c("gaussian", "binomial")) {
  y_family <- match.arg(y_family)
  out <- lapply(blocks, function(M) setNames(rep(NA_real_, ncol(M)), colnames(M)))

  active_design <- character(0)
  active_map <- data.frame(view = character(), mediator = character(), design_name = character())

  for (view in names(blocks)) {
    b <- coef_by_view[[view]]
    active <- names(b)[is.finite(b) & b != 0]
    if (length(active) == 0) next

    dn <- paste0(view, "::", active)
    active_design <- c(active_design, dn)
    active_map <- rbind(
      active_map,
      data.frame(view = view, mediator = active, design_name = dn, stringsAsFactors = FALSE)
    )
  }
  if (length(active_design) == 0) return(out)

  active_views <- unique(active_map$view)
  active_blocks <- lapply(active_views, function(view) {
    blocks[[view]][, active_map$mediator[active_map$view == view], drop = FALSE]
  })
  names(active_blocks) <- active_views

  pref <- prefix_multiview_blocks(active_blocks)
  M_active <- pref$X[, intersect(colnames(pref$X), active_design), drop = FALSE]

  Z <- cbind(`(Intercept)` = 1, X = X)
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, as.matrix(C))
  Z <- cbind(Z, M_active)

  if (identical(y_family, "gaussian")) {
    fit <- lm.fit(x = as.matrix(Z), y = as.numeric(Y))
    rank <- fit$rank
    df <- nrow(Z) - rank
    if (df <= 0) return(out)

    rss <- sum(fit$residuals^2)
    sigma2 <- rss / df
    R <- try(qr.R(fit$qr)[seq_len(rank), seq_len(rank), drop = FALSE], silent = TRUE)
    if (inherits(R, "try-error")) return(out)
    xtx_inv <- try(chol2inv(R), silent = TRUE)
    if (inherits(xtx_inv, "try-error")) return(out)

    beta <- fit$coefficients
    se <- rep(NA_real_, length(beta))
    se[fit$qr$pivot[seq_len(rank)]] <- sqrt(diag(xtx_inv) * sigma2)
    names(se) <- names(beta)
    pvals <- 2 * stats::pt(abs(beta / se), df = df, lower.tail = FALSE)
  } else {
    dat <- as.data.frame(Z, check.names = TRUE)
    safe_names <- colnames(dat)
    fit <- try(suppressWarnings(glm(Y ~ . - 1, data = cbind(Y = Y, dat), family = binomial())), silent = TRUE)
    if (inherits(fit, "try-error")) return(out)
    sm <- summary(fit)$coefficients
    pvals <- sm[, ncol(sm)]
    names(pvals) <- colnames(Z)[match(rownames(sm), safe_names)]
  }

  for (i in seq_len(nrow(active_map))) {
    view <- active_map$view[i]
    med <- active_map$mediator[i]
    dn <- active_map$design_name[i]
    if (dn %in% names(pvals)) out[[view]][med] <- as.numeric(pvals[dn])
  }

  out
}

pick_lambda_from_cv <- function(cvfit, lambda_choice = c("lambda.1se", "lambda.min")) {
  lambda_choice <- match.arg(lambda_choice)
  lam <- suppressWarnings(as.numeric(cvfit[[lambda_choice]]))
  if (!is.finite(lam) || lam <= 0) {
    lam_grid <- suppressWarnings(as.numeric(cvfit$lambda))
    lam <- if (length(lam_grid) > 0) max(min(lam_grid), 1e-06) else 1e-04
  }
  lam
}

debiased_lasso_pvals_multiview_gaussian <- function(
  Y,
  X,
  blocks,
  C = NULL,
  coef_by_view,
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_targets = 200L
) {
  # HIMA-style de-biased lasso for the B-stage.
  #
  # This uses the whole SIS-selected multiview design, then de-biases the
  # nonzero lasso targets. It is more defensible than active-set refit because
  # it corrects lasso shrinkage, but it is still not a complete selective
  # inference theorem for the full screening + fusion pipeline.
  lambda_choice <- match.arg(lambda_choice)

  out <- lapply(blocks, function(M) setNames(rep(NA_real_, ncol(M)), colnames(M)))

  active_map <- data.frame(view = character(), mediator = character(), design_name = character())
  for (view in names(blocks)) {
    b <- coef_by_view[[view]]
    active <- names(b)[is.finite(b) & b != 0]
    if (length(active) == 0) next
    active_map <- rbind(
      active_map,
      data.frame(
        view = view,
        mediator = active,
        design_name = paste0(view, "::", active),
        stringsAsFactors = FALSE
      )
    )
  }
  if (nrow(active_map) == 0) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_no_active_targets"))
  }

  pref <- prefix_multiview_blocks(blocks)
  dat <- data.frame(y = as.numeric(Y), X = as.numeric(X), check.names = FALSE)
  x_cols <- "X"
  if (!is.null(C) && ncol(C) > 0) {
    C_df <- as.data.frame(C, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(C_df) <- paste0("C::", colnames(C_df))
    dat <- cbind(dat, C_df)
    x_cols <- c(x_cols, colnames(C_df))
  }
  dat <- cbind(dat, as.data.frame(pref$X, check.names = FALSE, stringsAsFactors = FALSE))

  cc <- complete.cases(dat)
  dat <- dat[cc, , drop = FALSE]
  if (nrow(dat) < 20) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_too_few_samples"))
  }

  y <- as.numeric(dat$y)
  Z_raw <- as.matrix(dat[, setdiff(colnames(dat), "y"), drop = FALSE])
  n <- nrow(Z_raw)
  p <- ncol(Z_raw)
  if (p < 2) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_design_too_small"))
  }

  z_center <- colMeans(Z_raw)
  z_scale <- apply(Z_raw, 2, stats::sd)
  z_scale[!is.finite(z_scale) | z_scale < 1e-08] <- 1
  Z <- sweep(sweep(Z_raw, 2, z_center, "-"), 2, z_scale, "/")
  y0 <- y - mean(y)

  pf <- rep(1, p)
  names(pf) <- colnames(Z)
  pf[colnames(Z) %in% x_cols] <- 0

  cv_main <- glmnet::cv.glmnet(
    x = Z,
    y = y0,
    family = "gaussian",
    alpha = 1,
    penalty.factor = pf,
    standardize = FALSE,
    intercept = FALSE
  )
  lam_main <- pick_lambda_from_cv(cv_main, lambda_choice = lambda_choice)
  fit_main <- glmnet::glmnet(
    x = Z,
    y = y0,
    family = "gaussian",
    alpha = 1,
    lambda = lam_main,
    penalty.factor = pf,
    standardize = FALSE,
    intercept = FALSE
  )
  beta_hat <- as.numeric(glmnet::coef.glmnet(fit_main, s = lam_main))[-1]
  names(beta_hat) <- colnames(Z)

  resid <- y0 - as.numeric(Z %*% beta_hat)
  df_eff <- sum(abs(beta_hat) > 0) + 1L
  sigma2_hat <- sum(resid^2) / max(1, n - df_eff)
  sigma2_hat <- max(sigma2_hat, 1e-10)

  target_names <- active_map$design_name[active_map$design_name %in% colnames(Z)]
  if (length(target_names) == 0) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_no_valid_targets"))
  }

  max_targets <- max(1L, as.integer(max_targets))
  if (length(target_names) > max_targets) {
    ord <- order(abs(beta_hat[target_names]), decreasing = TRUE, na.last = TRUE)
    target_names <- target_names[ord[seq_len(max_targets)]]
  }

  p_map <- setNames(rep(NA_real_, length(target_names)), target_names)
  for (nm in target_names) {
    j <- match(nm, colnames(Z))
    if (is.na(j)) next

    z_j <- Z[, j]
    z_m <- Z[, -j, drop = FALSE]
    if (ncol(z_m) == 0) next

    cv_j <- glmnet::cv.glmnet(
      x = z_m,
      y = z_j,
      family = "gaussian",
      alpha = 1,
      standardize = FALSE,
      intercept = FALSE
    )
    lam_j <- pick_lambda_from_cv(cv_j, lambda_choice = lambda_choice)
    fit_j <- glmnet::glmnet(
      x = z_m,
      y = z_j,
      family = "gaussian",
      alpha = 1,
      lambda = lam_j,
      standardize = FALSE,
      intercept = FALSE
    )
    gamma_j <- as.numeric(glmnet::coef.glmnet(fit_j, s = lam_j))[-1]
    r_j <- z_j - as.numeric(z_m %*% gamma_j)
    tau_j <- mean(z_j * r_j)
    if (!is.finite(tau_j) || abs(tau_j) < 1e-08) next

    w_j <- r_j / tau_j
    beta_db <- beta_hat[[nm]] + mean(w_j * resid)
    se_db <- sqrt(sigma2_hat * mean(w_j^2) / n)
    if (!is.finite(se_db) || se_db <= 0) next

    p_map[[nm]] <- 2 * stats::pnorm(-abs(beta_db / se_db))
  }

  for (i in seq_len(nrow(active_map))) {
    view <- active_map$view[i]
    med <- active_map$mediator[i]
    dn <- active_map$design_name[i]
    if (dn %in% names(p_map)) out[[view]][med] <- as.numeric(p_map[[dn]])
  }

  list(p_b = out, method = "debiased_lasso_gaussian")
}

infer_p_b_multiview <- function(
  Y,
  X,
  blocks,
  C = NULL,
  coef_by_view,
  y_family = c("gaussian", "binomial"),
  method = c("debiased_lasso", "refit"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_debias_targets = 200L
) {
  y_family <- match.arg(y_family)
  method <- match.arg(method)
  lambda_choice <- match.arg(lambda_choice)

  if (identical(method, "debiased_lasso") && identical(y_family, "gaussian")) {
    db_try <- try(
      debiased_lasso_pvals_multiview_gaussian(
        Y = Y,
        X = X,
        blocks = blocks,
        C = C,
        coef_by_view = coef_by_view,
        lambda_choice = lambda_choice,
        max_targets = max_debias_targets
      ),
      silent = TRUE
    )
    if (!inherits(db_try, "try-error")) return(db_try)

    msg <- conditionMessage(attr(db_try, "condition"))
    warning("Debiased-lasso B-path inference failed; falling back to refit p-values. Reason: ", msg)
  } else if (identical(method, "debiased_lasso") && !identical(y_family, "gaussian")) {
    warning("debiased_lasso currently supports gaussian Y only; falling back to refit p-values.")
  }

  list(
    p_b = infer_p_b_multiview_refit(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      coef_by_view = coef_by_view,
      y_family = y_family
    ),
    method = "active_set_refit"
  )
}

assemble_mediator_table <- function(screen_tab, b_vec, p_b_vec) {
  a_map <- setNames(screen_tab$a_hat, screen_tab$mediator)
  p_a_map <- setNames(screen_tab$p_a, screen_tab$mediator)
  q_a_map <- setNames(screen_tab$q_a, screen_tab$mediator)
  sel_map <- setNames(screen_tab$selected, screen_tab$mediator)
  med <- names(a_map)
  selected <- as.logical(sel_map[med])

  out <- data.frame(
    mediator = med,
    a = a_map[med],
    p_a = p_a_map[med],
    q_a = q_a_map[med],
    b = b_vec[med],
    p_b = p_b_vec[med],
    selected_by_screen = selected,
    selected_by_sis = selected,
    stringsAsFactors = FALSE
  )

  out$score <- out$a * out$b
  out$abs_score <- abs(out$score)
  both_legs <- is.finite(out$p_a) & is.finite(out$p_b)
  out$joint_p_ab <- NA_real_
  out$joint_p_ab[both_legs] <- pmax(out$p_a[both_legs], out$p_b[both_legs])
  out$p_primary <- out$joint_p_ab
  out$q_primary <- p.adjust(out$p_primary, method = "BH")

  out <- out[order(out$abs_score, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}

estimate_reduced_x_effect <- function(Y, X, C = NULL, y_family = c("gaussian", "binomial")) {
  y_family <- match.arg(y_family)

  W <- cbind(`(Intercept)` = 1, X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) W <- cbind(W, as.matrix(C))
  W <- as.matrix(W)

  fit <- try(
    if (identical(y_family, "gaussian")) {
      lm.fit(x = W, y = as.numeric(Y))
    } else {
      suppressWarnings(glm.fit(x = W, y = as.numeric(Y), family = binomial(), intercept = FALSE))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NA_real_)

  out <- as.numeric(fit$coefficients["X"])
  if (!is.finite(out)) NA_real_ else out
}

compute_effect_decomposition_multiview <- function(combined_mediators, direct_effect, Y, X, C = NULL, y_family = c("gaussian", "binomial")) {
  y_family <- match.arg(y_family)
  total_reduced <- estimate_reduced_x_effect(Y = Y, X = X, C = C, y_family = y_family)

  active <- if (!is.null(combined_mediators) && nrow(combined_mediators) > 0) {
    is.finite(combined_mediators$score) & is.finite(combined_mediators$b) & combined_mediators$b != 0
  } else {
    logical(0)
  }

  out <- data.frame(
    direct_effect = as.numeric(direct_effect),
    indirect_total = if (length(active) > 0) sum(combined_mediators$score[active], na.rm = TRUE) else 0,
    total_reduced = total_reduced,
    effect_scale = if (identical(y_family, "binomial")) "log_odds" else "linear_outcome",
    stringsAsFactors = FALSE
  )
  out$total_decomp <- out$direct_effect + out$indirect_total
  out$prop_mediated <- if (is.finite(total_reduced) && abs(total_reduced) > 1e-08) out$indirect_total / total_reduced else NA_real_
  out$decomposition_gap <- if (is.finite(total_reduced)) total_reduced - out$total_decomp else NA_real_
  out$n_active_total <- if (length(active) > 0) sum(active) else 0L

  views <- if (!is.null(combined_mediators) && nrow(combined_mediators) > 0) unique(as.character(combined_mediators$omics)) else character(0)
  for (view in views) {
    safe <- make.names(view)
    idx <- active & combined_mediators$omics == view
    out[[paste0("indirect_", safe)]] <- sum(combined_mediators$score[idx], na.rm = TRUE)
    out[[paste0("n_active_", safe)]] <- sum(idx)
  }

  out
}

make_bootstrap_sample <- function(blocks, pheno_df, bootstrap_id = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(pheno_df)
  if (n < 2) stop("Need at least two samples for bootstrap.")

  if (!is.null(bootstrap_id)) {
    if (!is.character(bootstrap_id) || length(bootstrap_id) != 1L) {
      stop("bootstrap_id must be NULL or a single column name in pheno_df.")
    }
    if (!(bootstrap_id %in% colnames(pheno_df))) {
      stop("bootstrap_id column not found in pheno_df: ", bootstrap_id)
    }
    unit <- as.character(pheno_df[[bootstrap_id]])
    unit[is.na(unit) | !nzchar(unit)] <- paste0("missing_unit_", seq_len(n))[is.na(unit) | !nzchar(unit)]
    levels_u <- unique(unit)
    sampled_units <- sample(levels_u, size = length(levels_u), replace = TRUE)
    idx <- unlist(lapply(sampled_units, function(u) which(unit == u)), use.names = FALSE)
  } else {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
  }

  new_ids <- paste0(rownames(pheno_df)[idx], "__boot", seq_along(idx))
  pheno_b <- pheno_df[idx, , drop = FALSE]
  rownames(pheno_b) <- new_ids

  blocks_b <- lapply(blocks, function(M) {
    out <- M[idx, , drop = FALSE]
    rownames(out) <- new_ids
    out
  })
  names(blocks_b) <- names(blocks)

  list(blocks = blocks_b, pheno_df = pheno_b)
}

extract_bootstrap_vectors <- function(fit, key_ref) {
  tab <- fit$combined_mediators
  key <- paste(tab$omics, tab$mediator, sep = "::")
  score <- setNames(rep(NA_real_, length(key_ref)), key_ref)
  active <- setNames(rep(NA_real_, length(key_ref)), key_ref)

  hit <- match(key_ref, key)
  ok <- !is.na(hit)
  if (any(ok)) {
    score[ok] <- tab$score[hit[ok]]
    active[ok] <- as.numeric(is.finite(tab$b[hit[ok]]) & tab$b[hit[ok]] != 0)
  }

  list(score = score, active = active)
}

add_bootstrap_columns <- function(combined, score_mat, active_mat, ci_level = 0.95) {
  if (is.null(score_mat) || nrow(score_mat) == 0 || nrow(combined) == 0) return(combined)

  alpha <- 1 - ci_level
  combined$score_boot_mean <- apply(score_mat, 2, mean, na.rm = TRUE)
  combined$score_boot_sd <- apply(score_mat, 2, stats::sd, na.rm = TRUE)
  combined$score_boot_low <- apply(score_mat, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
  combined$score_boot_high <- apply(score_mat, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  combined$score_boot_selection_freq <- apply(active_mat, 2, mean, na.rm = TRUE)
  combined
}

add_effect_bootstrap_columns <- function(effect_decomposition, effect_boot, ci_level = 0.95) {
  if (is.null(effect_boot) || nrow(effect_boot) == 0 || nrow(effect_decomposition) == 0) {
    return(effect_decomposition)
  }

  alpha <- 1 - ci_level
  numeric_cols <- colnames(effect_boot)[vapply(effect_boot, is.numeric, logical(1))]
  for (nm in numeric_cols) {
    vals <- effect_boot[[nm]]
    effect_decomposition[[paste0(nm, "_boot_mean")]] <- mean(vals, na.rm = TRUE)
    effect_decomposition[[paste0(nm, "_boot_sd")]] <- stats::sd(vals, na.rm = TRUE)
    effect_decomposition[[paste0(nm, "_boot_low")]] <- stats::quantile(vals, probs = alpha / 2, na.rm = TRUE, names = FALSE)
    effect_decomposition[[paste0(nm, "_boot_high")]] <- stats::quantile(vals, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  }

  effect_decomposition
}

bootstrap_multiview_fit <- function(
  blocks,
  pheno_df,
  x_var,
  y_var,
  covariates,
  residualize,
  sis_n,
  sis_rank,
  screen_method,
  fusion_mode,
  y_family,
  lambda_choice,
  b_inference,
  debias_max_targets,
  coop_rho,
  coop_maxit,
  coop_tol,
  key_ref,
  bootstrap_repeats,
  bootstrap_ci_level,
  bootstrap_seed,
  bootstrap_id
) {
  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  if (bootstrap_repeats < 1L) return(NULL)

  score_mat <- matrix(NA_real_, nrow = bootstrap_repeats, ncol = length(key_ref))
  active_mat <- matrix(NA_real_, nrow = bootstrap_repeats, ncol = length(key_ref))
  colnames(score_mat) <- key_ref
  colnames(active_mat) <- key_ref
  effect_rows <- list()
  failures <- character(0)

  for (b in seq_len(bootstrap_repeats)) {
    bt <- try(
      make_bootstrap_sample(
        blocks = blocks,
        pheno_df = pheno_df,
        bootstrap_id = bootstrap_id,
        seed = bootstrap_seed + b
      ),
      silent = TRUE
    )
    if (inherits(bt, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(bt, "condition"))))
      next
    }

    fit_b <- try(
      fit_multiview_parallel_zentangler_blocks(
        blocks = bt$blocks,
        pheno_df = bt$pheno_df,
        x_var = x_var,
        y_var = y_var,
        covariates = covariates,
        residualize = residualize,
        sis_n = sis_n,
        sis_rank = sis_rank,
        screen_method = screen_method,
        fusion_mode = fusion_mode,
        y_family = y_family,
        lambda_choice = lambda_choice,
        b_inference = b_inference,
        debias_max_targets = debias_max_targets,
        coop_rho = coop_rho,
        coop_maxit = coop_maxit,
        coop_tol = coop_tol,
        seed = bootstrap_seed + b,
        return_fits = FALSE,
        bootstrap_repeats = 0L
      ),
      silent = TRUE
    )
    if (inherits(fit_b, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(fit_b, "condition"))))
      next
    }

    vec <- extract_bootstrap_vectors(fit_b, key_ref = key_ref)
    score_mat[b, ] <- vec$score
    active_mat[b, ] <- vec$active
    effect_rows[[length(effect_rows) + 1L]] <- fit_b$effect_decomposition
  }

  effect_boot <- if (length(effect_rows) > 0) {
    out <- do.call(rbind, effect_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame()
  }

  list(
    score_matrix = score_mat,
    active_matrix = active_mat,
    effect_decomposition = effect_boot,
    ci_level = bootstrap_ci_level,
    repeats_requested = bootstrap_repeats,
    repeats_successful = length(effect_rows),
    failures = failures
  )
}

# -----------------------------------------------------------------------------
# Public K-view fit
# -----------------------------------------------------------------------------

fit_multiview_parallel_zentangler_blocks <- function(
  blocks,
  pheno_df,
  x_var,
  y_var,
  covariates = NULL,
  residualize = FALSE,
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  fusion_mode = c("early", "intermediate", "late"),
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  b_inference = c("debiased_lasso", "refit"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04,
  seed = 1,
  return_fits = FALSE,
  bootstrap_repeats = 0L,
  bootstrap_ci_level = 0.95,
  bootstrap_seed = seed + 10000L,
  bootstrap_id = NULL
) {
  sis_rank <- match.arg(sis_rank)
  screen_method <- match.arg(screen_method)
  fusion_mode <- match.arg(fusion_mode)
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  b_inference <- match.arg(b_inference)
  set.seed(seed)

  needed <- unique(c(x_var, y_var, covariates, bootstrap_id))
  al <- align_samples_multiview_blocks(blocks, pheno_df, needed_pheno_cols = needed)
  blocks0 <- al$blocks
  X_raw <- al$pheno[[x_var]]
  Y_raw <- al$pheno[[y_var]]
  if (!is.numeric(X_raw) || !is.numeric(Y_raw)) stop("x_var and y_var must be numeric.")

  C <- make_covariate_matrix(al$pheno, covariates)
  blocks_model <- blocks0
  X_model <- X_raw
  Y_model <- Y_raw

  if (residualize && !is.null(C) && ncol(C) > 0) {
    blocks_model <- lapply(blocks0, residualize_matrix, C = C)
    X_model <- as.numeric(residualize_matrix(matrix(X_raw, ncol = 1), C))
    Y_model <- as.numeric(residualize_matrix(matrix(Y_raw, ncol = 1), C))
    C_model <- NULL
  } else {
    C_model <- C
  }

  screens <- lapply(names(blocks_model), function(view) {
    hima_screen_mediators(
      blocks_model[[view]],
      X_model,
      C = C_model,
      sis_n = sis_n,
      rank_by = sis_rank,
      screen_method = screen_method
    )
  })
  names(screens) <- names(blocks_model)

  selected_blocks <- lapply(names(blocks_model), function(view) {
    blocks_model[[view]][, screens[[view]]$selected, drop = FALSE]
  })
  names(selected_blocks) <- names(blocks_model)

  yfit <- fit_y_multiview_stage(
    Y = Y_model,
    X = X_model,
    blocks = selected_blocks,
    C = C_model,
    fusion_mode = fusion_mode,
    y_family = y_family,
    lambda_choice = lambda_choice,
    coop_rho = coop_rho,
    coop_maxit = coop_maxit,
    coop_tol = coop_tol
  )

  p_b_infer <- infer_p_b_multiview(
    Y = Y_model,
    X = X_model,
    blocks = selected_blocks,
    C = C_model,
    coef_by_view = yfit$coef_by_view,
    y_family = y_family,
    method = b_inference,
    lambda_choice = lambda_choice,
    max_debias_targets = debias_max_targets
  )
  p_b_selected <- p_b_infer$p_b

  view_tables <- lapply(names(blocks_model), function(view) {
    b_full <- setNames(rep(0, ncol(blocks_model[[view]])), colnames(blocks_model[[view]]))
    hit_b <- intersect(names(yfit$coef_by_view[[view]]), names(b_full))
    b_full[hit_b] <- yfit$coef_by_view[[view]][hit_b]

    p_full <- setNames(rep(NA_real_, ncol(blocks_model[[view]])), colnames(blocks_model[[view]]))
    hit_p <- intersect(names(p_b_selected[[view]]), names(p_full))
    p_full[hit_p] <- p_b_selected[[view]][hit_p]

    assemble_mediator_table(screens[[view]]$table, b_vec = b_full, p_b_vec = p_full)
  })
  names(view_tables) <- names(blocks_model)

  combined <- do.call(
    rbind,
    lapply(names(view_tables), function(view) {
      data.frame(omics = view, view_tables[[view]], stringsAsFactors = FALSE)
    })
  )
  combined$p_primary <- combined$joint_p_ab
  combined$q_primary <- p.adjust(combined$p_primary, method = "BH")
  combined <- combined[order(combined$abs_score, decreasing = TRUE), , drop = FALSE]
  rownames(combined) <- NULL

  effect_decomposition <- compute_effect_decomposition_multiview(
    combined_mediators = combined,
    direct_effect = yfit$x_coef,
    Y = Y_model,
    X = X_model,
    C = C_model,
    y_family = y_family
  )

  key_ref <- paste(combined$omics, combined$mediator, sep = "::")
  bootstrap <- NULL
  if (as.integer(bootstrap_repeats) > 0L) {
    bootstrap <- bootstrap_multiview_fit(
      blocks = blocks0,
      pheno_df = al$pheno,
      x_var = x_var,
      y_var = y_var,
      covariates = covariates,
      residualize = residualize,
      sis_n = sis_n,
      sis_rank = sis_rank,
      screen_method = screen_method,
      fusion_mode = fusion_mode,
      y_family = y_family,
      lambda_choice = lambda_choice,
      b_inference = b_inference,
      debias_max_targets = debias_max_targets,
      coop_rho = coop_rho,
      coop_maxit = coop_maxit,
      coop_tol = coop_tol,
      key_ref = key_ref,
      bootstrap_repeats = bootstrap_repeats,
      bootstrap_ci_level = bootstrap_ci_level,
      bootstrap_seed = bootstrap_seed,
      bootstrap_id = bootstrap_id
    )
    combined <- add_bootstrap_columns(
      combined = combined,
      score_mat = bootstrap$score_matrix,
      active_mat = bootstrap$active_matrix,
      ci_level = bootstrap_ci_level
    )
    effect_decomposition <- add_effect_bootstrap_columns(
      effect_decomposition = effect_decomposition,
      effect_boot = bootstrap$effect_decomposition,
      ci_level = bootstrap_ci_level
    )
  }

  views_out <- lapply(names(blocks_model), function(view) {
    list(
      screened_mediators = screens[[view]]$table,
      selected_mediators = screens[[view]]$selected,
      a = setNames(screens[[view]]$table$a_hat, screens[[view]]$table$mediator),
      p_a = setNames(screens[[view]]$table$p_a, screens[[view]]$table$mediator),
      b = setNames(view_tables[[view]]$b, view_tables[[view]]$mediator),
      p_b = setNames(view_tables[[view]]$p_b, view_tables[[view]]$mediator),
      mediators = view_tables[[view]]
    )
  })
  names(views_out) <- names(blocks_model)

  list(
    settings = list(
      x_var = x_var,
      y_var = y_var,
      view_names = names(blocks_model),
      n_views = length(blocks_model),
      covariates = covariates,
      residualize = residualize,
      sis_n = sis_n,
      sis_rank = sis_rank,
      screen_method = screen_method,
      a_stage_model = "lm_univariate",
      fusion_mode = fusion_mode,
      y_family = y_family,
      lambda_choice = lambda_choice,
      b_inference = b_inference,
      b_inference_method = p_b_infer$method,
      debias_max_targets = debias_max_targets,
      y_model_method = yfit$method,
      bootstrap_repeats = as.integer(bootstrap_repeats),
      bootstrap_ci_level = bootstrap_ci_level,
      bootstrap_seed = bootstrap_seed,
      bootstrap_id = bootstrap_id
    ),
    sample_ids = rownames(al$pheno),
    views = views_out,
    combined_mediators = combined,
    x_to_y_coef = yfit$x_coef,
    effect_decomposition = effect_decomposition,
    bootstrap = bootstrap,
    fits = if (return_fits) yfit else NULL
  )
}

fit_multiview_parallel_zentangler <- function(
  mae,
  x_var,
  y_var,
  view_names = NULL,
  assay_names = NULL,
  covariates = NULL,
  residualize = FALSE,
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  fusion_mode = c("early", "intermediate", "late"),
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  b_inference = c("debiased_lasso", "refit"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04,
  seed = 1,
  return_fits = FALSE,
  bootstrap_repeats = 0L,
  bootstrap_ci_level = 0.95,
  bootstrap_seed = seed + 10000L,
  bootstrap_id = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  # MAE-first public entry point.
  #
  # This is intentionally thin: all Bioconductor sample/assay bookkeeping is
  # handled here, then the internal matrix engine does the statistical work.
  # That keeps the package-facing API Bioconductor-native without duplicating
  # the modeling code.
  duplicate_primary <- match.arg(duplicate_primary)
  screen_method <- match.arg(screen_method)

  inputs <- zentangler_mae_to_blocks(
    mae = mae,
    view_names = view_names,
    assay_names = assay_names,
    duplicate_primary = duplicate_primary,
    fill_nonfinite_zero = fill_nonfinite_zero
  )

  fit <- fit_multiview_parallel_zentangler_blocks(
    blocks = inputs$blocks,
    pheno_df = inputs$pheno_df,
    x_var = x_var,
    y_var = y_var,
    covariates = covariates,
    residualize = residualize,
    sis_n = sis_n,
    sis_rank = sis_rank,
    screen_method = screen_method,
    fusion_mode = fusion_mode,
    y_family = y_family,
    lambda_choice = lambda_choice,
    b_inference = b_inference,
    debias_max_targets = debias_max_targets,
    coop_rho = coop_rho,
    coop_maxit = coop_maxit,
    coop_tol = coop_tol,
    seed = seed,
    return_fits = return_fits,
    bootstrap_repeats = bootstrap_repeats,
    bootstrap_ci_level = bootstrap_ci_level,
    bootstrap_seed = bootstrap_seed,
    bootstrap_id = bootstrap_id
  )

  fit$settings$input_container <- "MultiAssayExperiment"
  fit$settings$mae_view_map <- inputs$view_map
  fit
}
