# -----------------------------------------------------------------------------
# K-layer sequential mediation
# -----------------------------------------------------------------------------

zentangler_validate_stage_views <- function(stage_views, available_views) {
  if (is.null(stage_views)) {
    stage_views <- as.list(available_views)
    names(stage_views) <- paste0("stage", seq_along(stage_views))
  }
  if (!is.list(stage_views) || length(stage_views) < 2L) {
    stop("stage_views must be a list with at least two ordered mediator layers.", call. = FALSE)
  }

  stage_views <- lapply(stage_views, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) stop("Each stage_views entry must contain at least one view.", call. = FALSE)
    x
  })
  if (is.null(names(stage_views)) || any(!nzchar(names(stage_views)))) {
    names(stage_views) <- paste0("stage", seq_along(stage_views))
  }

  missing <- setdiff(unique(unlist(stage_views, use.names = FALSE)), available_views)
  if (length(missing) > 0L) {
    stop("stage_views contains views not present in the input: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  duplicated_views <- unique(unlist(stage_views, use.names = FALSE)[
    duplicated(unlist(stage_views, use.names = FALSE))
  ])
  if (length(duplicated_views) > 0L) {
    stop("Views can appear in only one sequential layer: ", paste(duplicated_views, collapse = ", "), call. = FALSE)
  }

  stage_views
}

zentangler_validate_path_templates <- function(path_templates, available_views) {
  if (is.null(path_templates)) return(NULL)
  if (!is.list(path_templates) || length(path_templates) == 0L) {
    stop("path_templates must be a non-empty list of character vectors.", call. = FALSE)
  }
  if (is.null(names(path_templates)) || any(!nzchar(names(path_templates)))) {
    names(path_templates) <- paste0("path", seq_along(path_templates))
  }
  path_templates <- lapply(path_templates, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) stop("Each path template must contain at least one view.", call. = FALSE)
    if (length(x) < 2L) {
      stop(
        "Sequential path templates must contain at least two views. ",
        "Use fit_multiview_parallel_zentangler() for pure parallel mediation.",
        call. = FALSE
      )
    }
    missing <- setdiff(x, available_views)
    if (length(missing) > 0L) {
      stop("path_templates contains views not present in the input: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    x
  })
  path_templates
}

zentangler_map_stage_views <- function(stage_views, view_map) {
  if (is.null(stage_views) || is.null(view_map) || nrow(view_map) == 0L) return(stage_views)
  original_to_safe <- setNames(as.character(view_map$view), as.character(view_map$experiment_name))
  safe_to_safe <- setNames(as.character(view_map$view), as.character(view_map$view))
  lapply(stage_views, function(x) {
    x <- as.character(x)
    mapped <- original_to_safe[x]
    missing <- is.na(mapped)
    mapped[missing] <- safe_to_safe[x[missing]]
    mapped
  })
}

zentangler_map_path_templates <- function(path_templates, view_map) {
  if (is.null(path_templates) || is.null(view_map) || nrow(view_map) == 0L) return(path_templates)
  original_to_safe <- setNames(as.character(view_map$view), as.character(view_map$experiment_name))
  safe_to_safe <- setNames(as.character(view_map$view), as.character(view_map$view))
  lapply(path_templates, function(x) {
    x <- as.character(x)
    mapped <- original_to_safe[x]
    missing <- is.na(mapped)
    mapped[missing] <- safe_to_safe[x[missing]]
    mapped
  })
}

zentangler_check_route_api <- function(path_templates) {
  if (is.null(path_templates)) return(invisible(NULL))
  route_lengths <- vapply(path_templates, length, integer(1))
  if (length(route_lengths) != 1L) {
    stop(
      "fit_sequential_zentangler() accepts exactly one sequential route. ",
      "Run multiple routes in separate calls.",
      call. = FALSE
    )
  }
  if (route_lengths[[1L]] < 2L) {
    stop(
      "fit_sequential_zentangler() requires at least two modalities in its route. ",
      "Use fit_multiview_parallel_zentangler() for pure parallel mediation.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

zentangler_route_set_type <- function(path_templates = NULL, stage_views = NULL) {
  if (!is.null(path_templates)) {
    route_lengths <- vapply(path_templates, length, integer(1))
    if (length(route_lengths) == 1L && route_lengths[[1L]] > 1L) return("sequential")
    return("sequential")
  }
  if (!is.null(stage_views)) {
    if (length(stage_views) == 1L) return("parallel")
    return("sequential")
  }
  NA_character_
}

zentangler_feature_key <- function(view, mediator) {
  paste(view, mediator, sep = "::")
}

zentangler_residualize_vector <- function(y, Z = NULL) {
  y <- as.numeric(y)
  if (is.null(Z) || ncol(Z) == 0L) return(y)
  fit <- lm.fit(x = cbind(`(Intercept)` = 1, as.matrix(Z)), y = y)
  as.numeric(fit$residuals)
}

zentangler_cor_pvalue <- function(r, n, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (!is.finite(r) || n < 4L || abs(r) >= 1) {
    if (is.finite(r) && abs(r) >= 1) return(0)
    return(NA_real_)
  }
  if (identical(method, "pearson")) {
    stat <- r * sqrt((n - 2) / max(1e-12, 1 - r^2))
    2 * stats::pt(abs(stat), df = n - 2, lower.tail = FALSE)
  } else {
    # Large-sample approximation used for speed in high-dimensional screens.
    stat <- r * sqrt((n - 2) / max(1e-12, 1 - r^2))
    2 * stats::pt(abs(stat), df = n - 2, lower.tail = FALSE)
  }
}

zentangler_standardize_for_cor <- function(M, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  M <- as.matrix(M)
  dn <- dimnames(M)
  if (identical(method, "spearman")) {
    M <- apply(M, 2, rank, ties.method = "average", na.last = "keep")
    M <- as.matrix(M)
    dimnames(M) <- dn
  }
  center <- colMeans(M, na.rm = TRUE)
  scale <- apply(M, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- NA_real_
  out <- sweep(sweep(M, 2, center, "-"), 2, scale, "/")
  dimnames(out) <- dn
  out
}

zentangler_screen_cor_links <- function(
  from_df,
  to_df,
  X,
  C = NULL,
  cor_method = c("spearman", "pearson"),
  cor_fdr_method = c("BH", "BY"),
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25,
  residualize_links = TRUE,
  max_links = 5000L
) {
  cor_method <- match.arg(cor_method)
  cor_fdr_method <- validate_fdr_method(cor_fdr_method)
  max_links <- max(1L, as.integer(max_links))

  if (nrow(from_df) != nrow(to_df)) stop("from_df and to_df must have the same samples.", call. = FALSE)
  if (ncol(from_df) == 0L || ncol(to_df) == 0L) return(data.frame())

  Z <- NULL
  if (isTRUE(residualize_links)) {
    Z <- cbind(X = as.numeric(X))
    if (!is.null(C) && ncol(C) > 0L) Z <- cbind(Z, as.matrix(C))
  }

  A <- as.matrix(from_df)
  B <- as.matrix(to_df)
  if (!is.null(Z)) {
    a_dn <- dimnames(A)
    b_dn <- dimnames(B)
    A <- apply(A, 2, zentangler_residualize_vector, Z = Z)
    B <- apply(B, 2, zentangler_residualize_vector, Z = Z)
    A <- as.matrix(A)
    B <- as.matrix(B)
    dimnames(A) <- a_dn
    dimnames(B) <- b_dn
  }

  A <- zentangler_standardize_for_cor(A, method = cor_method)
  B <- zentangler_standardize_for_cor(B, method = cor_method)
  keep_a <- colSums(!is.finite(A)) == 0L
  keep_b <- colSums(!is.finite(B)) == 0L
  A <- A[, keep_a, drop = FALSE]
  B <- B[, keep_b, drop = FALSE]
  if (ncol(A) == 0L || ncol(B) == 0L) return(data.frame())

  R <- crossprod(A, B) / max(1L, nrow(A) - 1L)
  idx <- which(abs(R) >= min_abs_cor, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(data.frame())

  cor_val <- R[idx]
  p_val <- vapply(cor_val, zentangler_cor_pvalue, numeric(1), n = nrow(A), method = cor_method)
  q_val <- stats::p.adjust(p_val, method = cor_fdr_method)
  out <- data.frame(
    from_key = colnames(A)[idx[, 1]],
    to_key = colnames(B)[idx[, 2]],
    cor = as.numeric(cor_val),
    abs_cor = abs(as.numeric(cor_val)),
    p_cor = as.numeric(p_val),
    q_cor = as.numeric(q_val),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$q_cor) & out$q_cor <= cor_q_threshold, , drop = FALSE]
  if (nrow(out) == 0L) return(out)
  out <- out[order(out$q_cor, -out$abs_cor), , drop = FALSE]
  out <- out[seq_len(min(max_links, nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

zentangler_estimate_link_effects <- function(links, data_by_key, X, C = NULL) {
  if (is.null(links) || nrow(links) == 0L) return(links)
  d_hat <- p_d <- rep(NA_real_, nrow(links))
  for (i in seq_len(nrow(links))) {
    from <- data_by_key[[links$from_key[i]]]
    to <- data_by_key[[links$to_key[i]]]
    dat <- data.frame(to = as.numeric(to), X = as.numeric(X), from = as.numeric(from), check.names = FALSE)
    if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
    fit <- try(stats::lm(to ~ ., data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) next
    co <- summary(fit)$coefficients
    if ("from" %in% rownames(co)) {
      d_hat[i] <- co["from", "Estimate"]
      p_d[i] <- co["from", "Pr(>|t|)"]
    }
  }
  links$d_hat <- d_hat
  links$p_d <- p_d
  links
}

zentangler_path_primary_p <- function(p_a, p_d = numeric(0), p_b, b = NA_real_) {
  if (!is.finite(b) || b == 0) return(NA_real_)
  vals <- c(p_a, p_d, p_b)
  if (length(vals) == 0L || any(!is.finite(vals))) return(NA_real_)
  max(vals)
}

validate_sequential_path_inference <- function(path_inference) {
  match.arg(path_inference, choices = c("model_based", "bootstrap_score"))
}

zentangler_subset_covariates <- function(C, idx) {
  if (is.null(C) || ncol(C) == 0L) return(NULL)
  C[idx, , drop = FALSE]
}

zentangler_estimate_a_one <- function(mediator, X, C = NULL) {
  dat <- data.frame(mediator = as.numeric(mediator), X = as.numeric(X), check.names = FALSE)
  if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
  fit <- try(stats::lm(mediator ~ ., data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(a = NA_real_, p_a = NA_real_))
  co <- summary(fit)$coefficients
  if (!("X" %in% rownames(co))) return(c(a = NA_real_, p_a = NA_real_))
  c(a = co["X", "Estimate"], p_a = co["X", "Pr(>|t|)"])
}

zentangler_bootstrap_sequential_paths <- function(
  paths,
  data_by_key,
  pheno_df,
  X,
  Y,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  fusion_mode = c("early", "intermediate", "late"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  b_inference = c("debiased_lasso", "debiased_logistic_lasso", "debiased_cox_lasso", "refit", "bootstrap"),
  causal_inference = c("none", "bootstrap"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04,
  repeats = 0L,
  ci_level = 0.95,
  seed = 10001L,
  bootstrap_id = NULL,
  effect_x0 = NULL,
  effect_x1 = NULL
) {
  y_family <- match.arg(y_family)
  fusion_mode <- match.arg(fusion_mode)
  lambda_choice <- match.arg(lambda_choice)
  b_inference <- match.arg(b_inference)
  causal_inference <- match.arg(causal_inference)
  debias_max_targets <- max(1L, as.integer(debias_max_targets))
  repeats <- max(0L, as.integer(repeats))
  if (repeats < 1L || is.null(paths) || nrow(paths) == 0L) return(NULL)
  if (!("key_path" %in% colnames(paths))) stop("paths must contain key_path.", call. = FALSE)

  key_ref <- as.character(paths$key_path)
  score_mat <- matrix(NA_real_, nrow = repeats, ncol = length(key_ref), dimnames = list(NULL, key_ref))
  active_mat <- matrix(NA_real_, nrow = repeats, ncol = length(key_ref), dimnames = list(NULL, key_ref))
  effect_rows <- list()
  failures <- character(0)

  pseudo_blocks <- list(.sequential = matrix(0, nrow = nrow(pheno_df), ncol = 1L, dimnames = list(rownames(pheno_df), ".dummy")))

  for (b in seq_len(repeats)) {
    bt <- try(
      make_bootstrap_sample(
        blocks = pseudo_blocks,
        pheno_df = pheno_df,
        bootstrap_id = bootstrap_id,
        seed = seed + b
      ),
      silent = TRUE
    )
    if (inherits(bt, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(bt, "condition"))))
      next
    }

    original_ids <- sub("__boot[0-9]+$", "", rownames(bt$pheno_df))
    idx <- match(original_ids, rownames(pheno_df))
    if (anyNA(idx)) {
      failures <- c(failures, paste0("bootstrap_", b, ": could not match bootstrap sample IDs."))
      next
    }

    X_b <- X[idx]
    Y_b <- if (inherits(Y, "Surv")) Y[idx, , drop = FALSE] else Y[idx]
    C_b <- zentangler_subset_covariates(C, idx)
    terminal_keys <- unique(as.character(paths$terminal_key %||% zentangler_feature_key(paths$terminal_view, paths$terminal_mediator)))
    terminal_keys <- intersect(terminal_keys, names(data_by_key))
    terminal_df <- as.data.frame(lapply(terminal_keys, function(k) data_by_key[[k]][idx]), check.names = FALSE)
    colnames(terminal_df) <- terminal_keys

    terminal_effects <- try(
      zentangler_fit_terminal_effects(
        Y = Y_b,
        X = X_b,
        terminal_df = terminal_df,
        C = C_b,
        y_family = y_family,
        fusion_mode = fusion_mode,
        lambda_choice = lambda_choice,
        glmnet_alpha = glmnet_alpha,
        b_inference = b_inference,
        debias_max_targets = debias_max_targets,
        coop_rho = coop_rho,
        coop_maxit = coop_maxit,
        coop_tol = coop_tol
      ),
      silent = TRUE
    )
    if (inherits(terminal_effects, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(terminal_effects, "condition"))))
      next
    }
    b_map <- setNames(terminal_effects$b, terminal_effects$terminal_key)

    for (i in seq_len(nrow(paths))) {
      keys <- strsplit(as.character(paths$key_path[i]), " -> ", fixed = TRUE)[[1L]]
      if (length(keys) == 0L || any(!(keys %in% names(data_by_key)))) next
      a_est <- zentangler_estimate_a_one(data_by_key[[keys[[1L]]]][idx], X_b, C_b)
      d_vals <- numeric(0)
      if (length(keys) > 1L) {
        d_vals <- vapply(seq_len(length(keys) - 1L), function(j) {
          dat <- data.frame(
            to = as.numeric(data_by_key[[keys[[j + 1L]]]][idx]),
            X = as.numeric(X_b),
            from = as.numeric(data_by_key[[keys[[j]]]][idx]),
            check.names = FALSE
          )
          if (!is.null(C_b) && ncol(C_b) > 0L) dat <- cbind(dat, as.data.frame(C_b, check.names = FALSE))
          fit <- try(stats::lm(to ~ ., data = dat), silent = TRUE)
          if (inherits(fit, "try-error")) return(NA_real_)
          co <- stats::coef(fit)
          as.numeric(co[["from"]] %||% NA_real_)
        }, numeric(1))
      }
      b_val <- b_map[[tail(keys, 1L)]] %||% NA_real_
      score <- as.numeric(a_est[["a"]]) * prod(d_vals) * as.numeric(b_val)
      score_mat[b, i] <- score
      active_mat[b, i] <- as.numeric(is.finite(b_val) && b_val != 0)
    }
    if (identical(causal_inference, "bootstrap")) {
      data_by_key_b <- lapply(data_by_key, function(x) x[idx])
      paths_b <- paths
      paths_b$sequential_score <- score_mat[b, ]
      term_b <- match(paths_b$terminal_key, terminal_effects$terminal_key)
      if ("terminal_key" %in% colnames(paths_b) && any(!is.na(term_b))) {
        paths_b$b[!is.na(term_b)] <- terminal_effects$b[term_b[!is.na(term_b)]]
        paths_b$p_b[!is.na(term_b)] <- terminal_effects$p_b[term_b[!is.na(term_b)]]
      }
      eff_b <- try(
        compute_effect_decomposition_sequential(
          sequential_paths = paths_b,
          Y = Y_b,
          X = X_b,
          data_by_key = data_by_key_b,
          C = C_b,
          y_family = y_family,
          x0 = effect_x0,
          x1 = effect_x1
        ),
        silent = TRUE
      )
      if (inherits(eff_b, "try-error")) {
        failures <- c(failures, paste0("bootstrap_", b, "_effects: ", conditionMessage(attr(eff_b, "condition"))))
      } else {
        effect_rows[[length(effect_rows) + 1L]] <- eff_b
      }
    }
  }

  effect_boot <- if (length(effect_rows) > 0L) {
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
    ci_level = ci_level,
    repeats_requested = repeats,
    repeats_successful = sum(rowSums(is.finite(score_mat)) > 0L),
    failures = failures
  )
}

zentangler_add_sequential_bootstrap_columns <- function(paths, bootstrap, fdr_method = c("BH", "BY"), use_primary = FALSE) {
  fdr_method <- validate_fdr_method(fdr_method)
  if (is.null(bootstrap) || is.null(bootstrap$score_matrix) || nrow(paths) == 0L) return(paths)
  score_mat <- bootstrap$score_matrix
  active_mat <- bootstrap$active_matrix
  alpha <- 1 - (bootstrap$ci_level %||% 0.95)
  paths$sequential_score_boot_mean <- apply(score_mat, 2, mean, na.rm = TRUE)
  paths$sequential_score_boot_sd <- apply(score_mat, 2, stats::sd, na.rm = TRUE)
  paths$sequential_score_boot_low <- apply(score_mat, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
  paths$sequential_score_boot_high <- apply(score_mat, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  paths$path_boot_selection_freq <- apply(active_mat, 2, mean, na.rm = TRUE)
  paths$p_sequential_score_boot <- bootstrap_score_pvalues(paths$sequential_score, score_mat)
  paths$p_sequential_score_boot_sign <- bootstrap_score_sign_pvalues(score_mat)
  if (isTRUE(use_primary)) {
    ok <- is.finite(paths$p_sequential_score_boot)
    paths$p_primary_model_based <- paths$p_primary
    paths$p_primary[ok] <- paths$p_sequential_score_boot[ok]
    paths$q_primary_global <- stats::p.adjust(paths$p_primary, method = fdr_method)
    paths$q_primary <- paths$q_primary_global
  }
  paths
}

zentangler_fit_terminal_effects <- function(
  Y,
  X,
  terminal_df,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  fusion_mode = c("early", "intermediate", "late"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  b_inference = c("debiased_lasso", "debiased_logistic_lasso", "debiased_cox_lasso", "refit", "bootstrap"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04
) {
  y_family <- match.arg(y_family)
  fusion_mode <- match.arg(fusion_mode)
  lambda_choice <- match.arg(lambda_choice)
  b_inference <- match.arg(b_inference)
  debias_max_targets <- max(1L, as.integer(debias_max_targets))
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)
  terminal_keys <- colnames(terminal_df)
  out <- data.frame(
    terminal_key = terminal_keys,
    b = NA_real_,
    p_b = NA_real_,
    fusion_mode = fusion_mode,
    stringsAsFactors = FALSE
  )
  if (length(terminal_keys) == 0L) return(out)

  key_parts <- strsplit(terminal_keys, "::", fixed = TRUE)
  terminal_views <- vapply(key_parts, `[`, character(1), 1L)
  terminal_mediators <- vapply(key_parts, function(x) paste(x[-1L], collapse = "::"), character(1))
  if (identical(fusion_mode, "intermediate") && length(unique(terminal_views)) < 2L) {
    stop(
      "fusion_mode='intermediate' requires at least two terminal views. ",
      "Use fusion_mode='early' or 'late' for a single terminal modality.",
      call. = FALSE
    )
  }
  blocks <- lapply(unique(terminal_views), function(view) {
    idx <- terminal_views == view
    M <- as.matrix(terminal_df[, idx, drop = FALSE])
    colnames(M) <- terminal_mediators[idx]
    M
  })
  names(blocks) <- unique(terminal_views)
  yfit <- fit_y_multiview_stage(
    Y = Y,
    X = X,
    blocks = blocks,
    C = C,
    fusion_mode = fusion_mode,
    y_family = y_family,
    lambda_choice = lambda_choice,
    glmnet_alpha = glmnet_alpha,
    coop_rho = coop_rho,
    coop_maxit = coop_maxit,
    coop_tol = coop_tol
  )
  p_infer <- infer_p_b_multiview(
    Y = Y,
    X = X,
    blocks = blocks,
    C = C,
    coef_by_view = yfit$coef_by_view,
    y_family = y_family,
    method = b_inference,
    lambda_choice = lambda_choice,
    max_debias_targets = debias_max_targets
  )
  for (view in names(blocks)) {
    for (med in colnames(blocks[[view]])) {
      key <- zentangler_feature_key(view, med)
      ii <- match(key, out$terminal_key)
      if (is.na(ii)) next
      out$b[ii] <- yfit$coef_by_view[[view]][[med]] %||% NA_real_
      out$p_b[ii] <- p_infer$p_b[[view]][[med]] %||% NA_real_
    }
  }
  out$b_inference_method <- p_infer$method %||% b_inference
  if (!any(is.finite(out$b) & out$b != 0)) {
    warning(
      "Terminal fusion selected no mediators, so all terminal b coefficients are zero. ",
      "Consider lambda_choice='lambda.min', a smaller glmnet_alpha, a larger SIS/correlation screen, ",
      "or a different fusion_mode.",
      call. = FALSE
    )
  }
  attr(out, "terminal_yfit") <- yfit
  out
}

zentangler_extend_paths <- function(paths, links, max_paths = 10000L) {
  if (length(paths) == 0L || is.null(links) || nrow(links) == 0L) return(list())
  by_from <- split(seq_len(nrow(links)), links$from_key)
  out <- list()
  ii <- 0L
  for (path in paths) {
    last_key <- tail(path$keys, 1L)
    hit <- by_from[[last_key]]
    if (is.null(hit)) next
    for (j in hit) {
      ii <- ii + 1L
      out[[ii]] <- list(
        keys = c(path$keys, links$to_key[j]),
        link_rows = c(path$link_rows, j)
      )
      if (ii >= max_paths) return(out)
    }
  }
  out
}

#' Fit one k-layer sequential Zentangler mediation route
#'
#' Fits a sequential mediation screen for ordered mediator layers such as
#' `X -> M1 -> M2 -> ... -> Mk -> Y`. The first layer is screened by association
#' with the exposure. Adjacent mediator layers are screened by mediator-to-
#' mediator correlation, then retained links are refit by regression to estimate
#' directed link coefficients.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param x_var Exposure column in `colData(mae)`.
#' @param y_var Outcome column in `colData(mae)` for Gaussian/binomial outcomes.
#' @param stage_views List of ordered mediator layers. Each element is a
#'   character vector of experiment/view names. If `NULL`, `view_names` order is
#'   used as one view per layer.
#' @param path_templates Optional one-route list, for example
#'   `list(route = c("species", "fecal_metabolites"))`. The route must contain
#'   at least two modalities. Run multiple routes in separate calls.
#' @param view_names Optional subset/order of MAE experiments.
#' @param assay_names Optional assay names passed to MAE extraction.
#' @param covariates Optional phenotype covariates.
#' @param y_family Outcome family: `"gaussian"`, `"binomial"`, or `"survival"`.
#' @param survival_time_var,survival_event_var Survival outcome columns.
#' @param sis_n Number of first-layer mediators to retain per first-layer view.
#' @param sis_rank Ranking for first-layer screen.
#' @param screen_method First-layer screen method.
#' @param cor_method Correlation method for adjacent mediator-layer screens.
#' @param min_abs_cor Minimum absolute correlation for adjacent-layer links.
#' @param cor_q_threshold FDR cutoff for adjacent-layer correlation links.
#' @param cor_fdr_method FDR method for adjacent-layer correlation links.
#' @param residualize_links Whether to residualize mediator-to-mediator
#'   correlations on exposure and covariates before screening.
#' @param max_links_per_transition Maximum retained links per adjacent transition.
#' @param max_paths Maximum full paths to keep while expanding across k layers.
#' @param fusion_mode Terminal mediator-to-outcome fusion mode. Uses the same
#'   early, intermediate, and late fusion vocabulary as the parallel API.
#' @param fdr_method FDR method for final path-level p-values.
#' @param fdr_scope Currently `"global"` only for path-level q-values.
#' @param seed Random seed.
#' @param lambda_choice Cross-validated glmnet lambda choice for terminal
#'   fusion.
#' @param glmnet_alpha Glmnet mixing parameter for early and late terminal
#'   fusion. Use 1 for lasso, 0.5 for elastic net, and 0 for ridge.
#' @param b_inference Terminal B-stage inference method. Uses the same choices
#'   as the parallel API: `"debiased_lasso"`, `"debiased_logistic_lasso"`,
#'   `"debiased_cox_lasso"`, `"refit"`, or `"bootstrap"`.
#' @param debias_max_targets Maximum active terminal mediators to use in
#'   debiased B-stage calculations.
#' @param coop_rho Agreement penalty strength for intermediate terminal fusion.
#' @param coop_maxit Maximum cooperative updates for intermediate terminal fusion.
#' @param coop_tol Convergence tolerance for intermediate terminal fusion.
#' @param path_inference Path-level inference. `"model_based"` uses the
#'   conservative refit evidence assembled from A-stage, transition-link, and
#'   terminal B-stage p-values. `"bootstrap_score"` uses bootstrap p-values for
#'   the complete sequential path score and requires `bootstrap_repeats > 0`.
#' @param bootstrap_repeats Number of path bootstrap refits. Use 0 to disable.
#' @param bootstrap_ci_level Bootstrap confidence level.
#' @param bootstrap_seed Bootstrap random seed.
#' @param bootstrap_id Optional phenotype column used for cluster bootstrap.
#' @param causal_inference Optional causal-effect uncertainty mode. Use
#'   `"bootstrap"` to add bootstrap summaries for direct, indirect, and total
#'   sequential effects.
#' @param effect_x0,effect_x1 Optional exposure contrast used for sequential
#'   causal effect summaries.
#' @param duplicate_primary,fill_nonfinite_zero MAE extraction options.
#'
#' @return A list with settings, transition links, terminal effects, and
#'   `sequential_paths`.
#' @export
fit_sequential_zentangler <- function(
  mae,
  x_var,
  y_var = NULL,
  stage_views = NULL,
  path_templates = NULL,
  view_names = NULL,
  assay_names = NULL,
  covariates = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  survival_time_var = NULL,
  survival_event_var = NULL,
  sis_n = 50L,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  cor_method = c("spearman", "pearson"),
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25,
  cor_fdr_method = c("BH", "BY"),
  residualize_links = TRUE,
  max_links_per_transition = 5000L,
  max_paths = 10000L,
  fusion_mode = c("early", "intermediate", "late"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  b_inference = c("debiased_lasso", "debiased_logistic_lasso", "debiased_cox_lasso", "refit", "bootstrap"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04,
  path_inference = c("model_based", "bootstrap_score"),
  bootstrap_repeats = 0L,
  bootstrap_ci_level = 0.95,
  bootstrap_seed = seed + 10000L,
  bootstrap_id = NULL,
  causal_inference = c("none", "bootstrap"),
  effect_x0 = NULL,
  effect_x1 = NULL,
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global"),
  seed = 1,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  start_time <- Sys.time()
  duplicate_primary <- match.arg(duplicate_primary)
  y_family <- match.arg(y_family)
  sis_rank <- match.arg(sis_rank)
  screen_method <- match.arg(screen_method)
  cor_method <- match.arg(cor_method)
  cor_fdr_method <- validate_fdr_method(cor_fdr_method)
  fusion_mode <- match.arg(fusion_mode)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)
  b_inference <- match.arg(b_inference)
  debias_max_targets <- max(1L, as.integer(debias_max_targets))
  path_inference <- validate_sequential_path_inference(path_inference)
  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  bootstrap_ci_level <- as.numeric(bootstrap_ci_level)
  causal_inference <- match.arg(causal_inference)
  fdr_method <- validate_fdr_method(fdr_method)
  fdr_scope <- match.arg(fdr_scope)
  set.seed(seed)

  if (is.null(x_var) || length(x_var) != 1L || !nzchar(as.character(x_var))) {
    stop("x_var is required.", call. = FALSE)
  }
  if (!identical(y_family, "survival") && (is.null(y_var) || length(y_var) != 1L || !nzchar(as.character(y_var)))) {
    stop("y_var is required unless y_family = 'survival'.", call. = FALSE)
  }

  inputs <- zentangler_mae_to_blocks(
    mae = mae,
    view_names = view_names,
    assay_names = assay_names,
    duplicate_primary = duplicate_primary,
    fill_nonfinite_zero = fill_nonfinite_zero
  )
  path_templates <- zentangler_map_path_templates(path_templates, inputs$view_map)
  path_templates <- zentangler_validate_path_templates(
    path_templates,
    names(inputs$blocks)
  )
  zentangler_check_route_api(path_templates)
  if (!is.null(path_templates)) {
    template <- path_templates[[1L]]
    stage_views <- as.list(template)
    names(stage_views) <- paste0("layer", seq_along(template))
  }
  stage_views <- zentangler_map_stage_views(stage_views, inputs$view_map)
  stage_views <- zentangler_validate_stage_views(stage_views, names(inputs$blocks))
  route_set_type <- zentangler_route_set_type(stage_views = stage_views)
  needed <- unique(c(x_var, y_var, survival_time_var, survival_event_var, covariates))
  al <- align_samples_multiview_blocks(inputs$blocks, inputs$pheno_df, needed_pheno_cols = needed)
  blocks <- al$blocks
  pheno_df <- al$pheno
  X <- pheno_df[[x_var]]
  if (!is.numeric(X)) stop("x_var must be numeric.", call. = FALSE)
  if (identical(y_family, "survival")) {
    if (is.null(survival_time_var) || is.null(survival_event_var)) {
      stop("survival_time_var and survival_event_var are required for survival outcomes.", call. = FALSE)
    }
    time_raw <- suppressWarnings(as.numeric(pheno_df[[survival_time_var]]))
    event_raw <- suppressWarnings(as.numeric(pheno_df[[survival_event_var]]))
    Y <- survival::Surv(time_raw, event_raw)
  } else {
    Y <- pheno_df[[y_var]]
    if (!is.numeric(Y)) stop("y_var must be numeric.", call. = FALSE)
  }
  C <- make_covariate_matrix(pheno_df, covariates)

  data_by_key <- list()
  feature_meta <- data.frame(
    key = character(), stage = integer(), stage_name = character(),
    view = character(), mediator = character(), stringsAsFactors = FALSE
  )
  for (stage_i in seq_along(stage_views)) {
    for (view in stage_views[[stage_i]]) {
      M <- blocks[[view]]
      keys <- zentangler_feature_key(view, colnames(M))
      for (j in seq_along(keys)) data_by_key[[keys[j]]] <- M[, j]
      feature_meta <- rbind(
        feature_meta,
        data.frame(
          key = keys,
          stage = stage_i,
          stage_name = names(stage_views)[stage_i],
          view = view,
          mediator = colnames(M),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  layer1_views <- stage_views[[1L]]
  first_screens <- lapply(layer1_views, function(view) {
    hima_screen_mediators(
      M = blocks[[view]],
      X = X,
      C = C,
      sis_n = sis_n,
      rank_by = sis_rank,
      screen_method = screen_method,
      fdr_method = fdr_method
    )
  })
  names(first_screens) <- layer1_views

  first_selected <- unlist(lapply(names(first_screens), function(view) {
    zentangler_feature_key(view, first_screens[[view]]$selected)
  }), use.names = FALSE)
  if (length(first_selected) == 0L) stop("No first-layer mediators were selected.", call. = FALSE)

  a_lookup <- do.call(rbind, lapply(names(first_screens), function(view) {
    tab <- first_screens[[view]]$table
    data.frame(
      key = zentangler_feature_key(view, tab$mediator),
      a = tab$a_hat,
      p_a = tab$p_a,
      q_a = tab$q_a,
      selected = tab$selected,
      stringsAsFactors = FALSE
    )
  }))

  transition_links <- list()
  current_keys <- first_selected
  paths <- lapply(first_selected, function(key) list(keys = key, link_rows = integer(0)))

  for (stage_i in seq_len(length(stage_views) - 1L)) {
    to_views <- stage_views[[stage_i + 1L]]
    to_keys <- feature_meta$key[feature_meta$stage == stage_i + 1L]
    from_df <- as.data.frame(data_by_key[current_keys], check.names = FALSE)
    to_df <- as.data.frame(data_by_key[to_keys], check.names = FALSE)
    rownames(from_df) <- rownames(pheno_df)
    rownames(to_df) <- rownames(pheno_df)

    links <- zentangler_screen_cor_links(
      from_df = from_df,
      to_df = to_df,
      X = X,
      C = C,
      cor_method = cor_method,
      cor_fdr_method = cor_fdr_method,
      min_abs_cor = min_abs_cor,
      cor_q_threshold = cor_q_threshold,
      residualize_links = residualize_links,
      max_links = max_links_per_transition
    )
    if (nrow(links) == 0L) {
      warning("No correlation-screened links retained for transition ", stage_i, " -> ", stage_i + 1L, ".")
      transition_links[[paste0("stage", stage_i, "_to_stage", stage_i + 1L)]] <- links
      paths <- list()
      break
    }

    from_meta <- feature_meta[match(links$from_key, feature_meta$key), , drop = FALSE]
    to_meta <- feature_meta[match(links$to_key, feature_meta$key), , drop = FALSE]
    links$transition <- stage_i
    links$from_stage <- stage_i
    links$to_stage <- stage_i + 1L
    links$from_view <- from_meta$view
    links$from_mediator <- from_meta$mediator
    links$to_view <- to_meta$view
    links$to_mediator <- to_meta$mediator
    links <- zentangler_estimate_link_effects(links, data_by_key = data_by_key, X = X, C = C)
    transition_links[[paste0("stage", stage_i, "_to_stage", stage_i + 1L)]] <- links

    paths <- zentangler_extend_paths(paths, links, max_paths = max_paths)
    current_keys <- unique(vapply(paths, function(path) tail(path$keys, 1L), character(1)))
    if (length(paths) == 0L || length(current_keys) == 0L) break
  }

  bootstrap <- NULL
  effect_decomposition <- data.frame()
  causal_stage <- list()
  if (length(paths) == 0L) {
    sequential_paths <- data.frame()
    terminal_effects <- data.frame()
  } else {
    terminal_keys <- unique(vapply(paths, function(path) tail(path$keys, 1L), character(1)))
    terminal_df <- as.data.frame(data_by_key[terminal_keys], check.names = FALSE)
    rownames(terminal_df) <- rownames(pheno_df)
    terminal_effects <- zentangler_fit_terminal_effects(
      Y = Y,
      X = X,
      terminal_df = terminal_df,
      C = C,
      y_family = y_family,
      fusion_mode = fusion_mode,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha,
      b_inference = b_inference,
      debias_max_targets = debias_max_targets,
      coop_rho = coop_rho,
      coop_maxit = coop_maxit,
      coop_tol = coop_tol
    )

    link_all <- do.call(rbind, transition_links)
    path_rows <- vector("list", length(paths))
    for (i in seq_along(paths)) {
      keys <- paths[[i]]$keys
      terminal_key <- tail(keys, 1L)
      if (length(keys) > 1L) {
        link_hits <- lapply(seq_len(length(keys) - 1L), function(j) {
          hit <- which(link_all$from_key == keys[j] & link_all$to_key == keys[j + 1L])
          if (length(hit) == 0L) NA_integer_ else hit[[1L]]
        })
        link_hits <- unlist(link_hits)
        link_dat <- link_all[link_hits, , drop = FALSE]
      } else {
        link_dat <- data.frame(d_hat = numeric(0), p_d = numeric(0), abs_cor = numeric(0), q_cor = numeric(0))
      }
      a_row <- a_lookup[match(keys[[1L]], a_lookup$key), , drop = FALSE]
      b_row <- terminal_effects[match(terminal_key, terminal_effects$terminal_key), , drop = FALSE]
      path_meta <- feature_meta[match(keys, feature_meta$key), , drop = FALSE]
      d_vals <- link_dat$d_hat
      score <- as.numeric(a_row$a) * prod(d_vals) * as.numeric(b_row$b)
      path_rows[[i]] <- data.frame(
        path_id = i,
        n_layers = length(keys),
        route_set_type = route_set_type,
        route_type = "sequential",
        mediator_path = paste(paste(path_meta$view, path_meta$mediator, sep = ":"), collapse = " -> "),
        key_path = paste(keys, collapse = " -> "),
        first_view = path_meta$view[[1L]],
        first_mediator = path_meta$mediator[[1L]],
        terminal_view = tail(path_meta$view, 1L),
        terminal_mediator = tail(path_meta$mediator, 1L),
        terminal_key = terminal_key,
        a = as.numeric(a_row$a),
        d_product = prod(d_vals),
        b = as.numeric(b_row$b),
        sequential_score = score,
        abs_sequential_score = abs(score),
        p_a = as.numeric(a_row$p_a),
        max_p_d = if (any(is.finite(link_dat$p_d))) max(link_dat$p_d, na.rm = TRUE) else NA_real_,
        p_b = as.numeric(b_row$p_b),
        fusion_mode = as.character(b_row$fusion_mode),
        p_primary = zentangler_path_primary_p(
          p_a = as.numeric(a_row$p_a),
          p_d = link_dat$p_d,
          p_b = as.numeric(b_row$p_b),
          b = as.numeric(b_row$b)
        ),
        min_abs_cor = if (any(is.finite(link_dat$abs_cor))) min(link_dat$abs_cor, na.rm = TRUE) else NA_real_,
        max_q_cor = if (any(is.finite(link_dat$q_cor))) max(link_dat$q_cor, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    sequential_paths <- do.call(rbind, path_rows)
    sequential_paths$q_primary_global <- stats::p.adjust(sequential_paths$p_primary, method = fdr_method)
    sequential_paths$q_primary <- sequential_paths$q_primary_global
    effect_decomposition <- compute_effect_decomposition_sequential(
      sequential_paths = sequential_paths,
      Y = Y,
      X = X,
      data_by_key = data_by_key,
      C = C,
      y_family = y_family,
      x0 = effect_x0,
      x1 = effect_x1
    )
    bootstrap <- NULL
    if (bootstrap_repeats > 0L) {
      bootstrap <- zentangler_bootstrap_sequential_paths(
        paths = sequential_paths,
        data_by_key = data_by_key,
        pheno_df = pheno_df,
        X = X,
        Y = Y,
        C = C,
        y_family = y_family,
        fusion_mode = fusion_mode,
        lambda_choice = lambda_choice,
        glmnet_alpha = glmnet_alpha,
        b_inference = b_inference,
        causal_inference = causal_inference,
        debias_max_targets = debias_max_targets,
        coop_rho = coop_rho,
        coop_maxit = coop_maxit,
        coop_tol = coop_tol,
        repeats = bootstrap_repeats,
        ci_level = bootstrap_ci_level,
        seed = bootstrap_seed,
        bootstrap_id = bootstrap_id,
        effect_x0 = effect_x0,
        effect_x1 = effect_x1
      )
      sequential_paths <- zentangler_add_sequential_bootstrap_columns(
        sequential_paths,
        bootstrap = bootstrap,
        fdr_method = fdr_method,
        use_primary = identical(path_inference, "bootstrap_score")
      )
      if (identical(causal_inference, "bootstrap")) {
        effect_decomposition <- add_effect_bootstrap_columns(
          effect_decomposition = effect_decomposition,
          effect_boot = bootstrap$effect_decomposition,
          ci_level = bootstrap_ci_level
        )
      }
    } else if (identical(path_inference, "bootstrap_score")) {
      warning("path_inference='bootstrap_score' requires bootstrap_repeats > 0; using model-based path p-values.")
    }
    sequential_paths <- sequential_paths[order(sequential_paths$abs_sequential_score, decreasing = TRUE), , drop = FALSE]
    rownames(sequential_paths) <- NULL
    causal_stage <- build_causal_stage_summary(
      effect_decomposition = effect_decomposition,
      combined_mediators = sequential_paths,
      causal_inference = causal_inference,
      bootstrap = bootstrap
    )
  }

  end_time <- Sys.time()
  out <- list(
    settings = list(
      input_container = "MultiAssayExperiment",
      model = route_set_type,
      route_set_type = route_set_type,
      x_var = x_var,
      y_var = y_var,
      y_family = y_family,
      survival_time_var = survival_time_var,
      survival_event_var = survival_event_var,
      stage_views = stage_views,
      covariates = covariates,
      sis_n = sis_n,
      sis_rank = sis_rank,
      screen_method = screen_method,
      cor_method = cor_method,
      min_abs_cor = min_abs_cor,
      cor_q_threshold = cor_q_threshold,
      cor_fdr_method = cor_fdr_method,
      residualize_links = residualize_links,
      fusion_mode = fusion_mode,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha,
      b_inference = b_inference,
      debias_max_targets = debias_max_targets,
      coop_rho = coop_rho,
      coop_maxit = coop_maxit,
      coop_tol = coop_tol,
      path_inference = path_inference,
      bootstrap_repeats = bootstrap_repeats,
      bootstrap_ci_level = bootstrap_ci_level,
      bootstrap_seed = bootstrap_seed,
      bootstrap_id = bootstrap_id,
      causal_inference = causal_inference,
      effect_x0 = effect_decomposition$x0 %||% effect_x0,
      effect_x1 = effect_decomposition$x1 %||% effect_x1,
      effect_method = effect_decomposition$effect_method %||% NA_character_,
      effect_scale = effect_decomposition$effect_scale %||% NA_character_,
      effects_summary = if (nrow(effect_decomposition) > 0L) as.list(effect_decomposition[1, , drop = FALSE]) else list(),
      causal_stage = causal_stage$stage %||% NA_character_,
      causal_stage_method = causal_stage$stage_method %||% NA_character_,
      fdr_method = fdr_method,
      fdr_scope = fdr_scope
    ),
    diagnostics = list(
      started = start_time,
      finished = end_time,
      runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
      n_samples = nrow(pheno_df),
      n_layers = length(stage_views),
      n_paths = nrow(sequential_paths),
      n_active_paths_q025 = sum(is.finite(sequential_paths$q_primary) & sequential_paths$q_primary <= 0.25, na.rm = TRUE),
      n_first_layer_selected = length(first_selected),
      n_transition_links = sum(vapply(transition_links, nrow, integer(1))),
      n_terminal_effects = nrow(terminal_effects),
      n_nonzero_terminal_b = if (nrow(terminal_effects) > 0L) sum(is.finite(terminal_effects$b) & terminal_effects$b != 0, na.rm = TRUE) else 0L,
      bootstrap_repeats_requested = bootstrap_repeats,
      bootstrap_repeats_successful = if (!is.null(bootstrap)) bootstrap$repeats_successful else 0L,
      bootstrap_failures = if (!is.null(bootstrap)) bootstrap$failures else character(0)
    ),
    sample_ids = rownames(pheno_df),
    feature_metadata = feature_meta,
    first_layer_screens = first_screens,
    transition_links = transition_links,
    terminal_effects = terminal_effects,
    sequential_paths = sequential_paths,
    paths = sequential_paths,
    causal_stage = causal_stage,
    effects = effect_decomposition,
    effect_decomposition = effect_decomposition,
    causal_effects = effect_decomposition,
    bootstrap = bootstrap
  )
  class(out) <- c("zentangler_sequential_fit", "list")
  out
}

#' Extract sequential Zentangler paths
#'
#' @param fit Output from `fit_sequential_zentangler()`.
#' @return Data frame of sequential paths.
#' @export
zentangler_sequential_paths <- function(fit) {
  out <- fit$sequential_paths %||% fit$paths
  if (is.null(out)) return(data.frame())
  rownames(out) <- NULL
  out
}

#' Extract sequential Zentangler transition edges
#'
#' @param fit Output from `fit_sequential_zentangler()`.
#' @return Data frame of retained adjacent mediator-to-mediator links.
#' @export
zentangler_sequential_edges <- function(fit) {
  if (!is.null(fit$transition_links) && is.data.frame(fit$transition_links)) {
    out <- fit$transition_links
    rownames(out) <- NULL
    return(out)
  }
  if (!is.null(fit$transition_links) && is.list(fit$transition_links)) {
    out <- do.call(rbind, fit$transition_links)
    if (is.null(out)) return(data.frame())
    rownames(out) <- NULL
    return(out)
  }
  data.frame()
}

#' Extract sequential Zentangler terminal effects
#'
#' @param fit Output from `fit_sequential_zentangler()`.
#' @return Data frame of terminal mediator-to-outcome coefficients.
#' @export
zentangler_sequential_terminals <- function(fit) {
  out <- fit$terminal_effects
  if (is.null(out)) return(data.frame())
  rownames(out) <- NULL
  out
}

#' Extract sequential Zentangler diagnostics
#'
#' @param fit Output from `fit_sequential_zentangler()`.
#' @return Data frame with run-level diagnostic counts.
#' @export
zentangler_sequential_diagnostics <- function(fit) {
  d <- fit$diagnostics %||% list()
  if (length(d) == 0L) return(data.frame())
  scalar <- vapply(d, function(x) length(x) == 1L, logical(1))
  out <- as.data.frame(d[scalar], stringsAsFactors = FALSE)
  if ("bootstrap_failures" %in% names(d)) {
    out$n_bootstrap_failures <- length(d$bootstrap_failures)
  }
  rownames(out) <- NULL
  out
}

zentangler_sequential_active_paths <- function(fit, q_threshold = 0.25, q_col = "q_primary") {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  tab <- zentangler_sequential_paths(fit)
  if (nrow(tab) == 0L) return(tab)
  if (!(q_col %in% colnames(tab))) stop("q_col not found in sequential path table: ", q_col)
  keep <- is.finite(tab[[q_col]]) & tab[[q_col]] <= q_threshold & is.finite(tab$b) & tab$b != 0
  out <- tab[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

zentangler_sequential_top_paths <- function(
  fit,
  n = 20L,
  order_by = c("abs_sequential_score", "q_primary", "p_primary", "sequential_score"),
  decreasing = NULL
) {
  order_by <- match.arg(order_by)
  tab <- zentangler_sequential_paths(fit)
  if (nrow(tab) == 0L) return(tab)
  if (!(order_by %in% colnames(tab))) stop("order_by not found in sequential path table: ", order_by)
  if (is.null(decreasing)) decreasing <- order_by %in% c("abs_sequential_score", "sequential_score")
  ord <- order(tab[[order_by]], decreasing = decreasing, na.last = TRUE)
  out <- tab[ord, , drop = FALSE]
  out <- out[seq_len(min(as.integer(n), nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

zentangler_sequential_threshold_summary <- function(
  fit,
  q_threshold = seq(0.05, 0.25, by = 0.05),
  q_cols = NULL
) {
  tab <- zentangler_sequential_paths(fit)
  if (nrow(tab) == 0L) return(data.frame())
  q_threshold <- zentangler_clean_q_thresholds(q_threshold)
  if (is.null(q_cols)) q_cols <- intersect(c("q_primary", "q_primary_global"), colnames(tab))
  q_cols <- intersect(q_cols, colnames(tab))
  if (length(q_cols) == 0L) stop("No valid q-value columns found in sequential path table.")

  rows <- list()
  ii <- 0L
  for (q_col in q_cols) {
    for (threshold in q_threshold) {
      active <- is.finite(tab[[q_col]]) & tab[[q_col]] <= threshold & is.finite(tab$b) & tab$b != 0
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        q_col = q_col,
        q_threshold = threshold,
        n_active_paths = sum(active, na.rm = TRUE),
        n_positive_paths = sum(active & is.finite(tab$sequential_score) & tab$sequential_score > 0, na.rm = TRUE),
        n_negative_paths = sum(active & is.finite(tab$sequential_score) & tab$sequential_score < 0, na.rm = TRUE),
        min_q = if (any(is.finite(tab[[q_col]]))) min(tab[[q_col]], na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

zentangler_sequential_model_summary <- function(fit, q_threshold = 0.25) {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  tab <- zentangler_sequential_paths(fit)
  active <- if (nrow(tab) > 0L) {
    is.finite(tab$q_primary) & tab$q_primary <= q_threshold & is.finite(tab$b) & tab$b != 0
  } else {
    logical(0)
  }
  settings <- fit$settings %||% list()
  diagnostics <- fit$diagnostics %||% list()
  effects <- settings$effects_summary %||% list()
  data.frame(
    input_container = settings$input_container %||% NA_character_,
    model = settings$model %||% settings$route_set_type %||% NA_character_,
    route_set_type = settings$route_set_type %||% NA_character_,
    n_samples = diagnostics$n_samples %||% NA_integer_,
    n_paths = nrow(tab),
    n_active_paths = sum(active, na.rm = TRUE),
    n_terminal_effects = diagnostics$n_terminal_effects %||% nrow(zentangler_sequential_terminals(fit)),
    n_nonzero_terminal_b = diagnostics$n_nonzero_terminal_b %||% NA_integer_,
    n_transition_links = diagnostics$n_transition_links %||% nrow(zentangler_sequential_edges(fit)),
    fusion_mode = settings$fusion_mode %||% NA_character_,
    lambda_choice = settings$lambda_choice %||% NA_character_,
    glmnet_alpha = settings$glmnet_alpha %||% NA_real_,
    b_inference = settings$b_inference %||% NA_character_,
    path_inference = settings$path_inference %||% "model_based",
    causal_inference = settings$causal_inference %||% NA_character_,
    causal_stage = settings$causal_stage %||% NA_character_,
    causal_stage_method = settings$causal_stage_method %||% NA_character_,
    effect_method = settings$effect_method %||% NA_character_,
    effect_scale = settings$effect_scale %||% NA_character_,
    nde = effects$nde %||% NA_real_,
    nie_total = effects$nie_total %||% NA_real_,
    te = effects$te %||% NA_real_,
    prop_mediated = effects$prop_mediated %||% NA_real_,
    fdr_method = settings$fdr_method %||% NA_character_,
    runtime_seconds = diagnostics$runtime_seconds %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Summarize a sequential Zentangler fit
#'
#' @param fit Output from `fit_sequential_zentangler()`.
#' @param q_threshold One or more q-value thresholds.
#' @return A list of model, threshold, diagnostics, path, edge, and terminal
#'   summaries.
#' @export
summarize_sequential_zentangler <- function(fit, q_threshold = 0.25) {
  list(
    model_summary = zentangler_sequential_model_summary(fit, q_threshold = max(zentangler_clean_q_thresholds(q_threshold))),
    threshold_summary = zentangler_sequential_threshold_summary(fit, q_threshold = q_threshold),
    diagnostics = zentangler_sequential_diagnostics(fit),
    top_paths = zentangler_sequential_top_paths(fit, n = 20L),
    active_paths = zentangler_sequential_active_paths(fit, q_threshold = max(zentangler_clean_q_thresholds(q_threshold))),
    edges = zentangler_sequential_edges(fit),
    terminals = zentangler_sequential_terminals(fit),
    causal_effects = zentangler_effects(fit),
    effects = zentangler_effects(fit)
  )
}
