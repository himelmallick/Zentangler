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
  combined$p_score_boot <- bootstrap_score_pvalues(observed = combined$score, score_mat = score_mat)
  combined$p_score_boot_sign <- bootstrap_score_sign_pvalues(score_mat = score_mat)
  combined$p_primary_bootstrap <- combined$p_score_boot
  combined
}

bootstrap_score_pvalues <- function(observed, score_mat) {
  out <- rep(NA_real_, ncol(score_mat))
  names(out) <- colnames(score_mat)
  for (j in seq_len(ncol(score_mat))) {
    obs <- as.numeric(observed[j])
    vals <- as.numeric(score_mat[, j])
    vals <- vals[is.finite(vals)]
    if (!is.finite(obs) || length(vals) < 2L) next

    # Null-centered bootstrap test for H0: score = 0. The centered bootstrap
    # distribution approximates estimator fluctuation under the null.
    centered <- vals - obs
    out[j] <- (1 + sum(abs(centered) >= abs(obs), na.rm = TRUE)) / (length(centered) + 1)
  }
  pmin(1, pmax(0, out))
}

bootstrap_score_sign_pvalues <- function(score_mat) {
  out <- rep(NA_real_, ncol(score_mat))
  names(out) <- colnames(score_mat)
  for (j in seq_len(ncol(score_mat))) {
    vals <- as.numeric(score_mat[, j])
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L) next
    p_low <- (1 + sum(vals <= 0, na.rm = TRUE)) / (length(vals) + 1)
    p_high <- (1 + sum(vals >= 0, na.rm = TRUE)) / (length(vals) + 1)
    out[j] <- min(1, 2 * min(p_low, p_high))
  }
  pmin(1, pmax(0, out))
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
  method_preset,
  covariates,
  residualize,
  sis_n,
  sis_rank,
  screen_method,
  a_stage_model,
  maaslin2_random_effect,
  maaslin2_normalization,
  maaslin2_transform,
  maaslin2_analysis_method,
  maaslin2_standardize,
  maaslin2_output_dir,
  fusion_mode,
  y_family,
  lambda_choice,
  glmnet_alpha,
  fdr_method,
  fdr_scope,
  primary_inference,
  b_inference,
  debias_max_targets,
  coop_rho,
  coop_maxit,
  coop_tol,
  key_ref,
  bootstrap_repeats,
  bootstrap_ci_level,
  bootstrap_seed,
  bootstrap_id,
  effect_x0,
  effect_x1
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

    maaslin2_output_dir_b <- if (!is.null(maaslin2_output_dir) && nzchar(maaslin2_output_dir)) {
      file.path(maaslin2_output_dir, paste0("bootstrap_", b))
    } else {
      NULL
    }

    fit_b <- try(
      fit_multiview_parallel_zentangler_blocks(
        blocks = bt$blocks,
        pheno_df = bt$pheno_df,
        x_var = x_var,
        y_var = y_var,
        method_preset = method_preset,
        covariates = covariates,
        residualize = residualize,
        sis_n = sis_n,
        sis_rank = sis_rank,
        screen_method = screen_method,
        a_stage_model = a_stage_model,
        maaslin2_random_effect = maaslin2_random_effect,
        maaslin2_normalization = maaslin2_normalization,
        maaslin2_transform = maaslin2_transform,
        maaslin2_analysis_method = maaslin2_analysis_method,
        maaslin2_standardize = maaslin2_standardize,
        maaslin2_output_dir = maaslin2_output_dir_b,
        fusion_mode = fusion_mode,
        y_family = y_family,
        lambda_choice = lambda_choice,
        glmnet_alpha = glmnet_alpha,
        fdr_method = fdr_method,
        fdr_scope = fdr_scope,
        primary_inference = "model_based",
        b_inference = b_inference,
        debias_max_targets = debias_max_targets,
        coop_rho = coop_rho,
        coop_maxit = coop_maxit,
        coop_tol = coop_tol,
        seed = bootstrap_seed + b,
        return_fits = FALSE,
        bootstrap_repeats = 0L,
        effect_x0 = effect_x0,
        effect_x1 = effect_x1
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

