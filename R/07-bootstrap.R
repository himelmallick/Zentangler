# -----------------------------------------------------------------------------
# Bootstrap
#
# Two different bootstraps are used, each matched to what it is inferring:
#
# 1. Pipeline bootstrap (b_inference = "bootstrap", primary_inference =
#    "bootstrap_score"). Each replicate re-runs the A-stage, screening and the
#    fused outcome model, so selection uncertainty is propagated. Mediators that
#    are not selected in a replicate contribute b* = 0. p-values invert the
#    percentile interval: p = 2 * min(P*(b* <= 0), P*(b* >= 0)), with the usual
#    +1 correction; a mediator that is never selected gets p = 1.
#
# 2. Fixed-set effect bootstrap (causal_inference = "bootstrap"). The effect
#    set and the exposure contrast are fixed at their full-data values, and the
#    unpenalized effect models are refit in each replicate. Intervals are
#    therefore conditional on the selected mediator set.
#
# For clustered data (bootstrap_id), whole clusters are resampled, each
# resampled copy gets a fresh cluster label, and cross-validation folds inside
# a replicate are assigned by original cluster so duplicated rows never sit on
# both sides of a fold.
# -----------------------------------------------------------------------------

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
    idx_list <- lapply(sampled_units, function(u) which(unit == u))
    idx <- unlist(idx_list, use.names = FALSE)
    copy_label <- rep(paste0(sampled_units, "__copy", seq_along(sampled_units)), times = lengths(idx_list))
    fold_group <- rep(sampled_units, times = lengths(idx_list))
  } else {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    copy_label <- NULL
    fold_group <- rownames(pheno_df)[idx]
  }

  new_ids <- paste0(rownames(pheno_df)[idx], "__boot", seq_along(idx))
  pheno_b <- pheno_df[idx, , drop = FALSE]
  rownames(pheno_b) <- new_ids
  if (!is.null(copy_label)) pheno_b[[bootstrap_id]] <- copy_label

  blocks_b <- lapply(blocks, function(M) {
    out <- M[idx, , drop = FALSE]
    rownames(out) <- new_ids
    out
  })
  names(blocks_b) <- names(blocks)

  list(blocks = blocks_b, pheno_df = pheno_b, idx = idx, fold_group = fold_group)
}

# Run expr with cross-validation folds grouped by `groups` (see
# zentangler_make_foldid()).
zentangler_with_fold_groups <- function(groups, expr) {
  old <- options(zentangler.fold_groups = groups)
  on.exit(options(old), add = TRUE)
  force(expr)
}

bootstrap_sign_pvalues <- function(mat) {
  # Percentile-interval inversion p-value for H0: parameter = 0. Exact zeros
  # (not selected in that replicate) count toward both tails.
  out <- rep(NA_real_, ncol(mat))
  names(out) <- colnames(mat)
  for (j in seq_len(ncol(mat))) {
    vals <- as.numeric(mat[, j])
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L) next
    p_low <- (1 + sum(vals <= 0)) / (length(vals) + 1)
    p_high <- (1 + sum(vals >= 0)) / (length(vals) + 1)
    out[j] <- min(1, 2 * min(p_low, p_high))
  }
  out
}

zentangler_boot_summary_columns <- function(tab, mat, prefix, ci_level) {
  alpha <- 1 - ci_level
  q <- function(v, pr) if (sum(is.finite(v)) > 0) stats::quantile(v, probs = pr, na.rm = TRUE, names = FALSE) else NA_real_
  tab[[paste0(prefix, "_boot_mean")]] <- apply(mat, 2, function(v) if (any(is.finite(v))) mean(v, na.rm = TRUE) else NA_real_)
  tab[[paste0(prefix, "_boot_sd")]] <- apply(mat, 2, function(v) if (sum(is.finite(v)) > 1) stats::sd(v, na.rm = TRUE) else NA_real_)
  tab[[paste0(prefix, "_boot_low")]] <- apply(mat, 2, q, pr = alpha / 2)
  tab[[paste0(prefix, "_boot_high")]] <- apply(mat, 2, q, pr = 1 - alpha / 2)
  tab
}

add_bootstrap_columns <- function(combined, score_mat, active_mat, ci_level = 0.95) {
  if (is.null(score_mat) || nrow(score_mat) == 0 || nrow(combined) == 0) return(combined)
  combined <- zentangler_boot_summary_columns(combined, score_mat, "score", ci_level)
  combined$score_boot_selection_freq <- apply(active_mat, 2, mean, na.rm = TRUE)
  combined$p_score_boot <- bootstrap_sign_pvalues(score_mat)
  combined
}

add_b_bootstrap_columns <- function(combined, b_mat, ci_level = 0.95) {
  if (is.null(b_mat) || nrow(b_mat) == 0 || nrow(combined) == 0) return(combined)
  combined <- zentangler_boot_summary_columns(combined, b_mat, "b", ci_level)
  combined$p_b_bootstrap <- bootstrap_sign_pvalues(b_mat)
  combined
}

add_effect_bootstrap_columns <- function(effect_decomposition, effect_boot, ci_level = 0.95) {
  if (is.null(effect_boot) || nrow(effect_boot) == 0 || nrow(effect_decomposition) == 0) {
    return(effect_decomposition)
  }
  alpha <- 1 - ci_level
  numeric_cols <- intersect(
    colnames(effect_boot)[vapply(effect_boot, is.numeric, logical(1))],
    colnames(effect_decomposition)
  )
  numeric_cols <- setdiff(numeric_cols, c("x0", "x1", "delta_x"))
  for (nm in numeric_cols) {
    vals <- effect_boot[[nm]]
    if (!any(is.finite(vals))) next
    effect_decomposition[[paste0(nm, "_boot_mean")]] <- mean(vals, na.rm = TRUE)
    effect_decomposition[[paste0(nm, "_boot_sd")]] <- stats::sd(vals, na.rm = TRUE)
    effect_decomposition[[paste0(nm, "_boot_low")]] <- stats::quantile(vals, probs = alpha / 2, na.rm = TRUE, names = FALSE)
    effect_decomposition[[paste0(nm, "_boot_high")]] <- stats::quantile(vals, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  }
  for (nm in intersect(c("nde", "nie_total", "te"), colnames(effect_boot))) {
    vals <- effect_boot[[nm]]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L) next
    p_low <- (1 + sum(vals <= 0)) / (length(vals) + 1)
    p_high <- (1 + sum(vals >= 0)) / (length(vals) + 1)
    effect_decomposition[[paste0("p_", nm, "_bootstrap")]] <- min(1, 2 * min(p_low, p_high))
  }
  effect_decomposition
}

zentangler_bind_rows_fill <- function(rows) {
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (length(rows) == 0L) return(data.frame())
  cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    miss <- setdiff(cols, colnames(x))
    for (m in miss) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

zentangler_check_bootstrap_success <- function(n_ok, n_req, label = "Bootstrap") {
  if (n_req > 0L && n_ok < 0.8 * n_req) {
    warning(
      label, ": only ", n_ok, " of ", n_req, " replicates succeeded. ",
      "Failures are rarely random, so intervals and p-values may be biased; inspect diagnostics$bootstrap_failures.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

bootstrap_multiview_fit <- function(blocks, pheno_df, opts, key_ref, bootstrap_repeats, bootstrap_seed, bootstrap_id = NULL, verbose = FALSE) {
  # Pipeline bootstrap: every replicate re-runs screening and the outcome model
  # through the same engine as the main fit (point estimates only).
  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  if (bootstrap_repeats < 1L) return(NULL)
  bootstrap_seed <- as.integer(bootstrap_seed %||% 10001L)
  score_mat <- matrix(NA_real_, nrow = bootstrap_repeats, ncol = length(key_ref), dimnames = list(NULL, key_ref))
  b_mat <- score_mat
  active_mat <- score_mat
  failures <- character(0)
  n_ok <- 0L
  rep_opts <- opts
  rep_opts[c("b_inference", "causal_inference", "primary_inference", "fdr_joint")] <- list("none", "none", "model_based", "maxp")
  rep_opts$bootstrap_repeats <- 0L
  rep_opts$return_fits <- FALSE
  rep_opts$verbose <- FALSE

  for (b in seq_len(bootstrap_repeats)) {
    bt <- try(make_bootstrap_sample(blocks, pheno_df, bootstrap_id = bootstrap_id, seed = bootstrap_seed + b), silent = TRUE)
    if (inherits(bt, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(bt, "condition"))))
      next
    }
    o_b <- rep_opts
    o_b$seed <- bootstrap_seed + b
    # Write each replicate's MaAsLin2 / Tweedieverse output to its own folder.
    o_b$a_stage <- zentangler_a_stage_bootstrap_dirs(o_b$a_stage, b)
    fit_b <- try(
      suppressMessages(suppressWarnings(zentangler_with_fold_groups(
        bt$fold_group,
        zentangler_parallel_engine(bt$blocks, bt$pheno_df, o_b, skip_effects = TRUE)
      ))),
      silent = TRUE
    )
    if (inherits(fit_b, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(fit_b, "condition"))))
      next
    }
    tab <- fit_b$mediators
    hit <- match(key_ref, paste(tab$omics, tab$mediator, sep = "::"))
    ok <- !is.na(hit)
    # Features dropped or not screened in a replicate count as not selected.
    b_row <- setNames(rep(0, length(key_ref)), key_ref)
    s_row <- b_row
    b_row[ok] <- tab$b_fusion[hit[ok]]
    s_row[ok] <- tab$a[hit[ok]] * tab$b_fusion[hit[ok]]
    s_row[!is.finite(s_row)] <- 0
    b_mat[b, ] <- b_row
    score_mat[b, ] <- s_row
    active_mat[b, ] <- as.numeric(b_row != 0)
    n_ok <- n_ok + 1L
    if (b %% zentangler_progress_every(bootstrap_repeats) == 0L || b == bootstrap_repeats) {
      zentangler_log(verbose, sprintf("bootstrap %d / %d (%d succeeded)", b, bootstrap_repeats, n_ok), indent = 1L)
    }
  }
  zentangler_check_bootstrap_success(n_ok, bootstrap_repeats, "Pipeline bootstrap")
  list(score_matrix = score_mat, b_matrix = b_mat, active_matrix = active_mat,
       repeats_requested = bootstrap_repeats, repeats_successful = n_ok, failures = failures)
}

bootstrap_effects_fixed_set <- function(
  blocks,
  pheno_df,
  X,
  Y,
  C,
  effect_set,
  y_family,
  x0,
  x1,
  bootstrap_repeats,
  bootstrap_seed,
  bootstrap_id = NULL,
  verbose = FALSE
) {
  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  if (bootstrap_repeats < 1L) return(NULL)
  if (is.null(bootstrap_seed) || !is.finite(as.numeric(bootstrap_seed))) bootstrap_seed <- 10001L
  bootstrap_seed <- as.integer(bootstrap_seed) + 500000L
  med_keys <- paste(effect_set$omics, effect_set$mediator, sep = "::")
  ind_mat <- matrix(NA_real_, nrow = bootstrap_repeats, ncol = length(med_keys), dimnames = list(NULL, med_keys))
  rows <- list()
  failures <- character(0)

  for (b in seq_len(bootstrap_repeats)) {
    bt <- try(make_bootstrap_sample(blocks, pheno_df, bootstrap_id = bootstrap_id, seed = bootstrap_seed + b), silent = TRUE)
    if (inherits(bt, "try-error")) {
      failures <- c(failures, paste0("effect_bootstrap_", b, ": ", conditionMessage(attr(bt, "condition"))))
      next
    }
    idx <- bt$idx
    eff_b <- try(
      suppressWarnings(compute_effect_decomposition_multiview(
        effect_set = effect_set,
        Y = zentangler_subset_y(Y, idx),
        X = X[idx],
        blocks = bt$blocks,
        C = if (is.null(C)) NULL else C[idx, , drop = FALSE],
        y_family = y_family,
        x0 = x0,
        x1 = x1
      )),
      silent = TRUE
    )
    if (inherits(eff_b, "try-error") || !any(is.finite(eff_b$effects$nie_total))) {
      msg <- if (inherits(eff_b, "try-error")) conditionMessage(attr(eff_b, "condition")) else "effect model failed"
      failures <- c(failures, paste0("effect_bootstrap_", b, ": ", msg))
      next
    }
    rows[[length(rows) + 1L]] <- eff_b$effects
    if (b %% zentangler_progress_every(bootstrap_repeats) == 0L || b == bootstrap_repeats) {
      zentangler_log(verbose, sprintf("effect bootstrap %d / %d (%d succeeded)", b, bootstrap_repeats, length(rows)), indent = 1L)
    }
    me <- eff_b$mediator_effects
    if (nrow(me) > 0L) {
      k <- paste(me$omics, me$mediator, sep = "::")
      ind_mat[b, match(k, med_keys)] <- me$indirect_effect
    }
  }
  zentangler_check_bootstrap_success(length(rows), bootstrap_repeats, "Causal-effect bootstrap")
  list(
    effect_decomposition = zentangler_bind_rows_fill(rows),
    indirect_matrix = ind_mat,
    repeats_requested = bootstrap_repeats,
    repeats_successful = length(rows),
    failures = failures
  )
}

add_mediator_effect_bootstrap_columns <- function(combined, indirect_mat, ci_level = 0.95) {
  if (is.null(indirect_mat) || ncol(indirect_mat) == 0L || is.null(combined) || nrow(combined) == 0) return(combined)
  key_tab <- paste(combined$omics, combined$mediator, sep = "::")
  full <- matrix(NA_real_, nrow = nrow(indirect_mat), ncol = length(key_tab), dimnames = list(NULL, key_tab))
  hit <- intersect(colnames(indirect_mat), key_tab)
  full[, hit] <- indirect_mat[, hit]
  combined <- zentangler_boot_summary_columns(combined, full, "indirect_effect", ci_level)
  combined$p_indirect_bootstrap <- bootstrap_sign_pvalues(full)
  combined
}
