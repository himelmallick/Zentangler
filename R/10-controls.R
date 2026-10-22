# -----------------------------------------------------------------------------
# Negative controls (outcome permutation) and direction checks
#
# Both work for parallel fits (route = NULL) and sequential fits (route given)
# and use the same per-unit view of a fit: mediators for parallel fits, paths
# for sequential fits.
# -----------------------------------------------------------------------------

zentangler_fit_units <- function(fit, q_threshold) {
  # One row per tested unit (mediator or path) with a common schema.
  if (inherits(fit, "zentangler_sequential_fit")) {
    tab <- fit$paths
    if (is.null(tab) || nrow(tab) == 0L) return(data.frame())
    id <- tab$key_path
    score <- tab$abs_sequential_score
  } else {
    tab <- fit$mediators
    if (is.null(tab) || nrow(tab) == 0L) return(data.frame())
    id <- paste(tab$omics, tab$mediator, sep = "::")
    score <- tab$abs_score
  }
  data.frame(
    unit_id = id,
    abs_score = score,
    q_primary = tab$q_primary,
    is_active = zentangler_is_active(tab, q_threshold),
    stringsAsFactors = FALSE
  )
}

zentangler_unit_summary <- function(units) {
  data.frame(
    n_units = nrow(units),
    n_active = if (nrow(units)) sum(units$is_active) else 0L,
    best_q = if (nrow(units) && any(is.finite(units$q_primary))) min(units$q_primary, na.rm = TRUE) else NA_real_,
    max_abs_score = if (nrow(units) && any(is.finite(units$abs_score))) max(units$abs_score, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

zentangler_perm_p <- function(null_vals, obs, upper = TRUE) {
  # Permutation p-value with the +1 correction (Phipson & Smyth 2010).
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0L) return(NA_real_)
  hits <- if (upper) sum(null_vals >= obs) else sum(null_vals <= obs)
  (1 + hits) / (1 + length(null_vals))
}

zentangler_permute_outcome <- function(mae, vars, eligible_ids, seed) {
  # Permute the outcome column(s) jointly (survival time and event together),
  # only among the samples the analysis uses, so null fits see the same
  # samples as the observed fit.
  set.seed(seed)
  cd <- MultiAssayExperiment::colData(mae)
  rows <- which(rownames(cd) %in% eligible_ids)
  if (length(rows) < 2L) stop("Fewer than two eligible samples to permute.", call. = FALSE)
  perm <- rows[sample.int(length(rows))]
  for (v in vars) {
    vals <- as.vector(cd[[v]])
    vals[rows] <- vals[perm]
    cd[[v]] <- vals
  }
  MultiAssayExperiment::colData(mae) <- cd
  mae
}

zentangler_eligible_sample_ids <- function(mae, needed_vars, views = NULL) {
  cd <- as.data.frame(zentangler_mae_coldata(mae))
  needed_vars <- intersect(unique(needed_vars), colnames(cd))
  ids <- rownames(cd)[stats::complete.cases(cd[, needed_vars, drop = FALSE])]
  sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  for (v in views %||% unique(as.character(sm$assay))) {
    ids <- intersect(ids, as.character(sm$primary[as.character(sm$assay) == v]))
  }
  ids
}

zentangler_run_fit <- function(mae, x_var, y_var, model, route, seed, dots) {
  do.call(zentangler, c(list(mae = mae, model = model, x_var = x_var, y_var = y_var, route = route, seed = seed), dots))
}

#' Outcome-permutation negative controls
#'
#' Fits the model once on the observed data and `n_perm` times with the
#' outcome permuted among the analysed samples (survival time and event are
#' permuted jointly). This calibrates discoveries against the global null of
#' no outcome association.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param x_var Exposure column.
#' @param y_var Outcome column (`NULL` for survival; pass `y_family`,
#'   `survival_time_var` and `survival_event_var` through `...`).
#' @param model `"parallel"` or `"sequential"`, as in [zentangler()].
#' @param route Character route for the sequential model.
#' @param n_perm Number of permutations.
#' @param q_threshold q-value cutoff defining active mediators/paths.
#' @param seed Base seed.
#' @param return_null_fits Keep the permuted fits.
#' @param verbose Log the observed fit and permutation progress; permuted
#'   fits themselves run quietly.
#' @param ... Further arguments for [zentangler()].
#' @return A list with `observed_fit`, `observed_summary`, `null_summaries`,
#'   `empirical_summary` (permutation p-values for the number of discoveries,
#'   the best q-value and the largest score) and `empirical_units`
#'   (per-mediator or per-path null selection frequency and score p-value).
#' @export
run_negative_controls <- function(mae, x_var, y_var = NULL, model = c("parallel", "sequential"), route = NULL, n_perm = 100L, q_threshold = 0.25, seed = 1, return_null_fits = FALSE, verbose = TRUE, ...) {
  model <- zentangler_resolve_model(model, route, missing(model))
  n_perm <- as.integer(n_perm)
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) stop("n_perm must be a single integer >= 1.", call. = FALSE)
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  restore_rng <- zentangler_save_rng()
  on.exit(restore_rng(), add = TRUE)
  dots <- list(...)
  survival <- identical(dots$y_family, "survival")
  perm_vars <- if (survival) c(dots$survival_time_var, dots$survival_event_var) else y_var
  if (length(perm_vars) == 0L) stop("Could not determine the outcome column(s) to permute.", call. = FALSE)
  views <- if (identical(model, "parallel")) dots$view_names else zentangler_resolve_route(route, mae)
  eligible <- zentangler_eligible_sample_ids(mae, c(x_var, perm_vars, dots$covariates, dots$bootstrap_id), views)

  zentangler_log(verbose, sprintf("Negative controls (%s model): observed fit, then %d outcome permutations among %d samples.", model, n_perm, length(eligible)))
  observed_fit <- zentangler_run_fit(mae, x_var, y_var, model, route, as.integer(seed), c(dots, list(verbose = verbose)))
  obs_units <- zentangler_fit_units(observed_fit, q_threshold)
  null_rows <- list()
  null_units <- list()
  null_fits <- list()
  for (i in seq_len(n_perm)) {
    mae_null <- zentangler_permute_outcome(mae, perm_vars, eligible, as.integer(seed) + 100000L + i)
    fit_i <- try(suppressWarnings(zentangler_run_fit(mae_null, x_var, y_var, model, route, as.integer(seed) + i, c(dots, list(verbose = FALSE)))), silent = TRUE)
    if (inherits(fit_i, "try-error")) {
      warning("Permutation ", i, " failed: ", conditionMessage(attr(fit_i, "condition")), call. = FALSE)
      next
    }
    u <- zentangler_fit_units(fit_i, q_threshold)
    null_rows[[length(null_rows) + 1L]] <- cbind(perm_id = i, zentangler_unit_summary(u))
    if (nrow(u)) null_units[[length(null_units) + 1L]] <- cbind(perm_id = i, u)
    if (isTRUE(return_null_fits)) null_fits[[as.character(i)]] <- fit_i
    if (i %% zentangler_progress_every(n_perm) == 0L || i == n_perm) {
      zentangler_log(verbose, sprintf("permutation %d / %d (%d succeeded)", i, n_perm, length(null_rows)), indent = 1L)
    }
  }
  null_summaries <- zentangler_bind_rows_fill(null_rows)
  null_units <- zentangler_bind_rows_fill(null_units)
  obs <- zentangler_unit_summary(obs_units)

  empirical_summary <- data.frame(
    n_null = nrow(null_summaries),
    observed_n_active = obs$n_active,
    null_mean_n_active = if (nrow(null_summaries)) mean(null_summaries$n_active) else NA_real_,
    p_n_active = zentangler_perm_p(null_summaries$n_active, obs$n_active),
    observed_best_q = obs$best_q,
    p_best_q = zentangler_perm_p(null_summaries$best_q, obs$best_q, upper = FALSE),
    observed_max_abs_score = obs$max_abs_score,
    p_max_abs_score = zentangler_perm_p(null_summaries$max_abs_score, obs$max_abs_score),
    stringsAsFactors = FALSE
  )
  empirical_units <- obs_units
  if (nrow(obs_units) && nrow(null_units)) {
    split_null <- split(null_units, null_units$unit_id)
    empirical_units$null_selection_frequency <- vapply(obs_units$unit_id, function(id) {
      nn <- split_null[[id]]
      if (is.null(nn)) 0 else sum(nn$is_active) / nrow(null_summaries)
    }, numeric(1))
    empirical_units$p_abs_score <- vapply(seq_len(nrow(obs_units)), function(k) {
      nn <- split_null[[obs_units$unit_id[k]]]
      zentangler_perm_p(if (is.null(nn)) numeric(0) else nn$abs_score, obs_units$abs_score[k])
    }, numeric(1))
  }

  zentangler_log(verbose, sprintf(
    "Observed %d active (q <= %g) vs %.2f on average under permutation; permutation p = %s.",
    obs$n_active, q_threshold, empirical_summary$null_mean_n_active, zentangler_fmt_p(empirical_summary$p_n_active)
  ))
  out <- list(
    observed_fit = observed_fit,
    observed_summary = obs,
    null_summaries = null_summaries,
    empirical_summary = empirical_summary,
    empirical_units = empirical_units,
    settings = list(model = model, n_perm = n_perm, q_threshold = q_threshold, seed = seed),
    null_fits = if (isTRUE(return_null_fits)) null_fits else NULL
  )
  class(out) <- c("zentangler_negative_control", "list")
  out
}

#' Forward vs reverse direction check
#'
#' Fits `forward_x -> mediators -> forward_y` and the reverse hypothesis
#' `forward_y -> mediators -> forward_x` on the same data. This is a
#' directional sensitivity analysis, not a causal proof. Discoveries are
#' compared with scale-free evidence only (counts and q-values), because
#' forward and reverse scores are in different units.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param forward_x,forward_y Forward exposure and outcome columns.
#' @param model `"parallel"` or `"sequential"`, as in [zentangler()].
#' @param route Character route for the sequential model.
#' @param reverse_route For sequential fits: reuse the route (`"same_order"`)
#'   or reverse it (`"reverse_order"`).
#' @param reverse_y_family Family of the reverse outcome (`forward_x`);
#'   inferred as binomial for 0/1 data, otherwise Gaussian.
#' @param q_threshold q-value cutoff for active units.
#' @param return_fits Keep both fits.
#' @param seed Base seed.
#' @param verbose Log both fits.
#' @param ... Further fit arguments (applied to both directions; outcome-family
#'   arguments apply to the forward fit only).
#' @return A list with `summary` (one row per direction), `comparison` and
#'   optionally `forward_fit` / `reverse_fit`.
#' @export
run_direction_check <- function(mae, forward_x, forward_y, model = c("parallel", "sequential"), route = NULL,
                                reverse_route = c("same_order", "reverse_order"), reverse_y_family = NULL,
                                q_threshold = 0.25, return_fits = FALSE, seed = 1, verbose = TRUE, ...) {
  model <- zentangler_resolve_model(model, route, missing(model))
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  reverse_route <- match.arg(reverse_route)
  dots <- list(...)
  if (identical(dots$y_family, "survival")) stop("Direction checks need a non-survival forward outcome.", call. = FALSE)
  cd <- zentangler_mae_coldata(mae)
  for (v in c(forward_x, forward_y)) if (!(v %in% names(cd))) stop("Column not found in colData: ", v, call. = FALSE)
  if (is.null(reverse_y_family)) {
    vals <- unique(stats::na.omit(suppressWarnings(as.numeric(cd[[forward_x]]))))
    reverse_y_family <- if (length(vals) <= 2L && all(vals %in% c(0, 1))) "binomial" else "gaussian"
  }
  rev_route <- if (!is.null(route) && identical(reverse_route, "reverse_order")) rev(route) else route
  reverse_dots <- dots[setdiff(names(dots), c("y_family", "survival_time_var", "survival_event_var"))]

  zentangler_log(verbose, sprintf("Direction check, forward fit: %s -> mediators -> %s", forward_x, forward_y))
  fit_fwd <- zentangler_run_fit(mae, forward_x, forward_y, model, route, as.integer(seed), c(dots, list(verbose = verbose)))
  zentangler_log(verbose, sprintf("Direction check, reverse fit: %s -> mediators -> %s (%s outcome)", forward_y, forward_x, reverse_y_family))
  fit_rev <- zentangler_run_fit(mae, forward_y, forward_x, model, rev_route, as.integer(seed) + 1L, c(reverse_dots, list(y_family = reverse_y_family, verbose = verbose)))
  fits <- list(forward = fit_fwd, reverse = fit_rev)
  units <- lapply(fits, zentangler_fit_units, q_threshold = q_threshold)
  summary_tab <- do.call(rbind, lapply(names(fits), function(d) {
    eff <- zentangler_effects(fits[[d]])
    u <- units[[d]]
    top <- if (nrow(u)) u$unit_id[which.max(u$abs_score)] else NA_character_
    cbind(
      data.frame(direction = d, x_var = if (d == "forward") forward_x else forward_y,
                 y_var = if (d == "forward") forward_y else forward_x, stringsAsFactors = FALSE),
      zentangler_unit_summary(u),
      data.frame(
        top_unit = top,
        nde = if (nrow(eff)) eff$nde[[1L]] else NA_real_,
        nie_total = if (nrow(eff)) eff$nie_total[[1L]] else NA_real_,
        te = if (nrow(eff)) eff$te[[1L]] else NA_real_,
        prop_mediated = if (nrow(eff)) eff$prop_mediated[[1L]] else NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }))
  f <- summary_tab[1, ]
  r <- summary_tab[2, ]
  wins <- c(
    forward = (f$n_active > r$n_active) + isTRUE(f$best_q < r$best_q),
    reverse = (r$n_active > f$n_active) + isTRUE(r$best_q < f$best_q)
  )
  preference <- if (f$n_active == 0L && r$n_active == 0L) "no_signal" else if (wins[["forward"]] > 0L && wins[["reverse"]] == 0L) "forward_stronger" else if (wins[["reverse"]] > 0L && wins[["forward"]] == 0L) "reverse_stronger" else "mixed"
  top20 <- lapply(units, function(u) utils::head(u$unit_id[order(u$abs_score, decreasing = TRUE)], 20L))
  comparison <- data.frame(
    overlap_active = length(intersect(units$forward$unit_id[units$forward$is_active], units$reverse$unit_id[units$reverse$is_active])),
    overlap_top20 = length(intersect(top20$forward, top20$reverse)),
    direction_preference = preference,
    stringsAsFactors = FALSE
  )
  zentangler_log(verbose, sprintf("Direction check: forward %d vs reverse %d active; preference: %s.",
                                  f$n_active, r$n_active, preference))
  out <- list(
    summary = summary_tab,
    comparison = comparison,
    settings = list(model = model, reverse_route = reverse_route, reverse_y_family = reverse_y_family, q_threshold = q_threshold),
    forward_fit = if (isTRUE(return_fits)) fits$forward else NULL,
    reverse_fit = if (isTRUE(return_fits)) fits$reverse else NULL
  )
  class(out) <- c("zentangler_direction_check", "list")
  out
}
