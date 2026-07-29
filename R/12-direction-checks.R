# -----------------------------------------------------------------------------
# Directional sensitivity / reverse-hypothesis checks
# -----------------------------------------------------------------------------

zentangler_infer_reverse_y_family <- function(y) {
  y <- suppressWarnings(as.numeric(y))
  y <- y[is.finite(y)]
  if (length(y) == 0L) stop("Reverse outcome contains no finite numeric values.", call. = FALSE)
  vals <- sort(unique(y))
  if (length(vals) <= 2L && all(vals %in% c(0, 1))) return("binomial")
  "gaussian"
}

zentangler_parallel_mediator_ids <- function(tab) {
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(character(0))
  if (all(c("omics", "mediator") %in% names(tab))) {
    return(paste(tab$omics, tab$mediator, sep = "::"))
  }
  if ("mediator" %in% names(tab)) return(as.character(tab$mediator))
  character(0)
}

zentangler_sequential_path_ids <- function(tab) {
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(character(0))
  if ("key_path" %in% names(tab)) return(as.character(tab$key_path))
  if ("mediator_path" %in% names(tab)) return(as.character(tab$mediator_path))
  character(0)
}

zentangler_safe_first <- function(x) {
  if (length(x) == 0L) return(NA_character_)
  x[[1L]]
}

zentangler_overlap_count <- function(a, b) {
  length(intersect(unique(as.character(a)), unique(as.character(b))))
}

zentangler_direction_preference <- function(
  forward_n,
  reverse_n,
  forward_best_q,
  reverse_best_q,
  forward_score,
  reverse_score
) {
  forward_signal <- is.finite(forward_n) && forward_n > 0L
  reverse_signal <- is.finite(reverse_n) && reverse_n > 0L
  if (!forward_signal && !reverse_signal) return("no_signal")

  forward_wins <- 0L
  reverse_wins <- 0L

  if (is.finite(forward_n) && is.finite(reverse_n)) {
    if (forward_n > reverse_n) forward_wins <- forward_wins + 1L
    if (reverse_n > forward_n) reverse_wins <- reverse_wins + 1L
  }
  if (is.finite(forward_best_q) && is.finite(reverse_best_q)) {
    if (forward_best_q < reverse_best_q) forward_wins <- forward_wins + 1L
    if (reverse_best_q < forward_best_q) reverse_wins <- reverse_wins + 1L
  }
  if (is.finite(forward_score) && is.finite(reverse_score)) {
    if (forward_score > reverse_score) forward_wins <- forward_wins + 1L
    if (reverse_score > forward_score) reverse_wins <- reverse_wins + 1L
  }

  if (forward_wins > 0L && reverse_wins == 0L) return("forward_stronger")
  if (reverse_wins > 0L && forward_wins == 0L) return("reverse_stronger")
  if (forward_wins == 0L && reverse_wins == 0L) return("mixed")
  "mixed"
}

zentangler_parallel_direction_summary <- function(
  fit,
  direction,
  x_var,
  y_var,
  y_family,
  q_threshold = 0.25
) {
  tab <- zentangler_all_mediators(fit)
  active <- zentangler_active_mediators(fit, q_threshold = q_threshold)
  top <- zentangler_top_mediators(fit, n = 1L)
  data.frame(
    direction = as.character(direction),
    x_var = as.character(x_var),
    y_var = as.character(y_var),
    y_family = as.character(y_family),
    n_mediators_tested = nrow(tab),
    n_active_mediators = nrow(active),
    best_q = zentangler_min_finite_or_na(zentangler_get_primary_q_values(tab)),
    max_abs_score = if ("abs_score" %in% names(tab)) zentangler_max_finite_or_na(tab$abs_score) else if ("score" %in% names(tab)) zentangler_max_finite_or_na(abs(tab$score)) else NA_real_,
    top_mediator = if (nrow(top) > 0L && "mediator" %in% names(top)) as.character(top$mediator[[1L]]) else NA_character_,
    top_view = if (nrow(top) > 0L && "omics" %in% names(top)) as.character(top$omics[[1L]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

zentangler_sequential_direction_summary <- function(
  fit,
  direction,
  x_var,
  y_var,
  route,
  q_threshold = 0.25
) {
  tab <- zentangler_sequential_paths(fit)
  active <- zentangler_sequential_active_paths(fit, q_threshold = q_threshold)
  top <- zentangler_sequential_top_paths(fit, n = 1L)
  data.frame(
    direction = as.character(direction),
    x_var = as.character(x_var),
    y_var = as.character(y_var),
    route = paste(as.character(route), collapse = " -> "),
    n_paths = nrow(tab),
    n_significant_paths = nrow(active),
    best_q = zentangler_min_finite_or_na(zentangler_get_primary_q_values(tab)),
    max_abs_sequential_score = if ("abs_sequential_score" %in% names(tab)) zentangler_max_finite_or_na(tab$abs_sequential_score) else NA_real_,
    top_path = if (nrow(top) > 0L) zentangler_safe_first(zentangler_sequential_path_ids(top)) else NA_character_,
    top_terminal_view = if (nrow(top) > 0L && "terminal_view" %in% names(top)) as.character(top$terminal_view[[1L]]) else NA_character_,
    top_terminal_mediator = if (nrow(top) > 0L && "terminal_mediator" %in% names(top)) as.character(top$terminal_mediator[[1L]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

zentangler_parallel_direction_comparison <- function(forward_fit, reverse_fit, forward_summary, reverse_summary, q_threshold = 0.25) {
  active_forward <- zentangler_parallel_mediator_ids(zentangler_active_mediators(forward_fit, q_threshold = q_threshold))
  active_reverse <- zentangler_parallel_mediator_ids(zentangler_active_mediators(reverse_fit, q_threshold = q_threshold))
  top_forward <- zentangler_parallel_mediator_ids(zentangler_top_mediators(forward_fit, n = 20L))
  top_reverse <- zentangler_parallel_mediator_ids(zentangler_top_mediators(reverse_fit, n = 20L))
  data.frame(
    forward_x = forward_summary$x_var[[1L]],
    forward_y = forward_summary$y_var[[1L]],
    reverse_x = reverse_summary$x_var[[1L]],
    reverse_y = reverse_summary$y_var[[1L]],
    forward_n_active_mediators = forward_summary$n_active_mediators[[1L]],
    reverse_n_active_mediators = reverse_summary$n_active_mediators[[1L]],
    forward_best_q = forward_summary$best_q[[1L]],
    reverse_best_q = reverse_summary$best_q[[1L]],
    forward_max_abs_score = forward_summary$max_abs_score[[1L]],
    reverse_max_abs_score = reverse_summary$max_abs_score[[1L]],
    overlap_active_mediators = zentangler_overlap_count(active_forward, active_reverse),
    overlap_top20_mediators = zentangler_overlap_count(top_forward, top_reverse),
    direction_preference = zentangler_direction_preference(
      forward_n = forward_summary$n_active_mediators[[1L]],
      reverse_n = reverse_summary$n_active_mediators[[1L]],
      forward_best_q = forward_summary$best_q[[1L]],
      reverse_best_q = reverse_summary$best_q[[1L]],
      forward_score = forward_summary$max_abs_score[[1L]],
      reverse_score = reverse_summary$max_abs_score[[1L]]
    ),
    stringsAsFactors = FALSE
  )
}

zentangler_sequential_direction_comparison <- function(
  forward_fit,
  reverse_fit,
  forward_summary,
  reverse_summary,
  q_threshold = 0.25,
  reverse_route_mode,
  reverse_route
) {
  active_forward <- zentangler_sequential_path_ids(zentangler_sequential_active_paths(forward_fit, q_threshold = q_threshold))
  active_reverse <- zentangler_sequential_path_ids(zentangler_sequential_active_paths(reverse_fit, q_threshold = q_threshold))
  top_forward <- zentangler_sequential_path_ids(zentangler_sequential_top_paths(forward_fit, n = 20L))
  top_reverse <- zentangler_sequential_path_ids(zentangler_sequential_top_paths(reverse_fit, n = 20L))
  terminals_forward <- zentangler_sequential_terminals(forward_fit)
  terminals_reverse <- zentangler_sequential_terminals(reverse_fit)
  term_forward_ids <- if (nrow(terminals_forward) > 0L && "terminal_key" %in% names(terminals_forward)) as.character(terminals_forward$terminal_key) else character(0)
  term_reverse_ids <- if (nrow(terminals_reverse) > 0L && "terminal_key" %in% names(terminals_reverse)) as.character(terminals_reverse$terminal_key) else character(0)
  data.frame(
    forward_x = forward_summary$x_var[[1L]],
    forward_y = forward_summary$y_var[[1L]],
    forward_route = forward_summary$route[[1L]],
    reverse_route_mode = as.character(reverse_route_mode),
    reverse_route = paste(as.character(reverse_route), collapse = " -> "),
    forward_n_significant_paths = forward_summary$n_significant_paths[[1L]],
    reverse_n_significant_paths = reverse_summary$n_significant_paths[[1L]],
    forward_best_q = forward_summary$best_q[[1L]],
    reverse_best_q = reverse_summary$best_q[[1L]],
    forward_max_abs_sequential_score = forward_summary$max_abs_sequential_score[[1L]],
    reverse_max_abs_sequential_score = reverse_summary$max_abs_sequential_score[[1L]],
    overlap_terminal_mediators = zentangler_overlap_count(term_forward_ids, term_reverse_ids),
    overlap_top_paths = zentangler_overlap_count(top_forward, top_reverse),
    overlap_active_paths = zentangler_overlap_count(active_forward, active_reverse),
    direction_preference = zentangler_direction_preference(
      forward_n = forward_summary$n_significant_paths[[1L]],
      reverse_n = reverse_summary$n_significant_paths[[1L]],
      forward_best_q = forward_summary$best_q[[1L]],
      reverse_best_q = reverse_summary$best_q[[1L]],
      forward_score = forward_summary$max_abs_sequential_score[[1L]],
      reverse_score = reverse_summary$max_abs_sequential_score[[1L]]
    ),
    stringsAsFactors = FALSE
  )
}

#' Run a bidirectional direction check for parallel Zentangler
#'
#' Fits the forward hypothesis `forward_x -> mediators -> forward_y`, then fits
#' the reverse hypothesis `forward_y -> mediators -> forward_x` on the same MAE.
#' The result is a directional sensitivity analysis, not a causal proof.
#'
#' Phase 1 currently supports Gaussian and binomial reverse outcomes only.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param forward_x Forward exposure column in `colData(mae)`.
#' @param forward_y Forward outcome column in `colData(mae)`.
#' @param reverse_y_family Optional reverse outcome family. If `NULL`, it is
#'   inferred from the forward exposure column as `"binomial"` for 0/1 data or
#'   `"gaussian"` otherwise.
#' @param q_threshold Q-value cutoff used for active mediator overlap summaries.
#' @param return_fits Logical. If `TRUE`, return the full forward and reverse
#'   fitted objects.
#' @param seed Base seed for the two fits.
#' @param ... Additional arguments passed to `fit_multiview_parallel_zentangler()`.
#'
#' @return A list containing forward and reverse summaries, a comparison table,
#'   and optionally the two fitted objects.
#' @export
run_parallel_direction_check <- function(
  mae,
  forward_x,
  forward_y,
  reverse_y_family = NULL,
  q_threshold = 0.25,
  return_fits = FALSE,
  seed = 1,
  ...
) {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  dots <- list(...)
  forward_y_family <- dots$y_family %||% "gaussian"
  if (identical(forward_y_family, "survival")) {
    stop("run_parallel_direction_check() currently does not support survival outcomes.", call. = FALSE)
  }
  cd <- zentangler_mae_coldata(mae)
  if (!(forward_x %in% names(cd))) stop("forward_x not found in MAE colData: ", forward_x, call. = FALSE)
  if (!(forward_y %in% names(cd))) stop("forward_y not found in MAE colData: ", forward_y, call. = FALSE)
  reverse_y_family <- reverse_y_family %||% zentangler_infer_reverse_y_family(cd[[forward_x]])

  forward_fit <- fit_multiview_parallel_zentangler(
    mae = mae,
    x_var = forward_x,
    y_var = forward_y,
    seed = as.integer(seed),
    ...
  )
  reverse_fit <- fit_multiview_parallel_zentangler(
    mae = mae,
    x_var = forward_y,
    y_var = forward_x,
    y_family = reverse_y_family,
    seed = as.integer(seed) + 1L,
    ...
  )

  forward_summary <- zentangler_parallel_direction_summary(
    fit = forward_fit,
    direction = "forward",
    x_var = forward_x,
    y_var = forward_y,
    y_family = forward_y_family,
    q_threshold = q_threshold
  )
  reverse_summary <- zentangler_parallel_direction_summary(
    fit = reverse_fit,
    direction = "reverse",
    x_var = forward_y,
    y_var = forward_x,
    y_family = reverse_y_family,
    q_threshold = q_threshold
  )
  comparison_summary <- zentangler_parallel_direction_comparison(
    forward_fit = forward_fit,
    reverse_fit = reverse_fit,
    forward_summary = forward_summary,
    reverse_summary = reverse_summary,
    q_threshold = q_threshold
  )

  out <- list(
    forward_summary = forward_summary,
    reverse_summary = reverse_summary,
    comparison_summary = comparison_summary,
    settings = list(
      model_type = "parallel",
      q_threshold = q_threshold,
      seed = as.integer(seed),
      return_fits = isTRUE(return_fits)
    )
  )
  if (isTRUE(return_fits)) {
    out$forward_fit <- forward_fit
    out$reverse_fit <- reverse_fit
  }
  class(out) <- c("zentangler_parallel_direction_check", "list")
  out
}

#' Run a bidirectional direction check for strict-route sequential Zentangler
#'
#' Fits the forward hypothesis `forward_x -> route -> forward_y`, then fits a
#' reverse hypothesis `forward_y -> route -> forward_x` using either the same
#' mediator order or the reversed route order. The result is a directional
#' sensitivity analysis, not a causal proof.
#'
#' Phase 1 currently supports Gaussian and binomial reverse outcomes only.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param forward_x Forward exposure column in `colData(mae)`.
#' @param forward_y Forward outcome column in `colData(mae)`.
#' @param forward_route Character vector of view names defining the forward
#'   route.
#' @param reverse_route Whether the reverse fit should use the same mediator
#'   order or the reversed order.
#' @param reverse_y_family Optional reverse outcome family. If `NULL`, it is
#'   inferred from the forward exposure column as `"binomial"` for 0/1 data or
#'   `"gaussian"` otherwise.
#' @param q_threshold Q-value cutoff used for active path overlap summaries.
#' @param return_fits Logical. If `TRUE`, return the full forward and reverse
#'   fitted objects.
#' @param seed Base seed for the two fits.
#' @param ... Additional arguments passed to `fit_sequential_zentangler()`.
#'
#' @return A list containing forward and reverse summaries, a comparison table,
#'   and optionally the two fitted objects.
#' @export
run_sequential_direction_check <- function(
  mae,
  forward_x,
  forward_y,
  forward_route,
  reverse_route = c("same_order", "reverse_order"),
  reverse_y_family = NULL,
  q_threshold = 0.25,
  return_fits = FALSE,
  seed = 1,
  ...
) {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  reverse_route <- match.arg(reverse_route)
  dots <- list(...)
  forward_y_family <- dots$y_family %||% "gaussian"
  if (identical(forward_y_family, "survival")) {
    stop("run_sequential_direction_check() currently does not support survival outcomes.", call. = FALSE)
  }
  forward_route <- as.character(forward_route)
  forward_route <- forward_route[!is.na(forward_route) & nzchar(forward_route)]
  if (length(forward_route) < 2L) {
    stop("forward_route must contain at least two views.", call. = FALSE)
  }
  cd <- zentangler_mae_coldata(mae)
  if (!(forward_x %in% names(cd))) stop("forward_x not found in MAE colData: ", forward_x, call. = FALSE)
  if (!(forward_y %in% names(cd))) stop("forward_y not found in MAE colData: ", forward_y, call. = FALSE)
  reverse_y_family <- reverse_y_family %||% zentangler_infer_reverse_y_family(cd[[forward_x]])
  reverse_route_views <- if (identical(reverse_route, "reverse_order")) rev(forward_route) else forward_route

  forward_fit <- fit_sequential_zentangler(
    mae = mae,
    x_var = forward_x,
    y_var = forward_y,
    path_templates = list(route = forward_route),
    seed = as.integer(seed),
    ...
  )
  reverse_fit <- fit_sequential_zentangler(
    mae = mae,
    x_var = forward_y,
    y_var = forward_x,
    path_templates = list(route = reverse_route_views),
    y_family = reverse_y_family,
    seed = as.integer(seed) + 1L,
    ...
  )

  forward_summary <- zentangler_sequential_direction_summary(
    fit = forward_fit,
    direction = "forward",
    x_var = forward_x,
    y_var = forward_y,
    route = forward_route,
    q_threshold = q_threshold
  )
  reverse_summary <- zentangler_sequential_direction_summary(
    fit = reverse_fit,
    direction = "reverse",
    x_var = forward_y,
    y_var = forward_x,
    route = reverse_route_views,
    q_threshold = q_threshold
  )
  comparison_summary <- zentangler_sequential_direction_comparison(
    forward_fit = forward_fit,
    reverse_fit = reverse_fit,
    forward_summary = forward_summary,
    reverse_summary = reverse_summary,
    q_threshold = q_threshold,
    reverse_route_mode = reverse_route,
    reverse_route = reverse_route_views
  )

  out <- list(
    forward_summary = forward_summary,
    reverse_summary = reverse_summary,
    comparison_summary = comparison_summary,
    settings = list(
      model_type = "sequential",
      reverse_route_mode = reverse_route,
      q_threshold = q_threshold,
      seed = as.integer(seed),
      return_fits = isTRUE(return_fits)
    )
  )
  if (isTRUE(return_fits)) {
    out$forward_fit <- forward_fit
    out$reverse_fit <- reverse_fit
  }
  class(out) <- c("zentangler_sequential_direction_check", "list")
  out
}
