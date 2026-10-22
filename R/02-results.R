# -----------------------------------------------------------------------------
# Results: accessors, summary() and print() for both fit classes
#
#   zentangler_mediators(fit)  parallel fits: one row per feature
#   zentangler_paths(fit)      sequential fits: one row per retained path
#   zentangler_effects(fit)    both: one-row causal effect summary
#   summary(fit)               model, per-view/threshold counts, top/active
#                              rows, effects and an interpretation
# -----------------------------------------------------------------------------

zentangler_clean_q_thresholds <- function(q_threshold) {
  q_threshold <- sort(unique(as.numeric(q_threshold)))
  q_threshold <- q_threshold[is.finite(q_threshold)]
  if (length(q_threshold) == 0) stop("q_threshold must contain at least one finite numeric cutoff.", call. = FALSE)
  q_threshold
}

zentangler_is_active <- function(tab, q_threshold, q_col = "q_primary") {
  if (!(q_col %in% colnames(tab))) stop("Column not found: ", q_col, call. = FALSE)
  is.finite(tab[[q_col]]) & tab[[q_col]] <= q_threshold & is.finite(tab$b) & tab$b != 0
}

zentangler_filter_rows <- function(tab, q_threshold, top, order_by, q_col) {
  if (is.null(tab) || nrow(tab) == 0L) return(data.frame())
  if (!is.null(q_threshold)) {
    tab <- tab[zentangler_is_active(tab, max(zentangler_clean_q_thresholds(q_threshold)), q_col), , drop = FALSE]
  }
  if (!(order_by %in% colnames(tab))) stop("order_by column not found: ", order_by, call. = FALSE)
  decreasing <- grepl("^abs_|score$|effect$", order_by)
  tab <- tab[order(tab[[order_by]], decreasing = decreasing, na.last = TRUE), , drop = FALSE]
  if (!is.null(top)) tab <- utils::head(tab, as.integer(top))
  rownames(tab) <- NULL
  tab
}

#' Mediator table of a parallel Zentangler fit
#'
#' @param fit A `zentangler_fit`.
#' @param q_threshold If given, keep only mediators with `q_primary` at or
#'   below it (and non-zero `b`).
#' @param top If given, keep only the first `top` rows after ordering.
#' @param order_by Column to order by (default `abs_score`, decreasing;
#'   q- and p-value columns are ordered increasing).
#' @param q_col q-value column used by `q_threshold`.
#' @return A data frame with one row per feature.
#' @export
zentangler_mediators <- function(fit, q_threshold = NULL, top = NULL, order_by = "abs_score", q_col = "q_primary") {
  if (!inherits(fit, "zentangler_fit")) stop("fit must be a zentangler_fit (parallel fit).", call. = FALSE)
  zentangler_filter_rows(fit$mediators, q_threshold, top, order_by, q_col)
}

#' Path table of a sequential Zentangler fit
#'
#' @param fit A `zentangler_sequential_fit`.
#' @inheritParams zentangler_mediators
#' @param order_by Column to order by (default `abs_sequential_score`).
#' @return A data frame with one row per retained path.
#' @export
zentangler_paths <- function(fit, q_threshold = NULL, top = NULL, order_by = "abs_sequential_score", q_col = "q_primary") {
  if (!inherits(fit, "zentangler_sequential_fit")) stop("fit must be a zentangler_sequential_fit.", call. = FALSE)
  zentangler_filter_rows(fit$paths, q_threshold, top, order_by, q_col)
}

#' Causal effect summary of a Zentangler fit
#'
#' @param fit A `zentangler_fit` or `zentangler_sequential_fit`.
#' @return One-row data frame with NDE, NIE, TE, proportion mediated, the
#'   reduced X -> Y effect and (if requested) bootstrap intervals.
#' @export
zentangler_effects <- function(fit) {
  eff <- fit$effects
  if (is.null(eff)) return(data.frame())
  rownames(eff) <- NULL
  eff
}

zentangler_threshold_counts <- function(tab, q_threshold, score_col) {
  if (is.null(tab) || nrow(tab) == 0L) return(data.frame())
  do.call(rbind, lapply(zentangler_clean_q_thresholds(q_threshold), function(q) {
    act <- zentangler_is_active(tab, q)
    data.frame(
      q_threshold = q,
      n_active = sum(act),
      n_positive = sum(act & tab[[score_col]] > 0, na.rm = TRUE),
      n_negative = sum(act & tab[[score_col]] < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

zentangler_interpret_effect_row <- function(eff) {
  if (!is.data.frame(eff) || nrow(eff) == 0L) return(NA_character_)
  nie <- eff$nie_total[[1L]]
  te <- eff$te[[1L]]
  if (!is.finite(nie) || !is.finite(te)) return("Effects are not available.")
  if (abs(nie) < 1e-10) return("The estimated indirect effect is near zero: the significant mediators explain little of the total effect.")
  if (sign(nie) == sign(te)) "Indirect and total effects point in the same direction (same-direction mediation)."
  else "Indirect and total effects point in opposite directions (opposing mediation or suppression)."
}

#' @export
summary.zentangler_fit <- function(object, q_threshold = 0.25, top = 20L, ...) {
  q_main <- max(zentangler_clean_q_thresholds(q_threshold))
  tab <- object$mediators
  s <- object$settings
  views <- intersect(s$view_names %||% unique(tab$omics), unique(tab$omics))
  out <- list(
    model = data.frame(
      n_samples = object$diagnostics$n_samples,
      n_views = length(views),
      n_features = nrow(tab),
      n_screened = sum(tab$selected_by_screen),
      n_active = sum(zentangler_is_active(tab, q_main)),
      q_threshold = q_main,
      y_family = s$y_family,
      fusion = zentangler_fusion_label(s),
      rho = s$rho_used,
      b_inference = s$b_inference,
      fdr_joint = s$fdr_joint_used,
      runtime_seconds = object$diagnostics$runtime_seconds,
      stringsAsFactors = FALSE
    ),
    views = do.call(rbind, lapply(views, function(v) {
      d <- tab[tab$omics == v, , drop = FALSE]
      data.frame(
        omics = v, a_stage_model = object$views[[v]]$a_stage_model %||% "lm",
        n_features = nrow(d), n_screened = sum(d$selected_by_screen),
        n_active = sum(zentangler_is_active(d, q_main)),
        min_q = suppressWarnings(min(d$q_primary, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    })),
    thresholds = zentangler_threshold_counts(tab, q_threshold, "score"),
    top = zentangler_mediators(object, top = top),
    active = zentangler_mediators(object, q_threshold = q_main),
    exposure_outcome = object$diagnostics$exposure_outcome,
    effects = zentangler_effects(object),
    interpretation = zentangler_interpret_effect_row(object$effects)
  )
  class(out) <- c("summary.zentangler_fit", "list")
  out
}

#' @export
summary.zentangler_sequential_fit <- function(object, q_threshold = 0.25, top = 20L, ...) {
  q_main <- max(zentangler_clean_q_thresholds(q_threshold))
  paths <- object$paths
  s <- object$settings
  out <- list(
    model = data.frame(
      route = paste(s$route, collapse = " -> "),
      n_samples = object$diagnostics$n_samples,
      n_paths = nrow(paths),
      n_active_paths = if (nrow(paths) > 0L) sum(zentangler_is_active(paths, q_main)) else 0L,
      q_threshold = q_main,
      y_family = s$y_family,
      fusion = zentangler_fusion_label(s),
      rho = s$rho_used,
      b_inference = s$b_inference,
      runtime_seconds = object$diagnostics$runtime_seconds,
      stringsAsFactors = FALSE
    ),
    thresholds = zentangler_threshold_counts(paths, q_threshold, "sequential_score"),
    top = if (nrow(paths) > 0L) zentangler_paths(object, top = top) else data.frame(),
    active = if (nrow(paths) > 0L) zentangler_paths(object, q_threshold = q_main) else data.frame(),
    edges = object$edges,
    terminals = object$terminals,
    exposure_outcome = object$diagnostics$exposure_outcome,
    effects = zentangler_effects(object),
    interpretation = zentangler_interpret_effect_row(object$effects)
  )
  class(out) <- c("summary.zentangler_sequential_fit", "list")
  out
}

zentangler_fusion_label <- function(s) {
  if (identical(s$fusion, "integrated_learner")) {
    return(paste0("integrated_learner (", s$integrated_learner_fusion, ")"))
  }
  if (is.finite(s$rho_used %||% NA_real_)) sprintf("multiview (rho = %g)", s$rho_used) else "multiview"
}

zentangler_print_exposure_outcome <- function(eo) {
  if (is.null(eo) || !is.data.frame(eo) || nrow(eo) == 0L) return(invisible(NULL))
  cat(sprintf("Exposure-outcome association (%s, Y ~ X%s): coefficient %.4g, p = %s%s\n",
              eo$model, if (eo$n_covariates > 0L) " + covariates" else "", eo$coefficient,
              zentangler_fmt_p(eo$p_value), if (isTRUE(eo$associated_p05)) "" else "  [not significant]"))
}

zentangler_print_effects <- function(eff) {
  if (!is.data.frame(eff) || nrow(eff) == 0L) return(invisible(NULL))
  cat(sprintf(
    "Effects (%s, x %g -> %g): NDE %.4g  NIE %.4g  TE %.4g  prop. mediated %.3g\n",
    eff$effect_scale[[1L]], eff$x0[[1L]], eff$x1[[1L]], eff$nde[[1L]], eff$nie_total[[1L]], eff$te[[1L]], eff$prop_mediated[[1L]]
  ))
}

#' @export
print.zentangler_fit <- function(x, q_threshold = 0.25, ...) {
  s <- x$settings
  tab <- x$mediators
  cat("Zentangler parallel mediation fit\n")
  cat(sprintf("  %d samples, %d views, %d features (%d screened)\n",
              x$diagnostics$n_samples, length(unique(tab$omics)), nrow(tab), sum(tab$selected_by_screen)))
  cat(sprintf("  outcome: %s | outcome model: %s | b: %s | FDR: %s\n",
              s$y_family, zentangler_fusion_label(s), s$b_inference, s$fdr_joint_used))
  cat(sprintf("  active mediators at q <= %g: %d\n", q_threshold, sum(zentangler_is_active(tab, q_threshold))))
  zentangler_print_exposure_outcome(x$diagnostics$exposure_outcome)
  zentangler_print_effects(x$effects)
  invisible(x)
}

#' @export
print.zentangler_sequential_fit <- function(x, q_threshold = 0.25, ...) {
  s <- x$settings
  paths <- x$paths
  cat("Zentangler sequential mediation fit\n")
  cat(sprintf("  route: X -> %s -> Y | %d samples\n", paste(s$route, collapse = " -> "), x$diagnostics$n_samples))
  cat(sprintf("  outcome: %s | terminal model: %s | b: %s\n",
              s$y_family, zentangler_fusion_label(s), s$b_inference))
  n_act <- if (nrow(paths) > 0L) sum(zentangler_is_active(paths, q_threshold)) else 0L
  cat(sprintf("  %d paths, %d active at q <= %g\n", nrow(paths), n_act, q_threshold))
  zentangler_print_exposure_outcome(x$diagnostics$exposure_outcome)
  zentangler_print_effects(x$effects)
  invisible(x)
}

#' @export
print.summary.zentangler_fit <- function(x, ...) {
  cat("Zentangler parallel fit summary\n\n")
  print(x$model, row.names = FALSE)
  cat("\nPer view:\n")
  print(x$views, row.names = FALSE)
  cat("\nActive mediators by q threshold:\n")
  print(x$thresholds, row.names = FALSE)
  if (nrow(x$active) > 0L) {
    cat("\nActive mediators:\n")
    cols <- intersect(c("omics", "mediator", "a", "b", "score", "indirect_effect", "q_primary"), colnames(x$active))
    print(utils::head(x$active[, cols, drop = FALSE], 20L), row.names = FALSE)
  }
  if (!is.null(x$exposure_outcome)) {
    cat("\n")
    zentangler_print_exposure_outcome(x$exposure_outcome)
  }
  if (nrow(x$effects) > 0L) {
    zentangler_print_effects(x$effects)
    cat(x$interpretation, "\n")
  }
  invisible(x)
}

#' @export
print.summary.zentangler_sequential_fit <- function(x, ...) {
  cat("Zentangler sequential fit summary\n\n")
  print(x$model, row.names = FALSE)
  if (nrow(x$thresholds) > 0L) {
    cat("\nActive paths by q threshold:\n")
    print(x$thresholds, row.names = FALSE)
  }
  if (nrow(x$active) > 0L) {
    cat("\nActive paths:\n")
    cols <- intersect(c("mediator_path", "a", "d_product", "b", "sequential_score", "q_primary"), colnames(x$active))
    print(utils::head(x$active[, cols, drop = FALSE], 20L), row.names = FALSE)
  }
  if (!is.null(x$exposure_outcome)) {
    cat("\n")
    zentangler_print_exposure_outcome(x$exposure_outcome)
  }
  if (nrow(x$effects) > 0L) {
    zentangler_print_effects(x$effects)
    cat(x$interpretation, "\n")
  }
  invisible(x)
}
