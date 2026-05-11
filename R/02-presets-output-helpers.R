# -----------------------------------------------------------------------------
# User-facing presets and output helpers
# -----------------------------------------------------------------------------

zentangler_method_presets <- function() {
  data.frame(
    method_preset = c("custom", "fast_lasso", "elastic_net", "longitudinal_maaslin2", "full_exploratory"),
    a_stage_model = c(NA, "lm", "lm", "maaslin2", "lm"),
    screen_method = c(NA, "sis", "sis", "sis", "none"),
    fusion_mode = c(NA, "early", "early", "early", "early"),
    lambda_choice = c(NA, "lambda.1se", "lambda.min", "lambda.1se", "lambda.min"),
    glmnet_alpha = c(NA, 1, 0.5, 1, 0.5),
    b_inference = c(NA, "debiased_lasso", "debiased_lasso", "debiased_lasso", "debiased_lasso"),
    description = c(
      "Do not override individual options.",
      "Fast default: LM A-stage, SIS screening, early-fusion lasso.",
      "Elastic-net power screen: LM A-stage, SIS screening, early fusion, lambda.min, alpha=0.5.",
      "Longitudinal default: MaAsLin2 A-stage with optional random effects, early-fusion lasso.",
      "Broad exploratory mode: LM A-stage, no hard screening, early-fusion elastic net."
    ),
    stringsAsFactors = FALSE
  )
}

zentangler_preset_values <- function(method_preset) {
  method_preset <- match.arg(
    method_preset,
    choices = c("custom", "fast_lasso", "elastic_net", "longitudinal_maaslin2", "full_exploratory")
  )
  switch(
    method_preset,
    custom = list(),
    fast_lasso = list(
      a_stage_model = "lm",
      screen_method = "sis",
      fusion_mode = "early",
      lambda_choice = "lambda.1se",
      glmnet_alpha = 1,
      b_inference = "debiased_lasso"
    ),
    elastic_net = list(
      a_stage_model = "lm",
      screen_method = "sis",
      fusion_mode = "early",
      lambda_choice = "lambda.min",
      glmnet_alpha = 0.5,
      b_inference = "debiased_lasso"
    ),
    longitudinal_maaslin2 = list(
      a_stage_model = "maaslin2",
      screen_method = "sis",
      fusion_mode = "early",
      lambda_choice = "lambda.1se",
      glmnet_alpha = 1,
      b_inference = "debiased_lasso"
    ),
    full_exploratory = list(
      a_stage_model = "lm",
      screen_method = "none",
      fusion_mode = "early",
      lambda_choice = "lambda.min",
      glmnet_alpha = 0.5,
      b_inference = "debiased_lasso"
    )
  )
}

apply_zentangler_preset <- function(method_preset, env) {
  vals <- zentangler_preset_values(method_preset)
  if (length(vals) == 0) return(invisible(NULL))
  for (nm in names(vals)) assign(nm, vals[[nm]], envir = env)
  invisible(NULL)
}

zentangler_all_mediators <- function(fit) {
  tab <- fit$combined_mediators
  if (is.null(tab)) tab <- fit$mediators_all
  if (is.null(tab)) return(data.frame())
  rownames(tab) <- NULL
  tab
}

zentangler_clean_q_thresholds <- function(q_threshold) {
  q_threshold <- sort(unique(as.numeric(q_threshold)))
  q_threshold <- q_threshold[is.finite(q_threshold)]
  if (length(q_threshold) == 0) stop("q_threshold must contain at least one finite numeric cutoff.")
  q_threshold
}

zentangler_active_mediators <- function(fit, q_threshold = 0.25, q_col = "q_primary") {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  tab <- zentangler_all_mediators(fit)
  if (nrow(tab) == 0) return(tab)
  if (!(q_col %in% colnames(tab))) stop("q_col not found in mediator table: ", q_col)
  keep <- is.finite(tab[[q_col]]) & tab[[q_col]] <= q_threshold & is.finite(tab$b) & tab$b != 0
  out <- tab[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

zentangler_top_mediators <- function(
  fit,
  n = 20L,
  order_by = c("abs_score", "q_primary", "p_primary", "score"),
  decreasing = NULL
) {
  order_by <- match.arg(order_by)
  tab <- zentangler_all_mediators(fit)
  if (nrow(tab) == 0) return(tab)
  if (!(order_by %in% colnames(tab))) stop("order_by not found in mediator table: ", order_by)
  if (is.null(decreasing)) decreasing <- order_by %in% c("abs_score", "score")
  ord <- order(tab[[order_by]], decreasing = decreasing, na.last = TRUE)
  out <- tab[ord, , drop = FALSE]
  out <- out[seq_len(min(as.integer(n), nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

zentangler_compute_view_summary <- function(tab, q_threshold = 0.25, q_col = "q_primary") {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  if (is.null(tab) || nrow(tab) == 0) return(data.frame())
  if (!("omics" %in% colnames(tab))) stop("Mediator table must contain an 'omics' column.")
  if (!(q_col %in% colnames(tab))) stop("q_col not found in mediator table: ", q_col)

  views <- unique(as.character(tab$omics))
  out <- lapply(views, function(view) {
    d <- tab[as.character(tab$omics) == view, , drop = FALSE]
    selected <- if ("selected_by_screen" %in% colnames(d)) as.logical(d$selected_by_screen) else rep(NA, nrow(d))
    active <- is.finite(d[[q_col]]) & d[[q_col]] <= q_threshold & is.finite(d$b) & d$b != 0
    data.frame(
      omics = view,
      n_tested = nrow(d),
      n_screened = sum(selected, na.rm = TRUE),
      n_active = sum(active, na.rm = TRUE),
      n_positive_score = sum(active & is.finite(d$score) & d$score > 0, na.rm = TRUE),
      n_negative_score = sum(active & is.finite(d$score) & d$score < 0, na.rm = TRUE),
      max_abs_score = if (any(is.finite(d$abs_score))) max(d$abs_score, na.rm = TRUE) else NA_real_,
      min_q_primary = if (any(is.finite(d[[q_col]]))) min(d[[q_col]], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

zentangler_view_summary <- function(fit, q_threshold = 0.25, q_col = "q_primary") {
  zentangler_compute_view_summary(zentangler_all_mediators(fit), q_threshold = q_threshold, q_col = q_col)
}

zentangler_compute_model_summary <- function(settings, diagnostics, tab, q_threshold = 0.25) {
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))
  active <- if (!is.null(tab) && nrow(tab) > 0) {
    is.finite(tab$q_primary) & tab$q_primary <= q_threshold & is.finite(tab$b) & tab$b != 0
  } else {
    logical(0)
  }
  data.frame(
    method_preset = settings$method_preset %||% NA_character_,
    input_container = settings$input_container %||% NA_character_,
    n_samples = diagnostics$n_samples %||% NA_integer_,
    n_views = settings$n_views %||% NA_integer_,
    n_features_after_filtering = sum(diagnostics$features_after_filtering %||% 0L),
    n_screened = sum(diagnostics$screened_counts %||% 0L),
    n_active_q025 = sum(active, na.rm = TRUE),
    a_stage_model = settings$a_stage_model %||% NA_character_,
    fusion_mode = settings$fusion_mode %||% NA_character_,
    lambda_choice = settings$lambda_choice %||% NA_character_,
    glmnet_alpha = settings$glmnet_alpha %||% NA_real_,
    fdr_method = settings$fdr_method %||% NA_character_,
    fdr_scope = settings$fdr_scope %||% NA_character_,
    primary_inference = settings$primary_inference %||% NA_character_,
    b_inference = settings$b_inference %||% NA_character_,
    b_inference_method = settings$b_inference_method %||% NA_character_,
    runtime_seconds = diagnostics$runtime_seconds %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

zentangler_model_summary <- function(fit, q_threshold = 0.25) {
  zentangler_compute_model_summary(
    settings = fit$settings %||% list(),
    diagnostics = fit$diagnostics %||% list(),
    tab = zentangler_all_mediators(fit),
    q_threshold = q_threshold
  )
}

zentangler_threshold_summary <- function(
  fit,
  q_threshold = seq(0.05, 0.25, by = 0.05),
  q_cols = NULL
) {
  tab <- zentangler_all_mediators(fit)
  if (nrow(tab) == 0) return(data.frame())
  q_threshold <- zentangler_clean_q_thresholds(q_threshold)
  if (is.null(q_cols)) {
    q_cols <- intersect(c("q_primary", "q_primary_global", "q_primary_within_view"), colnames(tab))
  }
  q_cols <- intersect(q_cols, colnames(tab))
  if (length(q_cols) == 0) stop("No valid q-value columns found in mediator table.")

  rows <- list()
  ii <- 0L
  for (q_col in q_cols) {
    for (threshold in q_threshold) {
      active <- is.finite(tab[[q_col]]) & tab[[q_col]] <= threshold & is.finite(tab$b) & tab$b != 0
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        q_col = q_col,
        q_threshold = threshold,
        n_active = sum(active, na.rm = TRUE),
        n_positive_score = sum(active & is.finite(tab$score) & tab$score > 0, na.rm = TRUE),
        n_negative_score = sum(active & is.finite(tab$score) & tab$score < 0, na.rm = TRUE),
        min_q = if (any(is.finite(tab[[q_col]]))) min(tab[[q_col]], na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_zentangler <- function(fit, q_threshold = 0.25) {
  list(
    model_summary = zentangler_model_summary(fit, q_threshold = q_threshold),
    view_summary = zentangler_view_summary(fit, q_threshold = q_threshold),
    threshold_summary = zentangler_threshold_summary(fit, q_threshold = q_threshold),
    top_mediators = zentangler_top_mediators(fit, n = 20L),
    active_mediators = zentangler_active_mediators(fit, q_threshold = q_threshold),
    causal_effects = fit$causal_effects %||% fit$effect_decomposition
  )
}

