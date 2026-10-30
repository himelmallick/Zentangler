# -----------------------------------------------------------------------------
# Parallel multi-view mediation: X -> M_j (view k) -> Y
#
# fit_multiview_parallel_zentangler() validates its arguments once into an
# options list, extracts the MultiAssayExperiment into sample-by-feature
# blocks, and runs zentangler_parallel_engine(). Bootstrap replicates call the
# same engine with a modified copy of the options.
# -----------------------------------------------------------------------------

zentangler_validate_binary_outcome <- function(y) {
  y <- as.numeric(y)
  vals <- sort(unique(y[is.finite(y)]))
  if (length(vals) != 2L || !all(vals %in% c(0, 1))) {
    stop(
      "y_family = 'binomial' requires the outcome to be coded 0/1. Observed values: ",
      paste(utils::head(vals, 10), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

zentangler_require_string <- function(x, label) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))) {
    stop(label, " is required.", call. = FALSE)
  }
  as.character(x)
}

zentangler_parallel_options <- function(a) {
  # `a` is the named argument list of fit_multiview_parallel_zentangler().
  o <- a
  o$y_family <- match.arg(a$y_family, c("gaussian", "binomial", "survival"))
  o$study_design <- validate_study_design(a$study_design)
  o$exposure_role <- validate_exposure_role(a$exposure_role)
  o$sis_rank <- match.arg(a$sis_rank, c("abs_a", "pvalue"))
  o$screen_method <- match.arg(a$screen_method, c("sis", "none"))
  # a_stage_model / maaslin2_args may be per view; they are resolved by
  # zentangler_resolve_a_stage() once the views are known.
  o$fusion <- match.arg(a$fusion, c("multiview", "integrated_learner"))
  o$integrated_learner_fusion <- match.arg(a$integrated_learner_fusion, c("late", "early"))
  o$engine <- zentangler_resolve_fusion(o$fusion, o$integrated_learner_fusion)
  o$lambda_choice <- match.arg(a$lambda_choice, c("lambda.1se", "lambda.min"))
  o$glmnet_alpha <- validate_glmnet_alpha(a$glmnet_alpha)
  o$b_inference <- zentangler_validate_b_inference(a$b_inference)
  o$fdr_joint <- validate_fdr_joint(a$fdr_joint)
  o$fdr_method <- validate_fdr_method(a$fdr_method)
  o$fdr_scope <- validate_fdr_scope(a$fdr_scope)
  o$primary_inference <- validate_primary_inference(a$primary_inference)
  o$causal_inference <- match.arg(a$causal_inference, c("none", "bootstrap"))
  o$effect_q_threshold <- as.numeric(a$effect_q_threshold)
  o$bootstrap_repeats <- max(0L, as.integer(a$bootstrap_repeats))
  o$debias_max_targets <- max(1L, as.integer(a$debias_max_targets))
  o$maaslin2_args <- as.list(a$maaslin2_args %||% list())
  o$tweedieverse_args <- as.list(a$tweedieverse_args %||% list())
  if (identical(o$y_family, "survival")) {
    o$survival_time_var <- zentangler_require_string(a$survival_time_var, "survival_time_var (y_family = 'survival')")
    o$survival_event_var <- zentangler_require_string(a$survival_event_var, "survival_event_var (y_family = 'survival')")
    o$y_var <- NULL
  } else {
    o$y_var <- zentangler_require_string(a$y_var, "y_var")
    o$survival_time_var <- NULL
    o$survival_event_var <- NULL
  }
  needs_boot <- c(
    b_inference = identical(o$b_inference, "bootstrap"),
    primary_inference = identical(o$primary_inference, "bootstrap_score"),
    causal_inference = identical(o$causal_inference, "bootstrap")
  )
  if (any(needs_boot) && o$bootstrap_repeats < 1L) {
    stop(paste(names(needs_boot)[needs_boot], collapse = ", "), " = bootstrap requires bootstrap_repeats > 0.", call. = FALSE)
  }
  if (is.null(o$bootstrap_seed)) o$bootstrap_seed <- as.integer(o$seed) + 10000L
  o
}

#' Fit multi-view parallel Zentangler mediation
#'
#' Estimates exposure -> mediator -> outcome paths for every feature of every
#' view in a `MultiAssayExperiment`. The A-path (`M ~ X + C`) is estimated for
#' all features, mediators are screened per view, the outcome model
#' `Y ~ X + C + screened mediators` is fitted with multiview cooperative
#' learning (`rho` in \[0, 1\] chosen by cross-validation), the B-path is
#' de-biased or bootstrapped, and the joint hypothesis `a = 0 or b = 0` is
#' tested with a composite-null FDR step. Causal direct/indirect/total effects
#' are estimated for the significant mediators.
#'
#' @param mae A `MultiAssayExperiment`; phenotype columns live in `colData(mae)`.
#' @param x_var Numeric exposure column (`study_design = "standard"`).
#' @param y_var Outcome column (Gaussian or 0/1 binomial outcomes).
#' @param covariates Optional covariate columns.
#' @param y_family `"gaussian"`, `"binomial"` or `"survival"`.
#' @param survival_time_var,survival_event_var Survival time and 0/1 event columns.
#' @param study_design `"standard"`, `"case_control"`, `"time"` or
#'   `"case_control_time"`; non-standard designs derive a 0/1 exposure.
#' @param case_var,case_level,control_level Case/control column and levels.
#' @param time_var,time_ref,time_compare Time column and the two levels compared.
#' @param exposure_role For `"case_control_time"`: `"time"`, `"case"` or `"interaction"`.
#' @param add_design_covariates Add the non-exposure design indicators as covariates.
#' @param view_names,assay_names Optional subset of MAE experiments and assay names.
#' @param duplicate_primary How to combine several assay columns mapping to one
#'   primary sample: `"mean"` or `"first"`.
#' @param fill_nonfinite_zero Replace non-finite assay values by 0.
#' @param transform Feature transformation per view: `"none"`, `"log"` or
#'   `"clr"` (centered log-ratio). A single value applies to every view; a
#'   vector named by view sets each view, e.g.
#'   `c(species = "clr", metabolites = "log")`.
#' @param pseudocount Pseudocount for `"log"` / `"clr"` (single value or named
#'   by view); default half the smallest positive value of each view.
#' @param a_stage_model A-path (exposure -> mediator) model: `"lm"`,
#'   `"maaslin2"` or `"tweedieverse"`. One value applies to every view; a
#'   vector named by view sets each view, e.g.
#'   `c(species = "tweedieverse", metabolites = "lm")` (views not named use
#'   `"lm"`). Only the A-stage (screening and `p_a`) is affected; the B-path
#'   and causal effects always use the (transformed) feature scale.
#'   `"tweedieverse"` fits per-feature Tweedie GLMs to the untransformed
#'   values of the view (even when `transform` is `"log"` or `"clr"`): counts,
#'   relative abundances or continuous data, zero-inflated or not. For
#'   non-negative data `a` is a log fold change per unit of the exposure;
#'   data with negative values are fitted with the Gaussian case (`p = 0`,
#'   identity link).
#' @param maaslin2_args Arguments passed to `Maaslin2::Maaslin2()` for the
#'   MaAsLin2 views: any of its arguments (e.g. `random_effects`,
#'   `normalization`, `transform`, `analysis_method`, `reference`,
#'   `min_prevalence`, `output`). Either one list shared by all MaAsLin2 views,
#'   or a list named by view whose elements are such lists, e.g.
#'   `list(species = list(transform = "LOG"), kos = list(normalization = "TSS"))`.
#'   Zentangler supplies the data and adds the exposure and covariates to
#'   `fixed_effects`; its own defaults are no normalization/transformation, LM,
#'   and no prevalence/abundance filtering.
#' @param tweedieverse_args Arguments passed to `Tweedieverse::Tweedieverse()`
#'   for the Tweedieverse views (e.g. `random_effects`, `normalization`,
#'   `domain`, `adjust_offset`, `tweedie_p`, `link`, `reference`, `BPPARAM`,
#'   `output`); one shared list or a list named by view, as for
#'   `maaslin2_args`. Zentangler supplies the data and adds the exposure and
#'   covariates to `fixed_effects`; its own defaults turn off Tweedieverse's
#'   prevalence/abundance filtering, metadata standardization, plots and
#'   file output.
#' @param sis_n Mediators kept per view by the screen (default `n / log(n)`).
#' @param sis_rank Screen ranking: standardized `|a|` (`"abs_a"`) or `"pvalue"`.
#' @param screen_method `"sis"` or `"none"`.
#' @param fusion Outcome model: `"multiview"` (cooperative learning, default)
#'   or the supplementary `"integrated_learner"`.
#' @param rho Fixed multiview agreement parameter in \[0, 1\] (0 = early,
#'   1 = late fusion); `NULL` (default) tunes it by cross-validation.
#' @param rho_grid Candidate `rho` values for cross-validation.
#' @param integrated_learner_fusion For `fusion = "integrated_learner"`:
#'   `"late"` (per-view models stacked) or `"early"` (one concatenated model).
#' @param lambda_choice `"lambda.1se"` or `"lambda.min"`.
#' @param glmnet_alpha Elastic-net mixing (1 = lasso).
#' @param b_inference `"debiased"` or `"bootstrap"`.
#' @param debias_max_targets Maximum screened mediators to de-bias.
#' @param fdr_joint `"hdmt"` (composite-null FDR, default) or `"maxp"`.
#' @param fdr_method `"BH"` or `"BY"` (used by `"maxp"` and bootstrap-score p-values).
#' @param fdr_scope `"global"` or `"within_view"`.
#' @param primary_inference `"model_based"` or `"bootstrap_score"`.
#' @param causal_inference `"none"` or `"bootstrap"` (fixed-set effect intervals).
#' @param effect_q_threshold Mediators with `q_primary` at or below this enter
#'   the causal effect stage.
#' @param effect_x0,effect_x1 Exposure contrast for causal effects (default:
#'   0 vs 1 for binary exposures, otherwise the quartiles).
#' @param bootstrap_repeats,bootstrap_ci_level,bootstrap_seed,bootstrap_id
#'   Bootstrap replicates, interval level, seed and optional cluster column.
#' @param seed Random seed (the caller's RNG state is restored on exit).
#' @param return_fits Keep the fitted outcome model in `fit$model`.
#' @param verbose Print timestamped progress messages (data, exposure-outcome
#'   association, A-stage, outcome model, inference, FDR, effects, bootstrap
#'   progress). Messages go through `message()`.
#'
#' @return An object of class `zentangler_fit`: `settings`, `diagnostics`,
#'   `mediators` (one row per feature), `views` (per-view screens), `effects`
#'   (one-row causal effect summary), `bootstrap` and optionally `model`. Use
#'   [zentangler_mediators()], [zentangler_effects()] and `summary()`.
#' @export
fit_multiview_parallel_zentangler <- function(
  mae,
  x_var = NULL,
  y_var = NULL,
  covariates = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  survival_time_var = NULL,
  survival_event_var = NULL,
  study_design = c("standard", "case_control", "time", "case_control_time"),
  case_var = NULL,
  case_level = NULL,
  control_level = NULL,
  time_var = NULL,
  time_ref = NULL,
  time_compare = NULL,
  exposure_role = c("time", "case", "interaction"),
  add_design_covariates = TRUE,
  view_names = NULL,
  assay_names = NULL,
  transform = "none",
  pseudocount = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE,
  a_stage_model = "lm",
  maaslin2_args = list(),
  tweedieverse_args = list(),
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  fusion = c("multiview", "integrated_learner"),
  rho = NULL,
  rho_grid = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  integrated_learner_fusion = c("late", "early"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  b_inference = c("debiased", "bootstrap"),
  debias_max_targets = 200L,
  fdr_joint = c("hdmt", "maxp"),
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global", "within_view"),
  primary_inference = c("model_based", "bootstrap_score"),
  causal_inference = c("none", "bootstrap"),
  effect_q_threshold = 0.25,
  effect_x0 = NULL,
  effect_x1 = NULL,
  bootstrap_repeats = 0L,
  bootstrap_ci_level = 0.95,
  bootstrap_seed = NULL,
  bootstrap_id = NULL,
  seed = 1,
  return_fits = FALSE,
  verbose = TRUE
) {
  args <- mget(setdiff(names(formals()), "mae"), envir = environment())
  opts <- zentangler_parallel_options(args)
  inputs <- zentangler_mae_to_blocks(
    mae = mae,
    view_names = view_names,
    assay_names = assay_names,
    duplicate_primary = match.arg(duplicate_primary),
    fill_nonfinite_zero = fill_nonfinite_zero,
    transform = transform,
    pseudocount = pseudocount
  )
  opts$a_stage <- zentangler_resolve_a_stage(opts$a_stage_model, opts$maaslin2_args, names(inputs$blocks), inputs$view_map, opts$tweedieverse_args)
  if (isTRUE(verbose)) {
    zentangler_log(TRUE, "Parallel mediation: ", length(inputs$blocks), " view(s) from the MultiAssayExperiment.")
    tr <- inputs$view_map$transform
    if (any(tr != "none")) {
      zentangler_log(TRUE, "Transforms: ", paste(paste0(inputs$view_map$view, " = ", tr)[tr != "none"], collapse = ", "), indent = 1L)
    }
  }
  blocks <- zentangler_attach_raw_blocks(inputs$blocks, inputs$raw_blocks, opts$a_stage)
  fit <- zentangler_muffle_cox_ties(zentangler_parallel_engine(blocks, inputs$pheno_df, opts))
  fit$settings$mae_view_map <- inputs$view_map
  fit
}

zentangler_parallel_engine <- function(blocks, pheno_df, opts, skip_effects = FALSE) {
  start_time <- Sys.time()
  restore_rng <- zentangler_save_rng()
  on.exit(restore_rng(), add = TRUE)
  set.seed(opts$seed)
  o <- opts

  design <- prepare_zentangler_study_design(
    pheno_df = pheno_df, x_var = o$x_var, covariates = o$covariates,
    study_design = o$study_design, case_var = o$case_var, case_level = o$case_level,
    control_level = o$control_level, time_var = o$time_var, time_ref = o$time_ref,
    time_compare = o$time_compare, exposure_role = o$exposure_role,
    add_design_covariates = o$add_design_covariates,
    effect_x0 = o$effect_x0, effect_x1 = o$effect_x1
  )
  x_var <- design$x_var
  covariates <- design$covariates

  needed <- unique(c(x_var, o$y_var, o$survival_time_var, o$survival_event_var, covariates, o$bootstrap_id, o$a_stage$random_effects))
  sp <- zentangler_split_raw_blocks(blocks)
  al <- align_samples_multiview_blocks(sp$blocks, design$pheno_df, needed_pheno_cols = needed)
  B <- al$blocks
  B_raw <- zentangler_align_raw_blocks(sp$raw, B)
  views <- names(B)
  X <- al$pheno[[x_var]]
  Y <- zentangler_make_outcome(al$pheno, o)
  if (!is.numeric(X)) stop("The exposure must be numeric.", call. = FALSE)
  C <- make_covariate_matrix(al$pheno, covariates)
  verbose <- isTRUE(o$verbose)
  zentangler_log(verbose, sprintf(
    "Data: %d of %d samples complete across views; outcome %s; exposure %s%s.",
    length(al$sample_ids), nrow(pheno_df), o$y_family,
    if (identical(o$study_design, "standard")) x_var else paste0(o$study_design, " indicator"),
    if (length(covariates)) paste0("; covariates: ", paste(covariates, collapse = ", ")) else ""
  ))
  for (v in views) zentangler_log(verbose, sprintf("%s: %d features after filtering", v, ncol(B[[v]])), indent = 1L)

  # Total-effect check (reported, not used as a gate).
  contrast <- infer_effect_contrast(X, x0 = design$effect_x0, x1 = design$effect_x1)
  exposure_outcome <- if (isTRUE(skip_effects)) NULL else zentangler_exposure_outcome_test(Y, X, C, o$y_family, contrast[["x0"]], contrast[["x1"]])
  zentangler_log_exposure_outcome(verbose, exposure_outcome)

  # A-stage and screening.
  a_stage <- o$a_stage %||% zentangler_resolve_a_stage(o$a_stage_model, o$maaslin2_args, views, tweedieverse_args = o$tweedieverse_args)
  screens <- lapply(views, function(v) {
    zentangler_a_stage_screen(
      view = v, M = B[[v]], X = X, C = C, pheno_df = al$pheno, x_var = x_var, covariates = covariates,
      model = a_stage$model[[v]], maaslin2_args = a_stage$maaslin2_args[[v]],
      tweedieverse_args = a_stage$tweedieverse_args[[v]], M_raw = B_raw[[v]],
      sis_n = o$sis_n, rank_by = o$sis_rank, screen_method = o$screen_method, fdr_method = o$fdr_method
    )
  })
  names(screens) <- views
  selected <- lapply(views, function(v) B[[v]][, screens[[v]]$selected, drop = FALSE])
  names(selected) <- views
  zentangler_log(verbose, "A-stage (exposure -> mediator) and screening:")
  for (v in views) {
    st <- screens[[v]]$table
    zentangler_log(verbose, sprintf(
      "%s [%s]: %d features, %d with p_a < 0.05, %d screened in",
      v, a_stage$model[[v]], nrow(st), sum(st$p_a < 0.05, na.rm = TRUE), length(screens[[v]]$selected)
    ), indent = 1L)
  }
  zentangler_log(verbose, "Outcome model: ", if (identical(o$engine, "multiview")) "multiview cooperative learning (tuning rho by CV)" else paste0("IntegratedLearner-style ", sub("il_", "", o$engine), " fusion"), " ...")

  # Outcome model and B-path inference.
  yfit <- fit_y_multiview_stage(
    Y = Y, X = X, blocks = selected, C = C, engine = o$engine, y_family = o$y_family,
    lambda_choice = o$lambda_choice, glmnet_alpha = o$glmnet_alpha,
    rho = o$rho, rho_grid = o$rho_grid
  )
  if (verbose) {
    if (!is.null(yfit$fit$cv)) {
      cvt <- yfit$fit$cv
      zentangler_log(TRUE, sprintf("rho chosen = %g (CV %s by rho: %s)", yfit$rho, cvt$type_measure[1] %||% "score", paste(sprintf("%g:%.4g", cvt$rho, cvt$cv_score %||% cvt$cv_min), collapse = ", ")), indent = 1L)
    }
    zentangler_log(TRUE, sprintf("%d of %d screened mediators have non-zero outcome-model coefficients",
                                 sum(vapply(yfit$coef_by_view, function(b) sum(b != 0), numeric(1))),
                                 sum(vapply(selected, ncol, integer(1)))), indent = 1L)
    zentangler_log(TRUE, "B-path inference: ", if (identical(o$b_inference, "debiased")) "de-biased lasso ..." else "full-pipeline bootstrap ...")
  }
  b_infer <- infer_p_b_multiview(
    Y = Y, X = X, blocks = selected, C = C, y_family = o$y_family, method = o$b_inference,
    lambda_choice = o$lambda_choice, max_debias_targets = o$debias_max_targets,
    beta_init = zentangler_fusion_init(yfit, selected)
  )

  tab <- do.call(rbind, lapply(views, function(v) {
    meds <- colnames(B[[v]])
    b_fusion <- setNames(rep(0, length(meds)), meds)
    hit <- intersect(names(yfit$coef_by_view[[v]]), meds)
    b_fusion[hit] <- yfit$coef_by_view[[v]][hit]
    na <- setNames(rep(NA_real_, length(meds)), meds)
    if (identical(o$b_inference, "debiased")) {
      b <- se <- p <- na
      hit <- intersect(names(b_infer$b[[v]]), meds)
      b[hit] <- b_infer$b[[v]][hit]
      se[hit] <- b_infer$se[[v]][hit]
      p[hit] <- b_infer$p_b[[v]][hit]
    } else {
      b <- b_fusion
      se <- p <- na
    }
    data.frame(omics = v, assemble_mediator_table(screens[[v]]$table, b, b_fusion, se, p, o$fdr_method), stringsAsFactors = FALSE)
  }))
  if (identical(o$fdr_joint, "hdmt") && !identical(o$b_inference, "none")) {
    pm <- zentangler_marginal_b_pvalues(Y, X, B, C, o$y_family)
    tab$p_b_marginal <- NA_real_
    for (v in views) {
      idx <- tab$omics == v
      tab$p_b_marginal[idx] <- pm[[v]][tab$mediator[idx]]
    }
  }

  # Pipeline bootstrap (b / score inference).
  bootstrap <- list()
  if (o$bootstrap_repeats > 0L && (identical(o$b_inference, "bootstrap") || identical(o$primary_inference, "bootstrap_score"))) {
    rep_opts <- o
    rep_opts[c("x_var", "covariates", "study_design")] <- list(x_var, covariates, "standard")
    rep_opts$rho <- yfit$rho %||% o$rho  # rho fixed at its full-data CV value
    rep_opts$bootstrap_id <- NULL
    rep_opts$verbose <- FALSE
    zentangler_log(verbose, sprintf("Pipeline bootstrap: %d replicates (rho fixed at %s) ...", o$bootstrap_repeats, format(yfit$rho %||% NA)))
    boot <- bootstrap_multiview_fit(
      verbose = verbose,
      blocks = zentangler_attach_raw_blocks(B, B_raw, o$a_stage %||% a_stage), pheno_df = al$pheno, opts = rep_opts,
      key_ref = paste(tab$omics, tab$mediator, sep = "::"),
      bootstrap_repeats = o$bootstrap_repeats, bootstrap_seed = o$bootstrap_seed, bootstrap_id = o$bootstrap_id
    )
    tab <- add_bootstrap_columns(tab, boot$score_matrix, boot$active_matrix, ci_level = o$bootstrap_ci_level)
    tab <- add_b_bootstrap_columns(tab, boot$b_matrix, ci_level = o$bootstrap_ci_level)
    if (identical(o$b_inference, "bootstrap")) {
      tab$p_b <- tab$p_b_bootstrap
      tab <- zentangler_finalize_joint_p(tab)
    }
    if (identical(o$primary_inference, "bootstrap_score")) {
      tab$p_primary <- tab$p_score_boot
      tab$p_primary[!is.finite(tab$p_primary)] <- 1
    }
    bootstrap$pipeline <- boot
  }
  tab <- apply_primary_fdr_scope(
    tab, fdr_method = o$fdr_method, fdr_scope = o$fdr_scope, fdr_joint = o$fdr_joint,
    joint_ok = !identical(o$primary_inference, "bootstrap_score")
  )
  fdr_joint_used <- attr(tab, "fdr_joint_used")
  if (verbose) {
    zentangler_log(TRUE, sprintf("%d mediators tested in the B-stage (%s); FDR over %d features (%s, %s scope).",
                                 sum(tab$tested_b_stage), b_infer$method, nrow(tab), fdr_joint_used, o$fdr_scope))
    act <- vapply(c(0.05, 0.1, 0.25), function(q) sum(zentangler_is_active(tab, q)), integer(1))
    zentangler_log(TRUE, sprintf("Active mediators: %d at q <= 0.05, %d at q <= 0.10, %d at q <= 0.25", act[1], act[2], act[3]), indent = 1L)
  }

  # Post-selection causal effects.
  effects <- data.frame()
  if (!isTRUE(skip_effects)) {
    effect_set <- zentangler_effect_set(tab, q_threshold = o$effect_q_threshold)
    eff <- compute_effect_decomposition_multiview(effect_set, Y, X, B, C, o$y_family, contrast[["x0"]], contrast[["x1"]])
    effects <- eff$effects
    tab <- add_mediator_effect_columns(tab, effects, eff$mediator_effects)
    zentangler_log(verbose, sprintf(
      "Causal effects (%d mediators with q <= %g; %s, X %.3g -> %.3g): NDE %.4g, NIE %.4g, TE %.4g, proportion mediated %.3g",
      nrow(effect_set), o$effect_q_threshold, gsub("_", " ", effects$effect_scale[[1L]]),
      contrast[["x0"]], contrast[["x1"]], effects$nde[[1L]], effects$nie_total[[1L]], effects$te[[1L]], effects$prop_mediated[[1L]]
    ))
    if (identical(o$causal_inference, "bootstrap")) {
      zentangler_log(verbose, sprintf("Causal-effect bootstrap: %d replicates (fixed mediator set) ...", o$bootstrap_repeats))
      eb <- bootstrap_effects_fixed_set(
        verbose = verbose,
        blocks = B, pheno_df = al$pheno, X = X, Y = Y, C = C, effect_set = effect_set,
        y_family = o$y_family, x0 = contrast[["x0"]], x1 = contrast[["x1"]],
        bootstrap_repeats = o$bootstrap_repeats, bootstrap_seed = o$bootstrap_seed, bootstrap_id = o$bootstrap_id
      )
      effects <- add_effect_bootstrap_columns(effects, eb$effect_decomposition, ci_level = o$bootstrap_ci_level)
      tab <- add_mediator_effect_bootstrap_columns(tab, eb$indirect_matrix, ci_level = o$bootstrap_ci_level)
      bootstrap$effects <- eb
    }
  }
  tab <- tab[order(tab$abs_score, decreasing = TRUE), , drop = FALSE]
  rownames(tab) <- NULL

  end_time <- Sys.time()
  settings <- o
  settings[c("x_var", "covariates")] <- list(x_var, covariates)
  settings <- utils::modifyList(settings, list(
    design_info = design$design_info,
    view_names = views,
    effect_x0 = contrast[["x0"]],
    effect_x1 = contrast[["x1"]],
    rho_used = yfit$rho %||% NA_real_,
    rho_cv = yfit$fit$cv,
    y_model_method = yfit$method,
    b_inference_method = b_infer$method,
    fdr_joint_used = fdr_joint_used
  ))
  out <- list(
    settings = settings,
    diagnostics = list(
      runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
      n_samples_input = nrow(pheno_df),
      n_samples = length(al$sample_ids),
      sample_ids = al$sample_ids,
      features_before_filtering = vapply(sp$blocks, function(x) ncol(as.data.frame(x)), integer(1)),
      features_after_filtering = vapply(B, ncol, integer(1)),
      screened_counts = vapply(screens, function(x) length(x$selected), integer(1)),
      y_model_nonzero = vapply(views, function(v) sum(tab$b_fusion[tab$omics == v] != 0), integer(1)),
      n_b_stage_tested = sum(tab$tested_b_stage),
      n_fdr_tests = nrow(tab),
      exposure_outcome = exposure_outcome,
      bootstrap_successful = c(pipeline = bootstrap$pipeline$repeats_successful %||% NA_integer_, effects = bootstrap$effects$repeats_successful %||% NA_integer_),
      bootstrap_failures = c(bootstrap$pipeline$failures, bootstrap$effects$failures)
    ),
    mediators = tab,
    views = lapply(screens, function(s) {
      v <- list(a_stage_model = s$a_stage_model %||% "lm", screen = s$table, selected = s$selected)
      if (!is.null(s$tweedie_index)) v$tweedie_index <- s$tweedie_index
      v
    }),
    effects = effects,
    bootstrap = if (length(bootstrap) > 0L) bootstrap else NULL,
    model = if (isTRUE(o$return_fits)) yfit else NULL
  )
  class(out) <- c("zentangler_fit", "list")
  zentangler_log(verbose, sprintf("Done in %.1f s.", out$diagnostics$runtime_seconds))
  out
}

zentangler_make_outcome <- function(pheno, o) {
  if (identical(o$y_family, "survival")) {
    time <- suppressWarnings(as.numeric(pheno[[o$survival_time_var]]))
    event <- suppressWarnings(as.numeric(pheno[[o$survival_event_var]]))
    if (any(!is.finite(time)) || any(time <= 0)) stop("survival_time_var must contain finite positive times.", call. = FALSE)
    if (any(!is.finite(event)) || !all(event %in% c(0, 1))) stop("survival_event_var must be coded 0/1.", call. = FALSE)
    return(survival::Surv(time = time, event = event))
  }
  Y <- pheno[[o$y_var]]
  if (!is.numeric(Y)) stop("y_var must be numeric.", call. = FALSE)
  if (identical(o$y_family, "binomial")) zentangler_validate_binary_outcome(Y)
  Y
}
