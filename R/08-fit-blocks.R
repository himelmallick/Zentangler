# -----------------------------------------------------------------------------
# Public K-view fit
# -----------------------------------------------------------------------------

fit_multiview_parallel_zentangler_blocks <- function(
  blocks,
  pheno_df,
  x_var = NULL,
  y_var = NULL,
  study_design = c("standard", "case_control", "time", "case_control_time"),
  case_var = NULL,
  case_level = NULL,
  control_level = NULL,
  time_var = NULL,
  time_ref = NULL,
  time_compare = NULL,
  exposure_role = c("time", "case", "interaction"),
  add_design_covariates = TRUE,
  method_preset = c("custom", "fast_lasso", "elastic_net", "longitudinal_maaslin2", "full_exploratory"),
  covariates = NULL,
  residualize = FALSE,
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  a_stage_model = c("lm", "maaslin2"),
  maaslin2_random_effect = NULL,
  maaslin2_normalization = "NONE",
  maaslin2_transform = "NONE",
  maaslin2_analysis_method = "LM",
  maaslin2_standardize = FALSE,
  maaslin2_output_dir = NULL,
  fusion_mode = c("early", "intermediate", "late"),
  y_family = c("gaussian", "binomial", "survival"),
  survival_time_var = NULL,
  survival_event_var = NULL,
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global", "within_view"),
  primary_inference = c("model_based", "bootstrap_score"),
  b_inference = c("debiased", "bootstrap"),
  causal_inference = c("none", "bootstrap"),
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
  effect_x0 = NULL,
  effect_x1 = NULL
) {
  start_time <- Sys.time()
  method_preset <- match.arg(method_preset)
  apply_zentangler_preset(method_preset, environment())
  sis_rank <- match.arg(sis_rank)
  screen_method <- match.arg(screen_method)
  a_stage_model <- match.arg(a_stage_model)
  fusion_mode <- match.arg(fusion_mode)
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)
  fdr_method <- validate_fdr_method(fdr_method)
  fdr_scope <- validate_fdr_scope(fdr_scope)
  primary_inference <- validate_primary_inference(primary_inference)
  b_inference <- normalize_b_inference_method(b_inference)
  causal_inference <- match.arg(causal_inference)
  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  set.seed(seed)
  if (!identical(y_family, "survival") && (is.null(y_var) || length(y_var) != 1L || !nzchar(as.character(y_var)))) {
    stop("y_var is required.", call. = FALSE)
  }
  y_var <- if (is.null(y_var)) NULL else as.character(y_var)
  if (identical(y_family, "survival")) {
    if (is.null(survival_time_var) || length(survival_time_var) != 1L || !nzchar(as.character(survival_time_var))) {
      stop("survival_time_var is required when y_family='survival'.", call. = FALSE)
    }
    if (is.null(survival_event_var) || length(survival_event_var) != 1L || !nzchar(as.character(survival_event_var))) {
      stop("survival_event_var is required when y_family='survival'.", call. = FALSE)
    }
    survival_time_var <- as.character(survival_time_var)
    survival_event_var <- as.character(survival_event_var)
    if (identical(fusion_mode, "intermediate")) {
      stop("Survival outcomes currently support fusion_mode='early' or 'late'.", call. = FALSE)
    }
  }

  design <- prepare_zentangler_study_design(
    pheno_df = pheno_df,
    x_var = x_var,
    covariates = covariates,
    study_design = study_design,
    case_var = case_var,
    case_level = case_level,
    control_level = control_level,
    time_var = time_var,
    time_ref = time_ref,
    time_compare = time_compare,
    exposure_role = exposure_role,
    add_design_covariates = add_design_covariates,
    effect_x0 = effect_x0,
    effect_x1 = effect_x1
  )
  pheno_df <- design$pheno_df
  x_var <- design$x_var
  covariates <- design$covariates
  effect_x0 <- design$effect_x0
  effect_x1 <- design$effect_x1

  feature_counts_before <- vapply(blocks, function(x) ncol(as.data.frame(x)), integer(1))

  needed <- unique(c(x_var, y_var, survival_time_var, survival_event_var, covariates, bootstrap_id, maaslin2_random_effect))
  al <- align_samples_multiview_blocks(blocks, pheno_df, needed_pheno_cols = needed)
  blocks0 <- al$blocks
  feature_counts_after <- vapply(blocks0, ncol, integer(1))
  X_raw <- al$pheno[[x_var]]
  if (identical(y_family, "survival")) {
    time_raw <- suppressWarnings(as.numeric(al$pheno[[survival_time_var]]))
    event_raw <- suppressWarnings(as.numeric(al$pheno[[survival_event_var]]))
    if (any(!is.finite(time_raw)) || any(time_raw <= 0)) {
      stop("survival_time_var must contain finite positive survival times after sample alignment.", call. = FALSE)
    }
    if (any(!is.finite(event_raw)) || !all(event_raw %in% c(0, 1))) {
      stop("survival_event_var must contain 0/1 event indicators after sample alignment.", call. = FALSE)
    }
    Y_raw <- survival::Surv(time = time_raw, event = event_raw)
  } else {
    Y_raw <- al$pheno[[y_var]]
  }
  if (!is.numeric(X_raw) || (!identical(y_family, "survival") && !is.numeric(Y_raw))) {
    stop("x_var and y_var must be numeric.", call. = FALSE)
  }

  C <- make_covariate_matrix(al$pheno, covariates)
  blocks_model <- blocks0
  X_model <- X_raw
  Y_model <- Y_raw

  if (residualize && !is.null(C) && ncol(C) > 0) {
    blocks_model <- lapply(blocks0, residualize_matrix, C = C)
    X_model <- as.numeric(residualize_matrix(matrix(X_raw, ncol = 1), C))
    if (identical(y_family, "survival")) {
      warning("residualize=TRUE does not residualize survival outcomes; covariates are modeled directly in the Cox Y-stage.")
      Y_model <- Y_raw
      C_model <- C
    } else {
      Y_model <- as.numeric(residualize_matrix(matrix(Y_raw, ncol = 1), C))
      C_model <- NULL
    }
  } else {
    C_model <- C
  }

  screens <- lapply(names(blocks_model), function(view) {
    if (identical(a_stage_model, "maaslin2")) {
      if (residualize) {
        warning("MaAsLin2 A-stage uses non-residualized X/mediators and models covariates/random effects directly.")
      }
      view_out_dir <- if (!is.null(maaslin2_output_dir) && nzchar(maaslin2_output_dir)) {
        file.path(maaslin2_output_dir, paste0("a_stage_", view))
      } else {
        NULL
      }
      hima_screen_mediators_maaslin2(
        M = blocks0[[view]],
        pheno_df = al$pheno,
        x_var = x_var,
        covariates = covariates,
        random_effects = maaslin2_random_effect,
        sis_n = sis_n,
        rank_by = sis_rank,
        screen_method = screen_method,
        fdr_method = fdr_method,
        normalization = maaslin2_normalization,
        transform = maaslin2_transform,
        analysis_method = maaslin2_analysis_method,
        standardize = maaslin2_standardize,
        output_dir = view_out_dir
      )
    } else {
      hima_screen_mediators(
        blocks_model[[view]],
        X_model,
        C = C_model,
        sis_n = sis_n,
        rank_by = sis_rank,
        screen_method = screen_method,
        fdr_method = fdr_method
      )
    }
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
    glmnet_alpha = glmnet_alpha,
    coop_rho = coop_rho,
    coop_maxit = coop_maxit,
    coop_tol = coop_tol
  )

  b_inference_requested <- b_inference
  if (identical(b_inference_requested, "bootstrap") && bootstrap_repeats < 1L) {
    warning("b_inference='bootstrap' requires bootstrap_repeats > 0; returning initial B-path estimates without bootstrap p-values.")
  }
  if (identical(causal_inference, "bootstrap") && bootstrap_repeats < 1L) {
    warning("causal_inference='bootstrap' requires bootstrap_repeats > 0; returning point estimates without bootstrap intervals.")
  }

  p_b_infer <- infer_p_b_multiview(
    Y = Y_model,
    X = X_model,
    blocks = selected_blocks,
    C = C_model,
    coef_by_view = yfit$coef_by_view,
    y_family = y_family,
    method = if (identical(b_inference_requested, "bootstrap")) "bootstrap" else "debiased",
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

    assemble_mediator_table(screens[[view]]$table, b_vec = b_full, p_b_vec = p_full, fdr_method = fdr_method)
  })
  names(view_tables) <- names(blocks_model)

  combined <- do.call(
    rbind,
    lapply(names(view_tables), function(view) {
      data.frame(omics = view, view_tables[[view]], stringsAsFactors = FALSE)
    })
  )
  combined$p_primary_model_based <- combined$joint_p_ab
  combined$p_primary <- combined$p_primary_model_based
  combined <- apply_primary_fdr_scope(combined, fdr_method = fdr_method, fdr_scope = fdr_scope)
  combined <- combined[order(combined$abs_score, decreasing = TRUE), , drop = FALSE]
  rownames(combined) <- NULL

  effect_decomposition <- compute_effect_decomposition_multiview(
    combined_mediators = combined,
    direct_effect = yfit$x_coef,
    Y = Y_model,
    X = X_model,
    blocks = blocks_model,
    C = C_model,
    y_family = y_family,
    x0 = effect_x0,
    x1 = effect_x1
  )
  combined <- add_mediator_effect_columns(combined, effect_decomposition)

  key_ref <- paste(combined$omics, combined$mediator, sep = "::")
  bootstrap <- NULL
  use_bootstrap <- bootstrap_repeats > 0L && (
    identical(b_inference_requested, "bootstrap") ||
      identical(causal_inference, "bootstrap") ||
      identical(primary_inference, "bootstrap_score")
  )
  if (use_bootstrap) {
    bootstrap <- bootstrap_multiview_fit(
      blocks = blocks0,
      pheno_df = al$pheno,
      x_var = x_var,
      y_var = y_var,
      survival_time_var = survival_time_var,
      survival_event_var = survival_event_var,
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
      maaslin2_output_dir = maaslin2_output_dir,
      fusion_mode = fusion_mode,
      y_family = y_family,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha,
      fdr_method = fdr_method,
      fdr_scope = fdr_scope,
      primary_inference = primary_inference,
      b_inference = b_inference,
      debias_max_targets = debias_max_targets,
      coop_rho = coop_rho,
      coop_maxit = coop_maxit,
      coop_tol = coop_tol,
      key_ref = key_ref,
      bootstrap_repeats = bootstrap_repeats,
      bootstrap_ci_level = bootstrap_ci_level,
      bootstrap_seed = bootstrap_seed,
      bootstrap_id = bootstrap_id,
      effect_x0 = effect_x0,
      effect_x1 = effect_x1
    )
    combined <- add_bootstrap_columns(
      combined = combined,
      score_mat = bootstrap$score_matrix,
      active_mat = bootstrap$active_matrix,
      ci_level = bootstrap_ci_level
    )
    combined <- add_b_bootstrap_columns(
      combined = combined,
      b_mat = bootstrap$b_matrix,
      ci_level = bootstrap_ci_level
    )
    combined <- add_mediator_effect_bootstrap_columns(
      combined = combined,
      score_mat = bootstrap$score_matrix,
      effect_decomposition = effect_decomposition,
      ci_level = bootstrap_ci_level
    )
    if (identical(b_inference_requested, "bootstrap") && "p_b_bootstrap" %in% colnames(combined)) {
      has_boot_p <- is.finite(combined$p_b_bootstrap)
      combined$p_b_initial <- combined$p_b
      combined$p_b_model_based <- combined$p_b
      combined$p_b[has_boot_p] <- combined$p_b_bootstrap[has_boot_p]
      combined$joint_p_ab <- NA_real_
      both_legs <- is.finite(combined$p_a) & is.finite(combined$p_b)
      combined$joint_p_ab[both_legs] <- pmax(combined$p_a[both_legs], combined$p_b[both_legs])
      combined$p_primary_model_based <- combined$joint_p_ab
      combined$p_primary <- combined$p_primary_model_based
      combined <- apply_primary_fdr_scope(combined, fdr_method = fdr_method, fdr_scope = fdr_scope)
      combined <- combined[order(combined$abs_score, decreasing = TRUE), , drop = FALSE]
      rownames(combined) <- NULL
      p_b_infer$method <- "bootstrap_b_path"
    }
    if (identical(causal_inference, "bootstrap")) {
      effect_decomposition <- add_effect_bootstrap_columns(
        effect_decomposition = effect_decomposition,
        effect_boot = bootstrap$effect_decomposition,
        ci_level = bootstrap_ci_level
      )
    }
  }
  causal_stage <- build_causal_stage_summary(
    effect_decomposition = effect_decomposition,
    combined_mediators = combined,
    causal_inference = causal_inference,
    bootstrap = bootstrap
  )

  primary_inference_used <- "model_based"
  if (identical(primary_inference, "bootstrap_score")) {
    if (is.null(bootstrap) || !("p_score_boot" %in% colnames(combined))) {
      warning("primary_inference='bootstrap_score' requires bootstrap_repeats > 0; using model-based p_primary.")
    } else {
      combined$p_primary <- combined$p_score_boot
      primary_inference_used <- "bootstrap_score"
      combined <- apply_primary_fdr_scope(combined, fdr_method = fdr_method, fdr_scope = fdr_scope)
      combined <- combined[order(combined$abs_score, decreasing = TRUE), , drop = FALSE]
      rownames(combined) <- NULL
    }
  }

  for (view in names(view_tables)) {
    idx <- combined$omics == view
    p_b_map <- setNames(combined$p_b[idx], combined$mediator[idx])
    joint_map <- setNames(combined$joint_p_ab[idx], combined$mediator[idx])
    hit_p <- match(view_tables[[view]]$mediator, names(p_b_map))
    ok_p <- !is.na(hit_p)
    view_tables[[view]]$p_b[ok_p] <- p_b_map[hit_p[ok_p]]
    hit_joint <- match(view_tables[[view]]$mediator, names(joint_map))
    ok_joint <- !is.na(hit_joint)
    view_tables[[view]]$joint_p_ab[ok_joint] <- joint_map[hit_joint[ok_joint]]
    view_tables[[view]]$p_primary <- view_tables[[view]]$joint_p_ab
  }

  views_out <- lapply(names(blocks_model), function(view) {
    list(
      screened_mediators = screens[[view]]$table,
      selected_mediators = screens[[view]]$selected,
      a = setNames(screens[[view]]$table$a_hat, screens[[view]]$table$mediator),
      p_a = setNames(screens[[view]]$table$p_a, screens[[view]]$table$mediator),
      b = setNames(view_tables[[view]]$b, view_tables[[view]]$mediator),
      p_b = setNames(view_tables[[view]]$p_b, view_tables[[view]]$mediator),
      mediators = view_tables[[view]],
      a_stage_model = if (!is.null(screens[[view]]$a_stage_model)) screens[[view]]$a_stage_model else a_stage_model,
      maaslin2_output_dir = screens[[view]]$output_dir,
      maaslin2_fixed_effects = screens[[view]]$fixed_effects,
      maaslin2_random_effects = screens[[view]]$random_effects
    )
  })
  names(views_out) <- names(blocks_model)

  screened_counts <- vapply(screens, function(x) length(x$selected), integer(1))
  active_counts <- vapply(view_tables, function(x) sum(is.finite(x$b) & x$b != 0, na.rm = TRUE), integer(1))
  end_time <- Sys.time()
  diagnostics <- list(
    package_version = tryCatch(as.character(utils::packageVersion("Zentangler")), error = function(e) NA_character_),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    started = start_time,
    finished = end_time,
    runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
    n_samples_input = nrow(pheno_df),
    n_samples = length(al$sample_ids),
    n_views = length(blocks_model),
    features_before_filtering = feature_counts_before,
    features_after_filtering = feature_counts_after,
    screened_counts = screened_counts,
    active_b_counts = active_counts,
    bootstrap_failures = if (!is.null(bootstrap)) bootstrap$failures else character(0)
  )
  settings_out <- list(
    method_preset = method_preset,
    input_container = "matrix_blocks",
    study_design = design$design_info$study_design,
    design_info = design$design_info,
    x_var = x_var,
    y_var = y_var,
    survival_time_var = survival_time_var,
    survival_event_var = survival_event_var,
    view_names = names(blocks_model),
    n_views = length(blocks_model),
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
    maaslin2_output_dir = maaslin2_output_dir,
    fusion_mode = fusion_mode,
    y_family = y_family,
    survival_time_var = survival_time_var,
    survival_event_var = survival_event_var,
    lambda_choice = lambda_choice,
    glmnet_alpha = glmnet_alpha,
    fdr_method = fdr_method,
    fdr_scope = fdr_scope,
    primary_inference = primary_inference,
    primary_inference_used = primary_inference_used,
    glmnet_penalty = glmnet_alpha_label(if (identical(fusion_mode, "intermediate")) 1 else glmnet_alpha),
    glmnet_alpha_used_in_y_stage = if (identical(fusion_mode, "intermediate")) 1 else glmnet_alpha,
    b_inference = b_inference_requested,
    b_inference_method = p_b_infer$method,
    b_inference_note = if (identical(b_inference_requested, "bootstrap") && bootstrap_repeats < 1L) {
      "Requested bootstrap B-path inference, but bootstrap_repeats was 0; initial B-path p-values were returned without bootstrap replacement."
    } else if (identical(b_inference_requested, "bootstrap")) {
      "B-path p-values use bootstrap refits of the B-stage; p_b_initial retains the initial single-fit p-value."
    } else if (identical(b_inference_requested, "debiased") && !isTRUE(all.equal(glmnet_alpha, 1))) {
      "B-path p-values use the current HIMA-style de-biased lasso routine; elastic-net alpha changes early/late Y-stage coefficient fitting, not the de-biased inference formula."
    } else {
      NA_character_
    },
    causal_inference = causal_inference,
    debias_max_targets = debias_max_targets,
    y_model_method = yfit$method,
    bootstrap_repeats = as.integer(bootstrap_repeats),
    bootstrap_ci_level = bootstrap_ci_level,
    bootstrap_seed = bootstrap_seed,
    bootstrap_id = bootstrap_id,
    effect_x0 = effect_decomposition$x0,
    effect_x1 = effect_decomposition$x1,
    effect_method = effect_decomposition$effect_method,
    effect_scale = effect_decomposition$effect_scale,
    effects_summary = as.list(effect_decomposition[1, , drop = FALSE]),
    causal_stage = causal_stage$stage,
    causal_stage_method = causal_stage$stage_method
  )
  mediators_active_default <- zentangler_active_mediators(list(combined_mediators = combined), q_threshold = 0.25)
  view_summary_default <- zentangler_compute_view_summary(combined, q_threshold = 0.25)
  model_summary_default <- zentangler_compute_model_summary(
    settings = settings_out,
    diagnostics = diagnostics,
    tab = combined,
    q_threshold = 0.25
  )

  list(
    settings = settings_out,
    sample_ids = rownames(al$pheno),
    diagnostics = diagnostics,
    views = views_out,
    combined_mediators = combined,
    mediators_all = combined,
    mediators_active = mediators_active_default,
    mediators_top = zentangler_top_mediators(list(combined_mediators = combined), n = 20L),
    view_summary = view_summary_default,
    model_summary = model_summary_default,
    x_to_y_coef = yfit$x_coef,
    causal_stage = causal_stage,
    effects = effect_decomposition,
    effect_decomposition = effect_decomposition,
    causal_effects = effect_decomposition,
    bootstrap = bootstrap,
    fits = if (return_fits) yfit else NULL
  )
}
