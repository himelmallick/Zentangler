fit_multiview_parallel_zentangler <- function(
  mae,
  x_var = NULL,
  y_var = NULL,
  view_names = NULL,
  assay_names = NULL,
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
  effect_x1 = NULL,
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
  method_preset <- match.arg(method_preset)
  fdr_method <- validate_fdr_method(fdr_method)
  fdr_scope <- validate_fdr_scope(fdr_scope)
  primary_inference <- validate_primary_inference(primary_inference)
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
    study_design = study_design,
    case_var = case_var,
    case_level = case_level,
    control_level = control_level,
    time_var = time_var,
    time_ref = time_ref,
    time_compare = time_compare,
    exposure_role = exposure_role,
    add_design_covariates = add_design_covariates,
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
    survival_time_var = survival_time_var,
    survival_event_var = survival_event_var,
    lambda_choice = lambda_choice,
    glmnet_alpha = glmnet_alpha,
    fdr_method = fdr_method,
    fdr_scope = fdr_scope,
    primary_inference = primary_inference,
    b_inference = b_inference,
    causal_inference = causal_inference,
    debias_max_targets = debias_max_targets,
    coop_rho = coop_rho,
    coop_maxit = coop_maxit,
    coop_tol = coop_tol,
    seed = seed,
    return_fits = return_fits,
    bootstrap_repeats = bootstrap_repeats,
    bootstrap_ci_level = bootstrap_ci_level,
    bootstrap_seed = bootstrap_seed,
    bootstrap_id = bootstrap_id,
    effect_x0 = effect_x0,
    effect_x1 = effect_x1
  )

  fit$settings$input_container <- "MultiAssayExperiment"
  fit$settings$mae_view_map <- inputs$view_map
  fit$model_summary <- zentangler_model_summary(fit)
  fit
}
