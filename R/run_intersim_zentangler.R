#' Run InterSIM/SIMMBA Zentangler Benchmarks
#'
#' Convenience wrapper for generating InterSIM/SIMMBA MAE simulations and
#' benchmarking Zentangler fusion modes against known mediator truth.
#'
#' @param nrep Number of simulation replicates.
#' @param nsample Number of samples per simulated dataset before train/test split.
#' @param fusion_modes Fusion modes to evaluate.
#' @param sis_n Number of mediators retained per view when \code{screen_method = "sis"}.
#' @param sis_rank Ranking criterion for SIS.
#' @param screen_method Screening method passed to \code{fit_multiview_parallel_zentangler()}.
#' @param a_stage_model A-stage model: \code{"lm"} for HIMA-like univariate linear models or \code{"maaslin2"} for MaAsLin2 fixed/random-effect screening.
#' @param q_threshold BH q-value threshold used to define active mediators.
#' @param top_n Number of top-ranked mediators used for top-k truth recovery.
#' @param seed Random seed passed to \code{gen_simmba()} and model fits.
#' @param p.train Training split proportion passed to \code{gen_simmba()}.
#' @param ygen.mode Outcome generation mode passed to \code{gen_simmba()}.
#' @param outcome.type Outcome type passed to \code{gen_simmba()}.
#' @param lambda_choice Cross-validated glmnet lambda choice.
#' @param glmnet_alpha Glmnet mixing parameter for early and late fusion. Use 1 for lasso, values between 0 and 1 for elastic net, and 0 for ridge. Intermediate fusion currently uses the cooperative lasso-style update.
#' @param maaslin2_random_effect Optional MaAsLin2 random-effect column name.
#' @param maaslin2_normalization MaAsLin2 normalization method.
#' @param maaslin2_transform MaAsLin2 transform method.
#' @param maaslin2_analysis_method MaAsLin2 analysis method.
#' @param maaslin2_standardize Logical passed to MaAsLin2.
#' @param maaslin2_output_dir Optional directory for MaAsLin2 A-stage outputs.
#' @param b_inference B-path inference method.
#' @param debias_max_targets Maximum active targets for de-biased lasso inference.
#' @param coop_rho Agreement penalty strength for intermediate cooperative fusion.
#' @param residualize Logical. Residualize X, Y, and mediators against covariates.
#' @param covariates Optional covariate names in \code{colData(mae)}.
#' @param bootstrap_repeats Number of bootstrap repeats inside each fit.
#' @param bootstrap_id Optional bootstrap cluster ID column.
#' @param verbose Logical. Print progress.
#' @param ... Additional arguments passed to \code{gen_simmba()}.
#'
#' @return A list with \code{summary}, \code{detail}, \code{fits}, \code{truth}, \code{errors}, and \code{settings}.
#' @export
run_intersim_zentangler <- function(
  nrep = 10,
  nsample = 600,
  fusion_modes = c("early", "intermediate", "late"),
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  a_stage_model = c("lm", "maaslin2"),
  q_threshold = 0.25,
  top_n = 50L,
  seed = 1234,
  p.train = 0.7,
  ygen.mode = c("LM", "Friedman", "Friedman2"),
  outcome.type = c("continuous", "binary"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  maaslin2_random_effect = NULL,
  maaslin2_normalization = "NONE",
  maaslin2_transform = "NONE",
  maaslin2_analysis_method = "LM",
  maaslin2_standardize = FALSE,
  maaslin2_output_dir = NULL,
  b_inference = c("debiased_lasso", "refit"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  residualize = FALSE,
  covariates = NULL,
  bootstrap_repeats = 0L,
  bootstrap_id = NULL,
  verbose = TRUE,
  ...
) {
  fusion_modes <- match.arg(fusion_modes, choices = c("early", "intermediate", "late"), several.ok = TRUE)
  sis_rank <- match.arg(sis_rank)
  screen_method <- match.arg(screen_method)
  a_stage_model <- match.arg(a_stage_model)
  ygen.mode <- match.arg(ygen.mode)
  outcome.type <- match.arg(outcome.type)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)
  b_inference <- match.arg(b_inference)

  if (identical(outcome.type, "binary")) {
    y_family <- "binomial"
  } else {
    y_family <- "gaussian"
  }

  sim <- gen_simmba(
    nsample = nsample,
    nrep = nrep,
    seed = seed,
    p.train = p.train,
    ygen.mode = ygen.mode,
    outcome.type = outcome.type,
    ...
  )

  detail_rows <- list()
  fit_list <- list()
  errors <- character(0)

  for (rep_i in seq_len(nrep)) {
    mae_i <- sim$trainMae[[rep_i]]
    truth_i <- sim$truthDat[[rep_i]]

    if (is.null(mae_i)) {
      msg <- paste0("rep_", rep_i, ": missing trainMae")
      errors <- c(errors, msg)
      next
    }

    view_names <- names(MultiAssayExperiment::experiments(mae_i))
    truth_key <- zentangler_truth_key(truth_i)

    for (fusion_mode in fusion_modes) {
      if (isTRUE(verbose)) {
        message("Running rep ", rep_i, " / ", nrep, ", fusion = ", fusion_mode,
                ", screen_method = ", screen_method,
                ", a_stage_model = ", a_stage_model,
                ", glmnet_alpha = ", glmnet_alpha)
      }

      fit_i <- try(
        fit_multiview_parallel_zentangler(
          mae = mae_i,
          x_var = "A",
          y_var = "Y",
          view_names = view_names,
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
          b_inference = b_inference,
          debias_max_targets = debias_max_targets,
          coop_rho = coop_rho,
          seed = seed + rep_i,
          return_fits = FALSE,
          bootstrap_repeats = bootstrap_repeats,
          bootstrap_id = bootstrap_id
        ),
        silent = TRUE
      )

      key <- paste0("rep_", rep_i, "__", fusion_mode)
      if (inherits(fit_i, "try-error")) {
        errors <- c(errors, paste0(key, ": ", conditionMessage(attr(fit_i, "condition"))))
        next
      }

      fit_list[[key]] <- fit_i
      detail_rows[[length(detail_rows) + 1L]] <- zentangler_score_truth_recovery(
        tab = fit_i$combined_mediators,
        truth_key = truth_key,
        rep = rep_i,
        fusion_mode = fusion_mode,
        q_threshold = q_threshold,
        top_n = top_n
      )
    }
  }

  detail <- if (length(detail_rows) > 0) {
    out <- do.call(rbind, detail_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame()
  }

  summary <- zentangler_summarize_truth_recovery(detail)

  list(
    summary = summary,
    detail = detail,
    fits = fit_list,
    truth = sim$truthDat,
    simulation = sim,
    errors = errors,
    settings = list(
      nrep = nrep,
      nsample = nsample,
      fusion_modes = fusion_modes,
      sis_n = sis_n,
      sis_rank = sis_rank,
      screen_method = screen_method,
      a_stage_model = a_stage_model,
      q_threshold = q_threshold,
      top_n = top_n,
      seed = seed,
      p.train = p.train,
      ygen.mode = ygen.mode,
      outcome.type = outcome.type,
      y_family = y_family,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha,
      maaslin2_random_effect = maaslin2_random_effect,
      maaslin2_normalization = maaslin2_normalization,
      maaslin2_transform = maaslin2_transform,
      maaslin2_analysis_method = maaslin2_analysis_method,
      maaslin2_standardize = maaslin2_standardize,
      maaslin2_output_dir = maaslin2_output_dir,
      b_inference = b_inference,
      debias_max_targets = debias_max_targets,
      coop_rho = coop_rho,
      residualize = residualize,
      covariates = covariates,
      bootstrap_repeats = bootstrap_repeats,
      bootstrap_id = bootstrap_id
    )
  )
}

zentangler_truth_key <- function(truth) {
  if (is.null(truth) || nrow(truth) == 0) {
    return(data.frame(key = character(0), is_true = logical(0), stringsAsFactors = FALSE))
  }

  truth <- as.data.frame(truth, stringsAsFactors = FALSE)
  if (!all(c("featureID", "featureType") %in% colnames(truth))) {
    stop("truth table must contain featureID and featureType columns.")
  }

  raw_types <- unique(as.character(truth$featureType))
  raw_types <- raw_types[!is.na(raw_types) & nzchar(raw_types)]
  safe_type_map <- stats::setNames(simmba_sanitize_view_names(raw_types), raw_types)
  feature_type <- unname(safe_type_map[as.character(truth$featureType)])
  is_true <- if ("TIE_true" %in% colnames(truth)) {
    is.finite(truth$TIE_true) & truth$TIE_true != 0
  } else if ("is_mediator" %in% colnames(truth)) {
    as.logical(truth$is_mediator)
  } else {
    rep(FALSE, nrow(truth))
  }

  data.frame(
    key = paste(feature_type, as.character(truth$featureID), sep = "::"),
    is_true = is_true,
    stringsAsFactors = FALSE
  )
}

zentangler_score_truth_recovery <- function(tab,
                                            truth_key,
                                            rep,
                                            fusion_mode,
                                            q_threshold,
                                            top_n) {
  if (is.null(tab) || nrow(tab) == 0) {
    return(data.frame(
      rep = rep,
      fusion_mode = fusion_mode,
      n_active = 0,
      true_active = 0,
      false_active = 0,
      n_true = sum(truth_key$is_true, na.rm = TRUE),
      precision = NA_real_,
      recall = 0,
      fdr = NA_real_,
      top50_true = 0,
      top50_precision = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  key <- paste(as.character(tab$omics), as.character(tab$mediator), sep = "::")
  truth_map <- stats::setNames(truth_key$is_true, truth_key$key)
  is_true <- as.logical(truth_map[key])
  is_true[is.na(is_true)] <- FALSE

  active <- is.finite(tab$q_primary) & tab$q_primary <= q_threshold &
    is.finite(tab$b) & tab$b != 0

  n_active <- sum(active, na.rm = TRUE)
  true_active <- sum(active & is_true, na.rm = TRUE)
  false_active <- sum(active & !is_true, na.rm = TRUE)
  n_true <- sum(truth_key$is_true, na.rm = TRUE)

  ord <- order(tab$abs_score, decreasing = TRUE, na.last = NA)
  top_idx <- head(ord, min(as.integer(top_n), length(ord)))
  top_true <- sum(is_true[top_idx], na.rm = TRUE)
  top_precision <- if (length(top_idx) > 0) top_true / length(top_idx) else NA_real_

  data.frame(
    rep = rep,
    fusion_mode = fusion_mode,
    n_active = n_active,
    true_active = true_active,
    false_active = false_active,
    n_true = n_true,
    precision = if (n_active > 0) true_active / n_active else NA_real_,
    recall = if (n_true > 0) true_active / n_true else NA_real_,
    fdr = if (n_active > 0) false_active / n_active else NA_real_,
    top50_true = top_true,
    top50_precision = top_precision,
    stringsAsFactors = FALSE
  )
}

zentangler_summarize_truth_recovery <- function(detail) {
  if (is.null(detail) || nrow(detail) == 0) return(data.frame())

  fusion_levels <- unique(as.character(detail$fusion_mode))
  out <- lapply(fusion_levels, function(fm) {
    d <- detail[as.character(detail$fusion_mode) == fm, , drop = FALSE]
    data.frame(
      fusion_mode = fm,
      n_active = mean(d$n_active, na.rm = TRUE),
      true_active = mean(d$true_active, na.rm = TRUE),
      false_active = mean(d$false_active, na.rm = TRUE),
      precision = mean(d$precision, na.rm = TRUE),
      recall = mean(d$recall, na.rm = TRUE),
      fdr = mean(d$fdr, na.rm = TRUE),
      top50_true = mean(d$top50_true, na.rm = TRUE),
      top50_precision = mean(d$top50_precision, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
