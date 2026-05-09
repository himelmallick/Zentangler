\name{fit_multiview_parallel_zentangler}
\alias{fit_multiview_parallel_zentangler}
\title{Fit Generic Multi-View Parallel Mediation With MAE Input}
\usage{
fit_multiview_parallel_zentangler(
  mae,
  x_var,
  y_var,
  view_names = NULL,
  assay_names = NULL,
  method_preset = c("custom", "fast_lasso", "elastic_net", "longitudinal_maaslin2",
    "full_exploratory"),
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
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global", "within_view"),
  b_inference = c("debiased_lasso", "refit"),
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
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
)
}
\arguments{
\item{mae}{A \code{MultiAssayExperiment}. Each experiment is treated as one omics view.}
\item{x_var}{Name of the exposure column in \code{colData(mae)}. Must be numeric.}
\item{y_var}{Name of the outcome column in \code{colData(mae)}. Must be numeric.}
\item{view_names}{Optional character vector of MAE experiment names to use. Defaults to all experiments.}
\item{assay_names}{Optional named character vector of assay names, one per view. Defaults to the first assay in each experiment.}
\item{method_preset}{Optional named preset. \code{"custom"} preserves individual arguments. Other presets set coherent defaults for common workflows. See \code{zentangler_method_presets()}.}
\item{covariates}{Optional covariate column names in \code{colData(mae)}.}
\item{residualize}{Logical. If \code{TRUE}, residualize exposure, outcome, and mediators against covariates before fitting.}
\item{sis_n}{Number of mediators to retain per view after screening. If \code{NULL}, a HIMA-style default is used.}
\item{sis_rank}{Screening rank criterion: absolute a-path estimate or a-path p-value.}
\item{screen_method}{Mediator screening method. \code{"sis"} keeps the top \code{sis_n} mediators per view; \code{"none"} keeps all mediators that survive numeric and variance filtering.}
\item{a_stage_model}{A-stage model. \code{"lm"} uses HIMA-like univariate linear models. \code{"maaslin2"} uses MaAsLin2 fixed/random-effect models for the exposure-to-mediator leg.}
\item{maaslin2_random_effect}{Optional MaAsLin2 random-effect column name in \code{colData(mae)}, for example a subject ID in longitudinal data.}
\item{maaslin2_normalization}{MaAsLin2 normalization method.}
\item{maaslin2_transform}{MaAsLin2 transform method.}
\item{maaslin2_analysis_method}{MaAsLin2 analysis method.}
\item{maaslin2_standardize}{Logical passed to MaAsLin2.}
\item{maaslin2_output_dir}{Optional directory where MaAsLin2 A-stage outputs should be written. If \code{NULL}, temporary directories are used.}
\item{fusion_mode}{Y-stage fusion mode: early, intermediate, or late.}
\item{y_family}{Outcome family: Gaussian or binomial. Intermediate cooperative fusion currently supports Gaussian outcomes.}
\item{lambda_choice}{Cross-validated glmnet lambda choice.}
\item{glmnet_alpha}{Glmnet mixing parameter for early and late fusion. Use \code{1} for lasso, values between \code{0} and \code{1} for elastic net, and \code{0} for ridge. Intermediate fusion currently uses the cooperative lasso-style update.}
\item{fdr_method}{Multiple-testing correction used to compute q-values. Use \code{"BH"} for Benjamini-Hochberg or \code{"BY"} for Benjamini-Yekutieli.}
\item{fdr_scope}{Final FDR scope. Use \code{"global"} to compute \code{q_primary} across all mediators from all views, or \code{"within_view"} to compute \code{q_primary} separately inside each view. Both \code{q_primary_global} and \code{q_primary_within_view} are retained in the output table.}
\item{b_inference}{B-path inference method.}
\item{debias_max_targets}{Maximum number of active B-path targets for de-biased lasso inference.}
\item{coop_rho}{Agreement penalty strength for intermediate cooperative fusion.}
\item{coop_maxit}{Maximum cooperative fusion iterations.}
\item{coop_tol}{Convergence tolerance for cooperative fusion.}
\item{seed}{Random seed.}
\item{return_fits}{Logical. If \code{TRUE}, include fitted model objects.}
\item{bootstrap_repeats}{Number of bootstrap refits. Use 0 to disable.}
\item{bootstrap_ci_level}{Bootstrap confidence level.}
\item{bootstrap_seed}{Bootstrap seed.}
\item{bootstrap_id}{Optional cluster or subject ID column in \code{colData(mae)} for cluster bootstrap.}
\item{duplicate_primary}{How to handle duplicate primary sample IDs after MAE mapping.}
\item{fill_nonfinite_zero}{Logical. If \code{TRUE}, replace non-finite assay values with zero before filtering.}
}
\value{
A list with run settings, diagnostics, aligned sample IDs, per-view mediator tables, cleaner mediator aliases
(\code{mediators_all}, \code{mediators_active}, \code{mediators_top}), view/model summaries,
exploratory effect summaries, optional bootstrap summaries, and optional model fits.
}
\description{
Fits the Zentangler parallel multi-view mediation workflow using a \code{MultiAssayExperiment} input.
For each mediator in each view, the method estimates an exposure-to-mediator a-path, optionally screens
mediators within view, fits the retained mediators jointly in a Y-stage model, and ranks mediators using \code{a * b}.
}
\examples{
\dontrun{
fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  sis_n = 50,
  fusion_mode = "early",
  glmnet_alpha = 0.5,
  b_inference = "debiased_lasso",
  bootstrap_repeats = 100
)
head(fit$combined_mediators)
}
}
