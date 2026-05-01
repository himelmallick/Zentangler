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
  covariates = NULL,
  residualize = FALSE,
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  fusion_mode = c("early", "intermediate", "late"),
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
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
\item{covariates}{Optional covariate column names in \code{colData(mae)}.}
\item{residualize}{Logical. If \code{TRUE}, residualize exposure, outcome, and mediators against covariates before fitting.}
\item{sis_n}{Number of mediators to retain per view after screening. If \code{NULL}, a HIMA-style default is used.}
\item{sis_rank}{Screening rank criterion: absolute a-path estimate or a-path p-value.}
\item{fusion_mode}{Y-stage fusion mode: early, intermediate, or late.}
\item{y_family}{Outcome family: Gaussian or binomial. Intermediate cooperative fusion currently supports Gaussian outcomes.}
\item{lambda_choice}{Cross-validated glmnet lambda choice.}
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
A list with run settings, aligned sample IDs, per-view mediator tables, a combined mediator ranking table,
exploratory effect summaries, optional bootstrap summaries, and optional model fits.
}
\description{
Fits the Zentangler parallel multi-view mediation workflow using a \code{MultiAssayExperiment} input.
For each mediator in each view, the method estimates an exposure-to-mediator a-path, screens mediators
within view, fits selected mediators jointly in a Y-stage model, and ranks mediators using \code{a * b}.
}
\examples{
\dontrun{
fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  sis_n = 50,
  fusion_mode = "early",
  b_inference = "debiased_lasso",
  bootstrap_repeats = 100
)
head(fit$combined_mediators)
}
}
