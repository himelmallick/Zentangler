\name{run_intersim_zentangler}
\alias{run_intersim_zentangler}
\title{Run InterSIM/SIMMBA Zentangler Benchmarks}
\usage{
run_intersim_zentangler(
  nrep = 10,
  nsample = 600,
  fusion_modes = c("early", "intermediate", "late"),
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  q_threshold = 0.25,
  top_n = 50L,
  seed = 1234,
  p.train = 0.7,
  ygen.mode = c("LM", "Friedman", "Friedman2"),
  outcome.type = c("continuous", "binary"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  b_inference = c("debiased_lasso", "refit"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  residualize = FALSE,
  covariates = NULL,
  bootstrap_repeats = 0L,
  bootstrap_id = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{nrep}{Number of simulation replicates.}
\item{nsample}{Number of samples per simulated dataset before train/test split.}
\item{fusion_modes}{Fusion modes to evaluate.}
\item{sis_n}{Number of mediators retained per view when \code{screen_method = "sis"}.}
\item{sis_rank}{Ranking criterion for SIS.}
\item{screen_method}{Screening method passed to \code{fit_multiview_parallel_zentangler()}. Use \code{"none"} to skip hard SIS.}
\item{q_threshold}{BH q-value threshold used to define active mediators.}
\item{top_n}{Number of top-ranked mediators used for top-k truth recovery.}
\item{seed}{Random seed.}
\item{p.train}{Training split proportion passed to \code{gen_simmba()}.}
\item{ygen.mode}{Outcome generation mode passed to \code{gen_simmba()}.}
\item{outcome.type}{Outcome type passed to \code{gen_simmba()}.}
\item{lambda_choice}{Cross-validated glmnet lambda choice.}
\item{b_inference}{B-path inference method.}
\item{debias_max_targets}{Maximum active targets for de-biased lasso inference.}
\item{coop_rho}{Agreement penalty strength for intermediate cooperative fusion.}
\item{residualize}{Logical. Residualize X, Y, and mediators against covariates.}
\item{covariates}{Optional covariate names in \code{colData(mae)}.}
\item{bootstrap_repeats}{Number of bootstrap repeats inside each fit.}
\item{bootstrap_id}{Optional bootstrap cluster ID column.}
\item{verbose}{Logical. Print progress.}
\item{...}{Additional arguments passed to \code{gen_simmba()}.}
}
\value{
A list with \code{summary}, \code{detail}, \code{fits}, \code{truth}, \code{simulation}, \code{errors}, and \code{settings}.
}
\description{
Generates InterSIM/SIMMBA \code{MultiAssayExperiment} simulations, fits Zentangler for one or more
fusion modes, and summarizes empirical truth recovery.
}
\examples{
\dontrun{
res <- run_intersim_zentangler(
  nrep = 3,
  nsample = 200,
  fusion_modes = c("early", "late"),
  sis_n = 100,
  screen_method = "sis",
  q_threshold = 0.25
)
res$summary
}
}
