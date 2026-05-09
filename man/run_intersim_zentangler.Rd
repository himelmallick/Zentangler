\name{run_intersim_zentangler}
\alias{run_intersim_zentangler}
\title{Run InterSIM/SIMMBA Zentangler Benchmarks}
\usage{
run_intersim_zentangler(
  nrep = 10,
  nsample = 600,
  fusion_modes = c("early", "intermediate", "late"),
  method_preset = c("custom", "fast_lasso", "elastic_net", "longitudinal_maaslin2",
    "full_exploratory"),
  sis_n = NULL,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  a_stage_model = c("lm", "maaslin2"),
  q_threshold = 0.25,
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global", "within_view", "both"),
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
)
}
\arguments{
\item{nrep}{Number of simulation replicates.}
\item{nsample}{Number of samples per simulated dataset before train/test split.}
\item{fusion_modes}{Fusion modes to evaluate.}
\item{method_preset}{Optional named preset. Use \code{"custom"} to keep individual options.}
\item{sis_n}{Number of mediators retained per view when \code{screen_method = "sis"}.}
\item{sis_rank}{Ranking criterion for SIS.}
\item{screen_method}{Screening method passed to \code{fit_multiview_parallel_zentangler()}. Use \code{"none"} to skip hard SIS.}
\item{a_stage_model}{A-stage model: \code{"lm"} for HIMA-like univariate linear models or \code{"maaslin2"} for MaAsLin2 fixed/random-effect screening.}
\item{q_threshold}{One or more q-value thresholds used to define active mediators. When multiple thresholds are supplied, Zentangler fits the model once and evaluates truth recovery at each threshold.}
\item{fdr_method}{Multiple-testing correction used to compute q-values. Use \code{"BH"} for Benjamini-Hochberg or \code{"BY"} for Benjamini-Yekutieli.}
\item{fdr_scope}{Final FDR scope used in simulation evaluation. Use \code{"global"} to adjust across all mediators from all views, \code{"within_view"} to adjust separately inside each omics view, or \code{"both"} to report both from the same fitted model.}
\item{top_n}{Number of top-ranked mediators used for top-k truth recovery.}
\item{seed}{Random seed.}
\item{p.train}{Training split proportion passed to \code{gen_simmba()}.}
\item{ygen.mode}{Outcome generation mode passed to \code{gen_simmba()}.}
\item{outcome.type}{Outcome type passed to \code{gen_simmba()}.}
\item{lambda_choice}{Cross-validated glmnet lambda choice.}
\item{glmnet_alpha}{Glmnet mixing parameter for early and late fusion. Use \code{1} for lasso, values between \code{0} and \code{1} for elastic net, and \code{0} for ridge. Intermediate fusion currently uses the cooperative lasso-style update.}
\item{maaslin2_random_effect}{Optional MaAsLin2 random-effect column name.}
\item{maaslin2_normalization}{MaAsLin2 normalization method.}
\item{maaslin2_transform}{MaAsLin2 transform method.}
\item{maaslin2_analysis_method}{MaAsLin2 analysis method.}
\item{maaslin2_standardize}{Logical passed to MaAsLin2.}
\item{maaslin2_output_dir}{Optional directory for MaAsLin2 A-stage outputs.}
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
  glmnet_alpha = 0.75,
  q_threshold = seq(0.05, 0.25, by = 0.05),
  fdr_scope = "both"
)
res$summary
}
}
