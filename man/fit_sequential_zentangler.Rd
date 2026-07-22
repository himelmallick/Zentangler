\name{fit_sequential_zentangler}
\alias{fit_sequential_zentangler}
\alias{zentangler_sequential_paths}
\alias{zentangler_sequential_edges}
\alias{zentangler_sequential_terminals}
\alias{zentangler_sequential_diagnostics}
\alias{summarize_sequential_zentangler}
\title{Fit sequential Zentangler routes}
\usage{
fit_sequential_zentangler(
  mae,
  x_var,
  y_var = NULL,
  stage_views = NULL,
  path_templates = NULL,
  view_names = NULL,
  assay_names = NULL,
  covariates = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  survival_time_var = NULL,
  survival_event_var = NULL,
  sis_n = 50L,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  cor_method = c("spearman", "pearson"),
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25,
  cor_fdr_method = c("BH", "BY"),
  residualize_links = TRUE,
  max_links_per_transition = 5000L,
  max_paths = 10000L,
  fusion_mode = c("early", "intermediate", "late"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  b_inference = c("debiased", "bootstrap"),
  debias_max_targets = 200L,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04,
  path_inference = c("model_based", "bootstrap_score"),
  bootstrap_repeats = 0L,
  bootstrap_ci_level = 0.95,
  bootstrap_seed = seed + 10000L,
  bootstrap_id = NULL,
  causal_inference = c("none", "bootstrap"),
  effect_x0 = NULL,
  effect_x1 = NULL,
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global"),
  seed = 1,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
)

zentangler_sequential_paths(fit)

zentangler_sequential_edges(fit)

zentangler_sequential_terminals(fit)

zentangler_sequential_diagnostics(fit)

summarize_sequential_zentangler(fit, q_threshold = 0.25)
}
\arguments{
\item{mae}{A \code{MultiAssayExperiment}.}
\item{x_var}{Exposure column in \code{colData(mae)}.}
\item{y_var}{Outcome column in \code{colData(mae)} for Gaussian/binomial outcomes.}
\item{stage_views}{List of ordered mediator layers. Each element is a character vector of experiment/view names. If \code{NULL}, \code{view_names} order is used as one view per layer.}
\item{path_templates}{Optional one-route list such as \code{list(route = c("species", "fecal_metabolites"))}. The route must contain at least two modalities. Run multiple routes in separate calls.}
\item{view_names}{Optional subset/order of MAE experiments.}
\item{assay_names}{Optional assay names passed to MAE extraction.}
\item{covariates}{Optional phenotype covariates.}
\item{y_family}{Outcome family: \code{"gaussian"}, \code{"binomial"}, or \code{"survival"}.}
\item{survival_time_var, survival_event_var}{Survival outcome columns.}
\item{sis_n}{Number of first-layer mediators to retain per first-layer view.}
\item{sis_rank}{Ranking for first-layer screen.}
\item{screen_method}{First-layer screen method.}
\item{cor_method}{Correlation method for adjacent mediator-layer screens.}
\item{min_abs_cor}{Minimum absolute correlation for adjacent-layer links.}
\item{cor_q_threshold}{FDR cutoff for adjacent-layer correlation links.}
\item{cor_fdr_method}{FDR method for adjacent-layer correlation links.}
\item{residualize_links}{Whether to residualize mediator-to-mediator correlations on exposure and covariates before screening.}
\item{max_links_per_transition}{Maximum retained links per adjacent transition.}
\item{max_paths}{Maximum full paths to keep while expanding across k layers.}
\item{fusion_mode}{Terminal mediator-to-outcome fusion mode: \code{"early"}, \code{"intermediate"}, or \code{"late"}.}
\item{lambda_choice}{Cross-validated glmnet lambda choice for terminal fusion.}
\item{glmnet_alpha}{Glmnet mixing parameter for early and late terminal fusion. Use 1 for lasso, 0.5 for elastic net, and 0 for ridge.}
\item{b_inference}{Terminal B-stage inference method. \code{"debiased"} uses the package's de-biased terminal B-stage approximation. \code{"bootstrap"} uses repeated B-stage refits and reports bootstrap terminal p-values.}
\item{debias_max_targets}{Maximum active terminal mediators to use in debiased B-stage calculations.}
\item{coop_rho}{Agreement penalty strength for intermediate terminal fusion.}
\item{coop_maxit}{Maximum cooperative updates for intermediate terminal fusion.}
\item{coop_tol}{Convergence tolerance for intermediate terminal fusion.}
\item{path_inference}{Path-level inference. \code{"model_based"} uses the conservative refit evidence assembled from the A-stage, transition-link, and terminal B-stage p-values. \code{"bootstrap_score"} uses bootstrap p-values for the complete sequential path score and requires \code{bootstrap_repeats > 0}.}
\item{bootstrap_repeats}{Number of path bootstrap refits. Use 0 to disable.}
\item{bootstrap_ci_level}{Bootstrap confidence level.}
\item{bootstrap_seed}{Bootstrap random seed.}
\item{bootstrap_id}{Optional phenotype column used for cluster bootstrap.}
\item{causal_inference}{Optional post-selection sequential causal-effect uncertainty mode. Use \code{"bootstrap"} to add bootstrap summaries for direct, indirect, and total sequential effects.}
\item{effect_x0,effect_x1}{Optional exposure contrast used for sequential causal effect summaries. If \code{NULL}, binary exposures use the lower/higher observed levels and continuous exposures use the 25th/75th percentiles.}
\item{fdr_method}{FDR method for final path-level p-values.}
\item{fdr_scope}{Currently \code{"global"} only for path-level q-values.}
\item{seed}{Random seed.}
\item{duplicate_primary, fill_nonfinite_zero}{MAE extraction options.}
\item{fit}{Output from \code{fit_sequential_zentangler()}.}
}
\description{
Fits one sequential mediation route such as \code{X -> M1 -> M2 -> ... -> Mk -> Y}.
The first layer is screened by association with the exposure.
Adjacent mediator layers are screened by mediator-to-mediator correlation, then
retained links are refit by regression to estimate directed link coefficients.
}
\details{
The returned path score is \code{a * prod(d) * b}, where \code{a} is the
first-layer exposure-to-mediator coefficient, \code{d} values are adjacent
mediator-to-mediator coefficients after correlation screening, and \code{b} is
the terminal mediator-to-outcome coefficient. The path-level \code{p_primary}
currently uses the maximum available leg p-value as a conservative screening
summary, followed by global BH/BY adjustment.

Terminal mediator-to-outcome p-values use \code{b_inference}. Complete path
evidence is controlled separately by \code{path_inference}. When
\code{bootstrap_repeats > 0}, Zentangler refits fixed selected paths on
bootstrap samples and adds path-score bootstrap means, intervals, selection
frequencies, and bootstrap p-values. Setting
\code{path_inference = "bootstrap_score"} uses those bootstrap path p-values as
the final \code{p_primary}; otherwise they are reported alongside the
model-based path evidence.

Sequential fits also expose a post-selection causal-effect stage, analogous to
the parallel API. After route discovery, Zentangler propagates the exposure
contrast through the retained ordered mediator chain, fits an outcome model on
\code{Y ~ X + active sequential mediators + covariates}, and reports
\code{effects}, \code{causal_effects}, and \code{effect_decomposition}
containing direct, indirect, and total-effect-style summaries. When
\code{causal_inference = "bootstrap"} and \code{bootstrap_repeats > 0}, those
sequential causal effects receive bootstrap means and intervals.

Zentangler distinguishes pure parallel mediation from one-route sequential
mediation. Use \code{fit_multiview_parallel_zentangler()} for \code{X -> M -> Y}
models without mediator-to-mediator links, and use
\code{fit_sequential_zentangler()} for an explicit ordered route with at least
two mediator modalities.

The terminal mediator-to-outcome stage uses \code{fusion_mode}. Early fusion
combines terminal mediators before fitting one penalized outcome model.
Intermediate fusion applies the package's cooperative learner across terminal
views and currently requires Gaussian outcomes with at least two terminal
views. Late fusion fits view-level terminal models and combines their scores.

The sequential API supports Gaussian, binomial, and survival outcomes. For
binomial outcomes, \code{y_var} should name a numeric
0/1 outcome and \code{b_inference = "debiased"} requests the
logistic B-stage approximation. For survival outcomes, set
\code{y_family = "survival"}, provide \code{survival_time_var} and
\code{survival_event_var}, and use early or late fusion with
\code{b_inference = "debiased"} or \code{"bootstrap"}.
}
\value{
A list with settings, diagnostics, first-layer screens, transition links,
terminal effects, and \code{sequential_paths}. \code{zentangler_sequential_paths()}
returns the sequential path table, \code{zentangler_sequential_edges()} returns
retained adjacent mediator links, \code{zentangler_sequential_terminals()}
returns terminal mediator-to-outcome coefficients, and
\code{summarize_sequential_zentangler()} returns model, threshold, diagnostic,
path, edge, terminal, and causal-effect summaries. Sequential fits also store
\code{effects}, \code{causal_effects}, and \code{effect_decomposition}.
}
\examples{
\dontrun{
fit_seq <- fit_sequential_zentangler(
  mae,
  x_var = "A",
  y_var = "Y",
  stage_views = list(
    microbiome = c("species", "kos"),
    metabolites = c("fecal_metabolites", "plasma_metabolites")
  ),
  sis_n = 50,
  cor_method = "spearman",
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25
)

head(zentangler_sequential_paths(fit_seq))

fit_binary <- fit_sequential_zentangler(
  mae,
  x_var = "A",
  y_var = "Y_binary",
  path_templates = list(route = c("species", "fecal_metabolites")),
  y_family = "binomial",
  b_inference = "debiased"
)

fit_survival <- fit_sequential_zentangler(
  mae,
  x_var = "A",
  y_var = NULL,
  path_templates = list(route = c("species", "fecal_metabolites")),
  y_family = "survival",
  survival_time_var = "time",
  survival_event_var = "status",
  b_inference = "debiased"
)
}
}
