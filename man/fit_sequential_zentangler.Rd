\name{fit_sequential_zentangler}
\alias{fit_sequential_zentangler}
\alias{zentangler_sequential_paths}
\title{Fit k-layer sequential Zentangler mediation}
\usage{
fit_sequential_zentangler(
  mae,
  x_var,
  y_var = NULL,
  stage_views = NULL,
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
  outcome_model = c("joint", "marginal"),
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global"),
  seed = 1,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
)

zentangler_sequential_paths(fit)
}
\arguments{
\item{mae}{A \code{MultiAssayExperiment}.}
\item{x_var}{Exposure column in \code{colData(mae)}.}
\item{y_var}{Outcome column in \code{colData(mae)} for Gaussian/binomial outcomes.}
\item{stage_views}{List of ordered mediator layers. Each element is a character vector of experiment/view names. If \code{NULL}, \code{view_names} order is used as one view per layer.}
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
\item{outcome_model}{\code{"joint"} fits terminal mediators jointly when feasible; \code{"marginal"} fits one terminal mediator at a time.}
\item{fdr_method}{FDR method for final path-level p-values.}
\item{fdr_scope}{Currently \code{"global"} only for path-level q-values.}
\item{seed}{Random seed.}
\item{duplicate_primary, fill_nonfinite_zero}{MAE extraction options.}
\item{fit}{Output from \code{fit_sequential_zentangler()}.}
}
\description{
Fits a sequential mediation screen for ordered mediator layers such as
\code{X -> M1 -> M2 -> ... -> Mk -> Y}. The first layer is screened by
association with the exposure. Adjacent mediator layers are screened by
mediator-to-mediator correlation, then retained links are refit by regression to
estimate directed link coefficients.
}
\details{
The returned path score is \code{a * prod(d) * b}, where \code{a} is the
first-layer exposure-to-mediator coefficient, \code{d} values are adjacent
mediator-to-mediator coefficients after correlation screening, and \code{b} is
the terminal mediator-to-outcome coefficient. The path-level \code{p_primary}
currently uses the maximum available leg p-value as a conservative screening
summary, followed by global BH/BY adjustment.
}
\value{
A list with settings, diagnostics, first-layer screens, transition links,
terminal effects, and \code{sequential_paths}. \code{zentangler_sequential_paths()}
returns the sequential path table.
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
}
}
