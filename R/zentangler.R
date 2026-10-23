# Zentangler: multi-view mediation analysis for MultiAssayExperiment data
#
# Exported API
#   zentangler(model = "parallel" | "sequential")   main entry point
#   fit_multiview_parallel_zentangler()  X -> M_j (view k) -> Y
#   fit_sequential_zentangler()          X -> M1 -> ... -> Mk -> Y along a route
#   zentangler_mediators(), zentangler_paths(), zentangler_effects()
#   summary() / print() methods for both fit classes
#   run_negative_controls(), run_direction_check()
#
# File map
#   zentangler.R           zentangler() entry point (this file)
#   00-utils.R             validation, study design, sample alignment, FDR (HDMT-style)
#   01-mae-input.R         MultiAssayExperiment extraction, log / CLR transforms
#   02-results.R           accessors, summary() and print() methods
#   03-a-stage-screening.R A-path estimation and SIS / MaAsLin2 / Tweedieverse screening
#   04-y-stage-fusion.R    outcome model: multiview cooperative learning (main,
#                          rho in [0, 1]); IntegratedLearner-style fusion (supplementary)
#   05-b-path-inference.R  de-biased B-path inference, mediator tables
#   06-causal-effects.R    g-computation causal effects (parallel and sequential)
#   07-bootstrap.R         pipeline and fixed-set effect bootstraps
#   08-fit-parallel.R      parallel options, engine and public fit
#   09-fit-sequential.R    sequential route discovery and public fit
#   10-controls.R          outcome-permutation negative controls, direction checks

#' Zentangler multi-view mediation analysis
#'
#' One entry point for both Zentangler models.
#'
#' * `model = "parallel"`: every feature of every view is a candidate mediator,
#'   `X -> M_j -> Y` (see [fit_multiview_parallel_zentangler()] for all options).
#' * `model = "sequential"`: mediators chained along an ordered route of views,
#'   `X -> M1 -> ... -> Mk -> Y` (see [fit_sequential_zentangler()]).
#'
#' @param mae A `MultiAssayExperiment`.
#' @param model `"parallel"` (default) or `"sequential"`. If `route` is given
#'   and `model` is not, the sequential model is used.
#' @param x_var Exposure column in `colData(mae)`.
#' @param y_var Outcome column (`NULL` for survival outcomes).
#' @param route Character vector of views for the sequential model.
#' @param ... Further arguments for the chosen model, e.g. `covariates`,
#'   `y_family`, `transform`, `fusion`, `rho`, `b_inference`, `fdr_joint`,
#'   `bootstrap_repeats`.
#' @return A `zentangler_fit` (parallel) or `zentangler_sequential_fit`.
#' @examples
#' \dontrun{
#' fit <- zentangler(mae, x_var = "X", y_var = "Y", covariates = "age",
#'                   transform = c(species = "clr"))
#' summary(fit)
#' seq_fit <- zentangler(mae, model = "sequential", x_var = "X", y_var = "Y",
#'                       route = c("species", "metabolites"))
#' zentangler_paths(seq_fit, q_threshold = 0.1)
#' }
#' @export
zentangler <- function(mae, model = c("parallel", "sequential"), x_var = NULL, y_var = NULL, route = NULL, ...) {
  model <- zentangler_resolve_model(model, route, missing(model))
  if (identical(model, "parallel")) {
    fit_multiview_parallel_zentangler(mae = mae, x_var = x_var, y_var = y_var, ...)
  } else {
    fit_sequential_zentangler(mae = mae, x_var = x_var, y_var = y_var, route = route, ...)
  }
}

zentangler_resolve_model <- function(model, route, model_missing) {
  if (model_missing && !is.null(route)) return("sequential")
  model <- match.arg(as.character(model)[1], c("parallel", "sequential"))
  if (identical(model, "parallel") && !is.null(route)) {
    stop("route is only used with model = 'sequential'.", call. = FALSE)
  }
  if (identical(model, "sequential") && is.null(route)) {
    stop("model = 'sequential' needs route = c(view1, view2, ...).", call. = FALSE)
  }
  model
}
