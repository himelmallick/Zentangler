\name{zentangler_output_helpers}
\alias{zentangler_method_presets}
\alias{zentangler_all_mediators}
\alias{zentangler_active_mediators}
\alias{zentangler_top_mediators}
\alias{zentangler_view_summary}
\alias{zentangler_model_summary}
\alias{summarize_zentangler}
\title{Zentangler Presets and Output Helpers}
\usage{
zentangler_method_presets()

zentangler_all_mediators(fit)

zentangler_active_mediators(fit, q_threshold = 0.25, q_col = "q_primary")

zentangler_top_mediators(
  fit,
  n = 20L,
  order_by = c("abs_score", "q_primary", "p_primary", "score"),
  decreasing = NULL
)

zentangler_view_summary(fit, q_threshold = 0.25, q_col = "q_primary")

zentangler_model_summary(fit, q_threshold = 0.25)

summarize_zentangler(fit, q_threshold = 0.25)
}
\arguments{
\item{fit}{A fitted object returned by \code{fit_multiview_parallel_zentangler()}.}
\item{q_threshold}{BH q-value threshold used to define active mediators.}
\item{q_col}{Name of the q-value column used for active mediator calls.}
\item{n}{Number of mediators to return.}
\item{order_by}{Column used to order mediators.}
\item{decreasing}{Logical ordering direction. If \code{NULL}, score columns are sorted decreasing and p/q columns increasing.}
}
\description{
Convenience helpers for inspecting Zentangler method presets and fitted model outputs.
}
\value{
\code{zentangler_method_presets()} returns a data frame describing named presets. Mediator helpers return data frames. \code{summarize_zentangler()} returns a list containing model summary, view summary, top mediators, and active mediators.
}
\examples{
\dontrun{
zentangler_method_presets()
fit <- fit_multiview_parallel_zentangler(mae, x_var = "X", y_var = "Y")
zentangler_model_summary(fit)
zentangler_view_summary(fit)
zentangler_top_mediators(fit, n = 10)
zentangler_active_mediators(fit, q_threshold = 0.25)
}
}
