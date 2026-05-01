\docType{package}
\name{Zentangler-package}
\alias{Zentangler}
\alias{Zentangler-package}
\title{Zentangler: Multi-View Parallel Mediation Analysis}
\description{
Zentangler provides MAE-first tools for generic multi-view parallel mediation analysis.
It works with \code{MultiAssayExperiment} objects, screens candidate mediators within each
assay view, fits selected mediators jointly, and returns ranked mediator evidence tables.
}
\details{
The main user-facing function is \code{fit_multiview_parallel_zentangler()}.
Simulation utilities are available through \code{gen_simmba()} and \code{trigger_InterSIM()}.
}
\seealso{
\code{fit_multiview_parallel_zentangler}, \code{gen_simmba}
}
