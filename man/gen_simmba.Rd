\name{gen_simmba}
\alias{gen_simmba}
\title{Generate SIMMBA MultiAssayExperiment Simulation Data}
\usage{
gen_simmba(
  nsample,
  snr = 1,
  p.train = 0.7,
  de.prob = rep(0.1, 3),
  de.downProb = rep(0.5, 3),
  de.facLoc = rep(1, 3),
  de.facScale = rep(0.4, 3),
  mediation.spike = TRUE,
  TIE = rep(1/3, 3),
  n.pathways = NULL,
  dirPathways_a = rep(1, 3),
  dirPathways_b = rep(1, 3),
  effect_type = rep("spiked", 3),
  ygen.mode = c("LM", "Friedman", "Friedman2"),
  outcome.type = c("continuous", "binary", "survival"),
  surv.hscale = 1,
  cens.lower = 1,
  cens.upper = 3,
  nrep = 100,
  seed = 1234,
  return.pcl = FALSE
)
}
\arguments{
\item{nsample}{Number of samples to simulate.}
\item{snr}{Signal-to-noise ratio for continuous outcomes.}
\item{p.train}{Training split proportion.}
\item{de.prob, de.downProb, de.facLoc, de.facScale}{Simulation controls inherited from the SIMMBA-style generator.}
\item{mediation.spike}{Logical. If \code{TRUE}, spike mediation effects into selected features.}
\item{TIE}{Requested total indirect effect per generated modality.}
\item{n.pathways}{Number of true mediating features per modality.}
\item{dirPathways_a, dirPathways_b}{Directions of exposure-to-mediator and mediator-to-outcome effects.}
\item{effect_type}{Mediation truth type for each modality.}
\item{ygen.mode}{Outcome generation mode.}
\item{outcome.type}{Outcome type.}
\item{surv.hscale, cens.lower, cens.upper}{Survival outcome controls.}
\item{nrep}{Number of simulation replicates.}
\item{seed}{Random seed.}
\item{return.pcl}{Logical. If \code{TRUE}, also return legacy InterSIM-like raw objects.}
}
\value{
A list containing \code{trainMae}, \code{testMae}, \code{truthDat}, and simulation settings.
}
\description{
Generates multi-view simulation data as \code{MultiAssayExperiment} objects for evaluating Zentangler.
The function uses InterSIM to create baseline multi-omics data, adds mediation-style spike-in effects,
and stores the result as MAE train/test objects.
}
\examples{
\dontrun{
sim <- gen_simmba(nsample = 100, nrep = 1)
mae_train <- sim$trainMae[[1]]
truth <- sim$truthDat[[1]]
}
}
