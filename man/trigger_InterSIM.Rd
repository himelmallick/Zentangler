\name{trigger_InterSIM}
\alias{trigger_InterSIM}
\title{Generate a Null InterSIM Multi-Omics Dataset}
\usage{
trigger_InterSIM(n)
}
\arguments{
\item{n}{Number of samples to generate.}
}
\value{
A list with \code{feature_table}, \code{sample_metadata}, and \code{feature_metadata}.
}
\description{
Calls \code{InterSIM::InterSIM()} to generate a null multi-omics dataset and formats it into a
simple list used internally by \code{gen_simmba()}.
}
\examples{
\dontrun{
pcl <- trigger_InterSIM(100)
}
}
