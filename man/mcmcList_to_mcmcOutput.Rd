\name{mcmcList_to_mcmcOutput}
\alias{mcmcList_to_mcmcOutput}
\title{Converts mcmcList to mcmcOutput}
\description{Converts mcmcList object (i.e., list of matrices containing mcmc samples from parallel chains) to a list with two elements, an mcmcOutput object and a summary table of all parameters.}
\usage{
mcmcList_to_mcmcOutput(mcmcList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcList}{A list whose elements are parallel matrices (must be equal in dimension), each of which contains saved mcmc samples for a series of model parameters. Rows are samples and columns are parameters. Column names should meaningfully identify the parameters for this function to be useful.}
}
\value{List of two elements: 1) an mcmcOutput object and 2) a summary table of parameters, including Rhat and n.effective.}
\author{
Quresh S. Latif
}
\note{R packages required: coda, mcmcOutput, dplyr}
