\name{mcmcOutputSubset}
\alias{mcmcOutputSubset}
\title{
  Subset mcmcOutput object
}
\description{
  Utility function that subsets mcmcOutput object. Useful to implement prior to summarizing MCMC samples stored as mcmcOutput objects when not all parameters sampled need to be summarized.
}
\usage{
  mcmcOutputSubset(mcmcOutput, par.keep = c(), par.drop = c())
}
\arguments{
  \item{mcmcOutput}{Object containing MCMC samples compiled as an mcmcOutput object.}
  \item{par.keep}{Character vector of parameter names to be summarized from mcmcOutput. Must be specified if par.ignore = NULL.}
  \item{par.drop}{Character vector of parameter names to be trimmed from mcmcOutput. Must be specified if par.summarize = NULL. Also, any parameter listed here that is also listed in par.summarize will be summarized.}
}
\value{Object of class 'mcmcOutput' that is a subsetted version of the input object.}
\author{
  Bryan L. Nuse and Quresh S. Latif, Bird Conservancy of the Rockies
}

\keyword{nimble}
\keyword{mcmcOutput}
