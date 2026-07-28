\name{mcmcOutputSubset}
\alias{mcmcOutputSubset}
\title{
  Subset mcmcOutput object
}
\description{
  Utility function that subsets an mcmcOutput object. Useful to implement prior to summarizing MCMC samples stored as mcmcOutput objects when not all parameters sampled need to be summarized.

  Parameter name matching is by substring, not exact match: a name in \code{par.keep}/\code{par.drop} matches any column whose name contains it (e.g. \code{"p"} matches both \code{p[1]} and \code{beta.p[1]}). If both are supplied, the result is the union of \code{par.keep} and everything not matched by \code{par.drop} - \code{par.keep} does not restrict the output to only those columns when \code{par.drop} is also given.
}
\usage{
  mcmcOutputSubset(mcmcOutput, par.keep = c(), par.drop = c())
}
\arguments{
  \item{mcmcOutput}{Object containing MCMC samples compiled as an mcmcOutput object.}
  \item{par.keep}{Character vector of parameter names to keep regardless of \code{par.drop}. At least one of \code{par.keep}/\code{par.drop} must be specified.}
  \item{par.drop}{Character vector of parameter names to remove, except any also matched by \code{par.keep}. At least one of \code{par.keep}/\code{par.drop} must be specified.}
}
\value{Object of class 'mcmcOutput' that is a subsetted version of the input object.}
\author{
  Bryan L. Nuse and Quresh S. Latif, Bird Conservancy of the Rockies
}

\keyword{nimble}
\keyword{mcmcOutput}
