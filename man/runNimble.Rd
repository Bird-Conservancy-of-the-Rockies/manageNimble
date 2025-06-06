\name{runNimble}
\alias{runNimble}
\title{
  NIMBLE wrapper with parallel processing and automated sampling until convergence
}
\description{
  Implements Bayesian model fitting using the NIMBLE R packaged in parallelized clusters across multiple processors. Additionally, convergence and sampling metrics (i.e., Rhat and n_effective) are checked periodically, and sampling continues until these metrics meet specified minimum levels. Saved samples are periodically dumped to the hard drive to avoid overloading RAM with successive rounds of resuming sampling.
  
  ***Note: This function relies on submitting Linux commands via the processx R package, and thus will not operate in Windows.
}
\usage{
  runNimble(model.path, inits, data, constants, parameters, max.samples.saved,
           par.ignore = c(), par.dontign = c(), par.fuzzy.track = c(), fuzzy.threshold = 0.05,
           nc = 2, ni = 2000, nb = 0.5, nt = 10, mod.nam = "mod",
           rtrn.model = F, sav.model = T, check.freq = NULL,
           Rht.required = 1.1, neff.required = 100, max.tries = 10,
           dump.path = "dump", SamplerSourcePath = NA, delete.blocks = TRUE)
}
\arguments{
  \item{model.path}{File path to Nimble model code, which should be set up to be sourced whereby code is contained within `model <<- nimble::nimbleCode({})`.}
  \item{inits}{List containing all values for initializing Nimble model. Names and dimensions of listed elements should match parameters in Nimble model code.}
  \item{data}{List containing all data objects required to fit specified NIMBLE model.}
  \item{constants}{List containing all constants required to fit specified NIMBLE model.}
  \item{parameters}{String vector containing names of all parameters for which estimates are to be retained from sampler.}
  \item{max.samples.saved}{Integer that sets max samples to save across chains in final model object. Without automated convergence checking enabled, model run ends once saved samples have met or exceeded this limit. With automated convergence checking, additional thinning of posterior samples (i.e., after samples are recovered from processing clusters) is implemented as necessary to limit the size of the model object to this level.}
  \item{par.ignore}{String vector listing parameters to be ignored when calculating convergence and sampling metrics. All parameters with names containing strings listed here will be ignored. Default is an empty vector, whereby all parameters will be considered for calculating Rhat. Ignored parameters will be excluded from summary output table.}
  \item{par.dontign}{String vector listing parameters that should not be ignored when calculating convergence and sampling metrics. All parameters with names containing strings listed here will not be ignored. This argument is unnecessary if 'par.ignore' is not specified. Only relevant if automated convergence checks are enabled.}
  \item{par.fuzzy.track}{String vector listing parameters to designate for fuzzy evaluation of convergence. If more than (fuzzy.threshold * number of parameters) Rhats > Rht.required, sampling will continue. Default is an empty vector, whereby no fuzzy evaluation of convergence will occur. Only relevant if automated convergence checks are enabled.}
  \item{fuzzy.threshold}{Threshold for fuzzy evaluation of convergence (only relevant if length(par.fuzzy.track) > 0). Only relevant if automated convergence checks are enabled.}
  \item{nc}{Number of MCMC chains to run. Note available processing cores must equal or exceed 'nc'. Default = 2.}
  \item{ni}{Number of iterations per sampling block. Default = 2000.}
  \item{nb}{If < 1, proportion of MCMC chains to discard as burn-in. If > 1, chain length before thinning to discard as burn-in.}
  \item{nt}{Base level of thinning to implement within parallelized cluster. Default = 10. Setting this argument will reduce load on RAM for large models.}
  \item{mod.nam}{Character string used as a file label for saved model outputs.}
  \item{rtrn.model}{Logical indicating whether to return model output to R environment. Default = FALSE. Must be set to TRUE if sav.model = FALSE.}
  \item{sav.model}{Logical indicating whether to save model output to disk. Default = TRUE. Must be set to TRUE if rtrn.model = FALSE.}
  \item{check.freq}{How long to wait in number of sampling blocks between successive checks of convergence criteria. If NULL (default), automated convergence checking is disabled, and sampling will be capped based on max.samples.saved.}
  \item{Rht.required}{Maximum numeric value for Rhat for any parameter required to end sampling.}
  \item{neff.required}{Minimum value for n_effective for any parameter required to end sampling.}
  \item{max.tries}{The maximum number of times to check convergence, after which runNimble will quit with warning that convergence was not reached before reaching the specified cap. Default is 10. It is possible but not recommended to set to NULL for no constraint on number of checks.}
  \item{dump.path}{Directory in which to save intermediate files and sample blocks generated by the function. Necessary to specify if multiple analyses are being run in parallel within the same user profile for parallelized nodes to avoid over-writing each other's files. Do not include '/' at the end.}
  \item{SamplerSourcePath}{Path to script file containing necessary lines to insert custom MCMC samplers. Must include necessary lines to remove default samplers and add desired samplers for target parameters. If this works, an example may be added below. Default is to not provide a script file and thus use default samplers for all parameters.}
  \item{delete.blocks}{Whether to delete dump.path directory containing saved sample blocks once run is complete. Default is to delete.}
}
\value{If assigned to object in global environment and rtrn.model = T, model output is returned. Else, model output is saved to file. Regardless, returned value consists of a list with three elements:
  \item{mcmcOutput}{posterior samples as returned by mcmcOutput::mcmcOutput}
  \item{summary}{data frame containing summaries of posterior estimates and convergence criteria (Rhat and n_effective) for all parameters in mcmcOutput}
  \item{mcmc.info}{A vector with elemets: nchains = number of MCMC sampling chains, niterations = realized pre-thinning chain length including burnin, burnin = realized pre-thinning length of each chain discarded as burn-in, nthin = realized thinning expressed as the ith sample retained. Note that nthin can be a non-whole number if max.samples.saved is specified, in which case it is the mean interval between retained samples.}
  }
\author{
  Quresh S. Latif and Bryan L. Nuse, Bird Conservancy of the Rockies
}

\keyword{nimble}
\keyword{mcmcOutput}
\keyword{parallel processing}
