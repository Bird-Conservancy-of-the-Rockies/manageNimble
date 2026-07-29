\name{runNimbleBlock}
\alias{runNimbleBlock}
\title{
  Run single block of a single chain for a NIMBLE model
}
\description{
  Base function called by \code{\link{runNimble}}'s generated worker script to initiate or continue sampling for one chain, one block at a time.
}
\usage{
  runNimbleBlock(mod.lst = NULL, comp.mcmc = NULL, n.iter = 1000, n.thin = 1,
            tmp.path, dump.file.path, SamplerSourcePath = NA,
            monitors2 = character(), n.thin2 = 1)
}
\arguments{
  \item{mod.lst}{A list containing: 1) code, output of nimbleCode(), 2) list of constants, 3) data list, 4) initial values list, 5) variables to monitor. Used on the first block of a chain. If this is NULL, 'comp.mcmc' argument must be provided.}
  \item{comp.mcmc}{Compiled Nimble MCMC object, already used for a prior block. Used on continuation blocks. If this is NULL, 'mod.lst' argument must be provided.}
  \item{n.iter}{Length of MCMC chain to run for this block.}
  \item{n.thin}{MCMC thinning rate. Only used (and only takes effect) on the first block of a chain, when the MCMC is configured; continuation blocks reuse the already-built sampler.}
  \item{tmp.path}{Temporary directory in which to save NIMBLE temporary files to facilitate cleanup.}
  \item{dump.file.path}{Save samples to this file. Written via a temp file and atomic rename, so a file visible under this name is always complete - callers reading block dumps never need to guard against a partial write.}
  \item{SamplerSourcePath}{Path to script file containing necessary lines to insert custom MCMC samplers. Must include necessary lines to remove default samplers and add desired samplers for target parameters. Default is to not provide a script file and thus use default samplers for all parameters. Requires the \pkg{nimbleHMC} package to be installed. Only used (and only takes effect) on the first block of a chain.}
  \item{monitors2}{Optional second set of monitored parameters, saved at their own (typically much coarser) thinning rate \code{n.thin2}. Configured on the first block of a chain; on continuation blocks it must still be supplied (identically) for the second set to be included in the dump, but \code{n.thin2} is ignored there (the rate is already baked into the compiled sampler). Default is no second set, in which case the dump file's contents are unchanged from before this feature existed.}
  \item{n.thin2}{Thinning rate for \code{monitors2}. Only used (and only takes effect) on the first block of a chain.}
}
\value{
  The compiled NIMBLE MCMC object (\code{comp.mcmc}), for use on the next call as this chain's \code{comp.mcmc} argument. As a side effect, writes \code{dump.file.path} containing \code{samp} (a matrix of posterior samples for \code{parameters}) and, if \code{monitors2} is non-empty, \code{samp2} (the second monitor set).
}
\author{
  Bryan L. Nuse and Quresh S. Latif, Bird Conservancy of the Rockies.

  Substantially revised with Claude (Anthropic) assistance during a full code review; see
  \code{NEWS.md} for what changed.
}

\keyword{nimble}
\keyword{runNimble}
