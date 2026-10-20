\name{runNimble}
\alias{runNimble}
\title{
  NIMBLE wrapper with parallel processing and automated sampling until convergence
}
\description{
  Implements Bayesian model fitting using the NIMBLE R package in parallelized clusters across multiple processors. Additionally, convergence and sampling metrics (i.e., Rhat and n_effective) are checked periodically, and sampling continues until these metrics meet specified minimum levels. Saved samples are periodically dumped to the hard drive to avoid overloading RAM with successive rounds of resuming sampling.

  ***Note: This function relies on submitting Linux commands via the processx R package and GNU \code{parallel}, and thus will not operate on Windows.
}
\usage{
  runNimble(model.path, inits, data, constants, parameters, max.samples.saved,
           par.ignore = c(), par.dontign = c(), par.fuzzy.track = c(), fuzzy.threshold = 0.05,
           nc = 2, ni = 2000, nb = 0.5, nt = 1, mod.nam = "mod",
           rtrn.model = F, sav.model = T, check.freq = NULL,
           Rht.required = 1.1, neff.required = 100, max.tries = 10,
           dump.path = "dump", SamplerSourcePath = NA, delete.blocks = TRUE,
           parameters2 = character(), nt2 = NULL,
           consolidate.monitors2 = TRUE, seed = NULL)
}
\arguments{
  \item{model.path}{File path to Nimble model code, which should be set up to be sourced whereby code is contained within `model <<- nimble::nimbleCode({})`.}
  \item{inits}{Initial values, in one of three forms:
    \enumerate{
      \item A single named list, used by every chain. Names and dimensions of its elements should match nodes in the NIMBLE model code.
      \item An unnamed list of \code{nc} named lists, one per chain, e.g. \code{list(inits1, inits2, inits3)}.
      \item A function returning a named list. It is called once per chain in your R session, before any worker launches - as \code{f(chn)}, with the chain number, if it takes an argument, otherwise as \code{f()} - so it can refer to objects in your session, e.g. \code{function() list(z = zst, mu = rnorm(1))}. With \code{seed}, each chain's call draws from its own reproducible random-number stream, separate from the one the chain samples from.
    }
    The form is read from the \code{names} attribute: a fully named list (or an empty list, or \code{NULL}) is a single inits list; a list with no names is per-chain and must have length \code{nc}; a partly named list is an error. Any inits list with a list-valued element is also an error, because NIMBLE would silently ignore it - the usual cause is a per-chain list passed with names, \code{list(chain1 = ..., chain2 = ...)}. All of this is checked before anything is written to \code{dump.path}. Per-chain and generated inits are recorded in the saved object (see Value).}
  \item{data}{List containing all data objects required to fit specified NIMBLE model.}
  \item{constants}{List containing all constants required to fit specified NIMBLE model.}
  \item{parameters}{String vector containing names of all parameters for which estimates are to be retained from sampler.}
  \item{max.samples.saved}{Integer that sets max samples to save \strong{per chain} in the final model object (the total across chains is \code{nc * max.samples.saved}). Without automated convergence checking enabled, model run ends once saved samples have met or exceeded this limit; it is required (not optional) in that mode, since it is the sole stopping criterion. With automated convergence checking, additional thinning of posterior samples (i.e., after samples are recovered from processing clusters) is implemented as necessary to limit the size of the model object to this level.}
  \item{par.ignore}{String vector listing parameters to be ignored when calculating convergence and sampling metrics. All parameters with names \emph{containing} strings listed here will be ignored - this can match more than intended (e.g. \code{par.ignore = "p"} also matches \code{beta.p} and \code{mu.p}); use \code{par.dontign} to rescue any of those. Default is an empty vector, whereby all parameters will be considered for calculating Rhat. Ignored parameters will be excluded from summary output table.}
  \item{par.dontign}{String vector listing parameters that should not be ignored when calculating convergence and sampling metrics, even though their names contain one of the strings in \code{par.ignore}. Unnecessary if \code{par.ignore} is not specified. Only relevant if automated convergence checks are enabled.}
  \item{par.fuzzy.track}{String vector listing parameters to designate for fuzzy evaluation of convergence. These parameters are excluded from the strict Rhat/n_effective requirement and instead governed by \code{fuzzy.threshold}: sampling continues only if more than \code{fuzzy.threshold} of them fail Rhat, fail n_effective, or are missing altogether (NA Rhat or n_effective is always fatal, for these parameters as for any other - see \code{\link{checkNimble}}). Default is an empty vector, whereby no fuzzy evaluation of convergence will occur. Only relevant if automated convergence checks are enabled.}
  \item{fuzzy.threshold}{Threshold for fuzzy evaluation of convergence (only relevant if length(par.fuzzy.track) > 0). Only relevant if automated convergence checks are enabled.}
  \item{nc}{Number of MCMC chains to run. Note available processing cores must equal or exceed 'nc'. Default = 2.}
  \item{ni}{Number of iterations per sampling block. Default = 2000.}
  \item{nb}{If < 1, proportion of MCMC chains to discard as burn-in. If >= 1, chain length before thinning to discard as burn-in. When proportional, burn-in is recomputed against the run's current length at every convergence check, so it grows (and previously-retained blocks can become burn-in) as sampling continues.}
  \item{nt}{Base level of thinning to implement within parallelized cluster. Default = 1. Increasing this argument reduces the size of retained samples and the load on RAM for large models.}
  \item{mod.nam}{Character string used as a file label for saved model outputs, and as the name under which the result is assigned in the global environment when \code{rtrn.model = TRUE}.}
  \item{rtrn.model}{Logical indicating whether to return model output to R environment (assigned under the name \code{mod.nam}). Default = FALSE. Must be set to TRUE if sav.model = FALSE.}
  \item{sav.model}{Logical indicating whether to save model output to disk. Default = TRUE. Must be set to TRUE if rtrn.model = FALSE.}
  \item{check.freq}{How long to wait in number of sampling blocks between successive checks of convergence criteria. If NULL (default), automated convergence checking is disabled, and sampling will be capped based on max.samples.saved.}
  \item{Rht.required}{Maximum numeric value for Rhat for any parameter required to end sampling.}
  \item{neff.required}{Minimum value for n_effective for any parameter required to end sampling.}
  \item{max.tries}{The maximum number of times to check convergence, after which runNimble will quit with warning that convergence was not reached before reaching the specified cap. Default is 10; must be at least 1. It is possible but not recommended to set to NULL for no constraint on number of checks.}
  \item{dump.path}{Directory in which to save intermediate files and sample blocks generated by the function. Necessary to specify if multiple analyses are being run in parallel within the same user profile for parallelized nodes to avoid over-writing each other's files. Do not include '/' at the end.}
  \item{SamplerSourcePath}{Path to script file containing necessary lines to insert custom MCMC samplers. Must include necessary lines to remove default samplers and add desired samplers for target parameters. Default is to not provide a script file and thus use default samplers for all parameters. Requires the \pkg{nimbleHMC} package to be installed.}
  \item{delete.blocks}{Whether to delete dump.path directory containing saved sample blocks once run is complete. Default is to delete. Worker logs (see Details) are copied out to \code{paste0(mod.nam, "_logs")} before deletion regardless of this setting.}
  \item{parameters2}{Optional second set of monitored parameters, saved at their own (typically much coarser) thinning rate \code{nt2} instead of \code{nt}. Intended for quantities too large to retain at the primary rate - latent state arrays, for example - which are needed for post hoc calculations but not for convergence diagnostics. They take no part in convergence checking, and are gathered once at the end of the run (before \code{delete.blocks} removes the sample blocks) and always written to \code{paste0(mod.nam, "_monitors2.rds")} as \code{list(out, iter.key)} (an \code{[ndraw, nparam]} matrix with chains stacked, plus a \code{data.frame(chn, iter)} aligned with its rows) - see also \code{consolidate.monitors2}. The second set's kept draws are derived directly from the primary set's own retained iterations, which guarantees the ratio of retained second-set rows to primary-set rows is always \code{nt2 / nt} (regardless of how many retry blocks a chain needed) and that every retained second-set row traces back to a specific primary-set row - join \code{mod$iter.key} (see \code{Value} below) against \code{mod$monitors2$iter.key} on \code{c("chn", "iter")} to recover that correspondence. Must not overlap with \code{parameters}. Default is no second set, in which case nothing about existing behaviour changes.}
  \item{nt2}{Thinning rate for \code{parameters2}. Defaults to \code{ni} (one retained draw per block per chain) when \code{parameters2} is non-empty; set smaller for more draws at proportionally more disk. Must be a positive integer multiple of \code{nt} when \code{parameters2} is non-empty (validated up front, before any workers are launched) - required for the guarantee above. Ignored (fixed at 1) when \code{parameters2} is empty.}
  \item{consolidate.monitors2}{When \code{TRUE} (default) and \code{parameters2} is non-empty, the gathered second monitor set is also attached to the saved/returned model object as \code{mod$monitors2}, so a single \code{R.utils::loadObject(mod.nam)} call (or \code{rtrn.model = TRUE}) gives access to both monitor sets together; the standalone \code{paste0(mod.nam, "_monitors2.rds")} file is still written either way. The tradeoff is memory/load time: loading the primary object now also loads the second monitor set, which can be large - the reason \code{parameters2} exists as a separate set in the first place. Set to \code{FALSE} to keep the second set out of the primary object, available only via its own file - e.g. for callers who want the smaller/faster load. Has no effect when \code{parameters2} is empty; \code{mod} never gets a \code{monitors2} element in that case.}
  \item{seed}{\code{NULL} (default) leaves the random-number generator alone: chains differ, but runs are not reproducible. A single whole number makes the run reproducible: one non-overlapping L'Ecuyer-CMRG random-number stream per chain is derived from it (with \code{parallel::nextRNGStream}), and each worker installs its chain's stream before the model file is sourced or the model is built, compiled or sampled; continuation blocks carry on in the same stream. Your session's own random-number state is left as it was. Recorded in the saved object (see Value). See Details for exactly what is reproducible.}
}
\details{
  Coordination between the parent process and the parallel workers happens via plain text files under \code{dump.path} (\code{runNimbleDirective.txt} and \code{block<n>Status.txt}), read and written atomically (temp file + rename) by both sides to avoid torn reads.

  Each worker's console output and exit status is captured under \code{dump.path/worker_logs} and \code{dump.path/worker_joblog.txt} (via GNU \code{parallel}'s \code{--results}/\code{--joblog}), and any R-level error inside a worker is additionally logged to \code{dump.path/chn<n>_error.log} with a timestamp and the underlying error message. The parent checks for these on every polling iteration and stops with a pointer to the relevant file(s) if a worker has died, rather than polling forever. All of this is copied to \code{paste0(mod.nam, "_logs")} before \code{dump.path} is (optionally) deleted, whether the run succeeds or fails.

  \strong{Reproducibility.} With \code{seed}, a run reproduces the saved model object exactly, given the same \code{seed}, model file, \code{data}, \code{constants}, \code{inits}, \code{parameters}, \code{parameters2}, \code{SamplerSourcePath} and sampling settings (\code{nc}, \code{ni}, \code{nb}, \code{nt}, \code{nt2}, \code{check.freq}, \code{max.tries}, \code{max.samples.saved}, \code{Rht.required}, \code{neff.required} and the \code{par.*} arguments), on the same platform with the same R and NIMBLE versions. It is not guaranteed across machines or package versions. NIMBLE's compiled samplers, including \pkg{nimbleHMC}'s NUTS, draw from the worker's R random-number stream, so seeding the worker is sufficient.

  This also requires the number of blocks gathered to be the same every time, which the design ensures independently of timing. In automated mode, each convergence check waits until every chain has finished the same number of blocks. In manual mode (\code{check.freq = NULL}), workers keep sampling until the parent's periodic check (every 5 minutes) sees enough draws, so they typically overshoot; a seeded run therefore gathers exactly the smallest number of blocks that satisfies \code{max.samples.saved} and ignores the rest. An unseeded manual-mode run still gathers every block available when sampling stops. There is no way to resume an interrupted run, so there is no resume path that could break reproducibility.
}
\value{If assigned to object in global environment and rtrn.model = T, model output is returned. Else, model output is saved to file. Regardless, returned value consists of a list with four elements, plus up to four more when applicable:
  \item{mcmcOutput}{posterior samples as returned by mcmcOutput::mcmcOutput}
  \item{iter.key}{A \code{data.frame(chn, iter)} with one row per retained draw of \code{mcmcOutput}, aligned with \code{as.matrix(mcmcOutput)}'s row order (chain-major, ascending iteration within chain): which chain and absolute NIMBLE iteration each row came from. \code{mcmcOutput}'s own row metadata carries no such information on its own (coda's \code{mcpar} is always a bare \code{1:nrow}, unrelated to the true iteration numbers, especially once \code{max.samples.saved} forces an uneven subsample) - this is the only way to recover it. See \code{parameters2} above for joining this against \code{monitors2$iter.key}.}
  \item{summary}{data frame containing summaries of posterior estimates and convergence criteria (Rhat and n_effective) for all parameters in mcmcOutput}
  \item{mcmc.info}{A vector with elemets: nchains = number of MCMC sampling chains, niterations = realized pre-thinning chain length including burnin, burnin = realized pre-thinning length of each chain discarded as burn-in, nthin = realized thinning expressed as the ith sample retained. Note that nthin can be a non-whole number if max.samples.saved is specified, in which case it is the mean interval between retained samples.}
  \item{monitors2}{Present only when \code{parameters2} is non-empty and \code{consolidate.monitors2 = TRUE} (the default): the second monitor set, as \code{list(out = <[ndraw, nparam] matrix with chains stacked>, iter.key = <data.frame(chn, iter)>)} - the same content written to \code{paste0(mod.nam, "_monitors2.rds")}. See \code{parameters2} above for how \code{iter.key} relates rows here back to rows of \code{mcmcOutput} (via the top-level \code{iter.key} above).}
  \item{warning}{Present only in automated mode when convergence targets were not met within \code{max.tries} checks: a message saying so.}
  \item{rng}{Present only when \code{seed} is supplied: \code{list(seed, kind, streams)} - the seed as an integer; \code{kind}, the generator, normal and sample kinds the streams use (L'Ecuyer-CMRG, Inversion, Rejection); and \code{streams}, an \code{nc} x 7 integer matrix whose row \code{c} is chain \code{c}'s random-number state at the start of sampling.}
  \item{inits}{Present only when \code{inits} was a per-chain list or a function: a list of the \code{nc} named lists of initial values actually used, in chain order.}
  }
\author{
  Quresh S. Latif and Bryan L. Nuse, Bird Conservancy of the Rockies.

  Substantially revised with Claude (Anthropic) assistance across multiple sessions, including
  a full code review and bug-fix pass; see \code{NEWS.md} for what changed, including which
  changes affect results for existing callers.
}

\keyword{nimble}
\keyword{mcmcOutput}
\keyword{parallel processing}
