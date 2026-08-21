runNimbleBlock <- function (mod.lst = NULL, comp.mcmc = NULL, n.iter = 1000,
                       n.thin = 1, tmp.path, dump.file.path,
                       SamplerSourcePath = NA,
                       monitors2 = character(), n.thin2 = 1) {
  # monitors2 / n.thin2: optional second set of monitored parameters, saved at
  # its own (typically much coarser) thinning rate. Intended for quantities too
  # large to retain at the primary thinning rate - latent state arrays, for
  # example - which are needed for post hoc calculations but not for convergence
  # diagnostics. Defaults to no second set, in which case behaviour and the
  # contents of the block dump are unchanged.
  #
  # On continuation blocks the MCMC is already configured, so monitors2 is used
  # here only to decide whether mvSamples2 should be written to the dump file.
  # It must still be supplied on those calls for the second set to be saved.
  require(nimble)
  if(!is.na(SamplerSourcePath)) {
    if(!requireNamespace("nimbleHMC", quietly = TRUE))
      stop("Package 'nimbleHMC' is required when SamplerSourcePath is supplied. ",
           "Install it with install.packages('nimbleHMC').")
    require(nimbleHMC)
  }
  stopifnot(sum(c(is.null(mod.lst), is.null(comp.mcmc))) == 1)

  if (!is.null(mod.lst)) {
    nm <- nimbleModel(code = mod.lst[[1]], constants = mod.lst[[2]], data = mod.lst[[3]], inits = mod.lst[[4]], calculate = FALSE)
    # cat(paste0("Initialization info:\n", nm$initializeInfo(), "\n")) # Worth doing outside function but pointless here.
    if(length(monitors2) > 0) {
      nm.conf <- configureMCMC(nm, monitors = mod.lst[[5]], thin = n.thin,
                               monitors2 = monitors2, thin2 = n.thin2)
    } else {
      nm.conf <- configureMCMC(nm, monitors = mod.lst[[5]], thin = n.thin)
    }
    # local = TRUE evaluates the sourced file in this call frame, where nm.conf
    # actually lives. Without it, source()'s default (local = FALSE) evaluates in
    # .GlobalEnv instead - a SamplerSourcePath file's nm.conf$removeSamplers()/
    # addSampler() calls would then either fail with "object 'nm.conf' not found"
    # (the normal case, since each worker is a fresh Rscript process with no
    # global nm.conf) or, if a same-named global object happened to exist, would
    # silently mutate that unrelated object while nm.conf here - the one actually
    # passed to buildMCMC() below - stayed untouched.
    if(!is.na(SamplerSourcePath)) source(SamplerSourcePath, local = TRUE)
    nm.mcmc <- buildMCMC(nm.conf)
    nm.C <- compileNimble(nm, dirName = tmp.path)
    comp.mcmc <- compileNimble(nm.mcmc, project = nm.C, dirName = tmp.path)
    # cat(paste0("Calculate check: ", comp.mcmc$calculate(), "\n")) # Worth doing outside function but pointless here.
    comp.mcmc$run(niter = n.iter)
  } else {
    comp.mcmc$run(niter = n.iter, reset = FALSE, resetMV = TRUE)
  }
  samp <- as.matrix(comp.mcmc$mvSamples)
  # Write to a temp file and rename into place. countNimbleBlocks()/gatherNimble()
  # list dump.file.path's final name to decide what is safe to read; an atomic
  # rename means a file visible under that name is always complete, removing the
  # need for the Linux-only lsof polling that used to guard against partial reads.
  tmp.path.file <- paste0(dump.file.path, ".tmp")
  if(length(monitors2) > 0) {
    samp2 <- as.matrix(comp.mcmc$mvSamples2)
    save(samp, samp2, file = tmp.path.file)
  } else {
    save(samp, file = tmp.path.file)
  }
  file.rename(tmp.path.file, dump.file.path)

  return(comp.mcmc)
}
