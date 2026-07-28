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
  if(!is.na(SamplerSourcePath)) require(nimbleHMC)
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
    if(!is.na(SamplerSourcePath)) source(SamplerSourcePath)
    nm.mcmc <- buildMCMC(nm.conf)
    nm.C <- compileNimble(nm, dirName = tmp.path)
    comp.mcmc <- compileNimble(nm.mcmc, project = nm.C, dirName = tmp.path)
    # cat(paste0("Calculate check: ", comp.mcmc$calculate(), "\n")) # Worth doing outside function but pointless here.
    comp.mcmc$run(niter = n.iter)
  } else {
    comp.mcmc$run(niter = n.iter, reset = FALSE, resetMV = TRUE)
  }
  samp <- as.matrix(comp.mcmc$mvSamples)
  if(length(monitors2) > 0) {
    samp2 <- as.matrix(comp.mcmc$mvSamples2)
    save(samp, samp2, file = dump.file.path)
  } else {
    save(samp, file = dump.file.path)
  }

  return(comp.mcmc)
}
