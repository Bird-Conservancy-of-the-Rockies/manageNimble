regex.escape <- function(x) {
  # Not exported. Escapes every extended-regex metacharacter in `x` so it can
  # be used as a pgrep/grep -f pattern and match only itself, literally. `-`
  # is deliberately not in this set: it is a metacharacter only inside `[...]`
  # bracket expressions, which this function never constructs.
  gsub("([.^$*+?()\\[\\]{}|\\\\])", "\\\\\\1", x, perl = TRUE)
}

liveNimbleWorkers <- function(dump.path) {
  # Not exported. Reports processes whose command line references dump.path -
  # this covers both GNU parallel's supervisor process and the individual
  # Rscript ModRunScript.R chain workers runNimble() launches, since both are
  # invoked with dump.path baked into their arguments (see runNimble()'s own
  # process$new() call). A single pattern is enough; forceFinishNimble() has
  # no knowledge of, and does not need, whatever outer process a caller's own
  # pipeline may have wrapped around runNimble() (e.g. a per-species driver
  # script) - that is the caller's concern, not this package's.
  if(nchar(Sys.which("pgrep")) == 0) {
    warning("pgrep not found on PATH - cannot verify whether a run is still live for '",
            dump.path, "'. Proceeding without this safety check.")
    return(character(0))
  }
  # pgrep -f matches its pattern as an extended regex, not a literal
  # substring. An unescaped dump.path containing "(", "+", "[", or other
  # metacharacters - all plausible, unexceptional directory names - can then
  # fail to match a genuinely live process at all, silently defeating this
  # safety check exactly when it matters most. Escape it first.
  hits <- suppressWarnings(system2("pgrep", args = c("-af", regex.escape(dump.path)), stdout = TRUE))
  # system2() shells out via `sh -c '...'`, and that transient shell's own command line
  # contains the search pattern just passed to it - pgrep then reports that shell process
  # as a match against itself (a guaranteed false positive: it has already exited by the
  # time anyone would go kill it). Filter out any hit that is just that self-referential
  # pgrep invocation.
  unique(hits[!str_detect(hits, "pgrep")])
}

forceFinishNimble <- function(dump.path, burnin, mod.nam,
                              par.ignore = c(), par.dontign = c(),
                              par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                              Rht.required = 1.1, neff.required = 100,
                              max.samples.saved = NULL,
                              commit = FALSE, sav.model = TRUE,
                              delete.blocks = FALSE, plot.result = TRUE,
                              check.live = TRUE, quiet = FALSE) {
  if(!dir.exists(dump.path)) stop("dump.path '", dump.path, "' does not exist.")
  nimble.objects.path <- file.path(dump.path, "NimbleObjects.RData")
  if(!file.exists(nimble.objects.path))
    stop("'", nimble.objects.path, "' not found - '", dump.path,
         "' does not look like a runNimble() dump directory.")

  # ni/nt are recovered from the dump directory's own NimbleObjects.RData
  # (written by runNimble() at the start of the run) rather than asked of the
  # caller, so there is no way to accidentally interpret the block files
  # under an inconsistent ni/nt.
  e <- new.env()
  load(nimble.objects.path, envir = e)
  if(!all(c("ni", "nt") %in% ls(e)))
    stop("'", nimble.objects.path, "' does not contain 'ni'/'nt' - unexpected dump format.")
  ni <- e$ni; nt <- e$nt

  live <- if(check.live) liveNimbleWorkers(dump.path) else character(0)
  if(!quiet) {
    if(length(live) > 0) {
      message("NOTE: process(es) referencing '", dump.path, "' are still running (this is fine ",
              "for a preview - reading completed blocks only):\n", paste(live, collapse = "\n"))
    } else if(check.live) {
      message("NOTE: no live process found referencing '", dump.path, "' - OK to preview or commit.")
    }
    message("ni (iterations per block) = ", ni, " | nt (thin) = ", nt)
  }

  ## Quick progress summary before gathering - how many blocks/iterations are on disk right
  ## now per chain (useful context when re-previewing against a live, growing run). ##
  blk.count <- countNimbleBlocks(dump.path, burnin = burnin, ni.block = ni)
  if(!quiet)
    message("Blocks on disk per chain (post-burnin, chains trimmed to shortest): ", blk.count$nblks,
            " | total iterations/chain so far: ", blk.count$nblks * ni,
            " | burn-in requested: ", burnin, ifelse(burnin < 1, " (fraction)", " (absolute iterations)"))
  if(nrow(blk.count$m) == 0)
    stop("Requested burnin (", burnin, ") exceeds the iterations accumulated so far for at least ",
         "one chain - reduce burnin or wait for more blocks.")

  ## Gather whatever blocks are on disk so far, applying the burn-in specified above. ##
  mod.out <- gatherNimble(read.path = dump.path, burnin = burnin, ni.block = ni, base.thin = nt,
                          max.samples.saved = max.samples.saved)

  chk <- checkNimble(mod.out$out, Rht.required = Rht.required, neff.required = neff.required,
                     par.ignore = par.ignore, par.dontign = par.dontign,
                     par.fuzzy.track = par.fuzzy.track, fuzzy.threshold = fuzzy.threshold,
                     spit.summary = TRUE, mod.nam = mod.nam)

  if(!quiet) {
    message("---- Convergence with burnin = ", burnin, " ----")
    message("Converged (Rhat <= ", Rht.required, " and neff >= ", neff.required, "): ", chk$result)
    # chk$s already excludes par.ignore'd parameters - checkNimble() drops them from
    # mcmcOutput itself (via mcmcOutputSubset()) before ever computing the summary - so
    # re-filtering here would be redundant. Reported the same way, and over the same set
    # of parameters, as runNimble()'s own "Max Rhat"/"min neff" status message.
    if(length(par.fuzzy.track) == 0) {
      message("Max Rhat: ", max(chk$s$Rhat), " | Min neff: ", min(chk$s$n.eff))
    } else {
      message("Max Rhat: ", max(chk$s$Rhat), " | Min neff: ", min(chk$s$n.eff),
              " | Proportion fuzzy parameters not converged: ",
              round(chk$prp.fuzzy.not.converged, digits = 2))
    }
  }
  if(plot.result) plot(mod.out$out)

  if(!commit) {
    if(!quiet)
      message("commit = FALSE - nothing saved. Adjust 'burnin' and call again to preview a ",
              "different value, or call again with commit = TRUE once this looks acceptable.")
    return(invisible(list(result = chk$result, summary = chk$s, mcmcOutput = mod.out$out,
                          nblks = mod.out$nblks, additional.thin.rate = mod.out$additional.thin.rate)))
  }

  ## Hard gate: committing is not allowed while the run is still live. Re-checked here (not just
  ## at the top of the function) in case the run was still live when this call started and has
  ## only been killed since - the caller may be re-running with commit = TRUE straight away. ##
  if(check.live) {
    live <- liveNimbleWorkers(dump.path)
    if(length(live) > 0)
      stop("Refusing to commit - process(es) referencing '", dump.path, "' are still running:\n",
           paste(live, collapse = "\n"),
           "\n\nKill the run first, then call forceFinishNimble() again with commit = TRUE. ",
           "Preview (commit = FALSE) does not require this. check.live = FALSE skips this ",
           "safety check entirely (not recommended).")
  }

  ## Build the same list structure runNimble() itself saves/returns. ##
  nc <- attr(mod.out$out, "nChains")
  nblks <- mod.out$nblks
  nb.now <- ifelse(burnin < 1, burnin * ni * nblks, burnin)
  mcmc.info <- c(nchains = nc, niterations = ni * nblks, burnin = nb.now,
                 nthin = nt * mod.out$additional.thin.rate)
  mod <- list(mcmcOutput = mod.out$out, summary = chk$s, mcmc.info = mcmc.info)
  if(!chk$result) {
    mod$warning <- paste0("Manually forced to finish with burnin = ", burnin,
                          " via forceFinishNimble() after slow convergence. Forced on ",
                          as.character(Sys.time()), ".")
  }

  if(sav.model) R.utils::saveObject(mod, mod.nam)
  if(delete.blocks) unlink(dump.path, recursive = TRUE)

  mod
}
