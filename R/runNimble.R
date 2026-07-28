runNimble <-
  function(model.path, inits, data, constants, parameters,
           max.samples.saved, par.ignore = c(), par.dontign = c(),
           par.fuzzy.track = c(), fuzzy.threshold = 0.05,
           nc = 2, ni = 2000, nb = 0.5, nt = 1, mod.nam = "mod",
           rtrn.model = F, sav.model = T, check.freq = NULL,
           Rht.required = 1.1, neff.required = 100,
           max.tries = 10, dump.path = "dump",
           SamplerSourcePath = NA, delete.blocks = TRUE,
           parameters2 = character(), nt2 = NULL) {
  # parameters2 / nt2: optional second set of monitored parameters saved at its
  # own thinning rate, for quantities too large to retain at nt - latent state
  # arrays, for example. They take no part in convergence checking; they are
  # gathered once at the end of the run and written to
  # paste0(mod.nam, "_monitors2.rds") as a plain [ndraw, nparam] matrix with
  # chains stacked.
  #
  # nt2 defaults to ni, i.e. one retained draw per block per chain. That suits
  # the block architecture: each block dump stays small, and draws accumulate to
  # roughly nblks * nc. Set it smaller for more draws at proportionally more
  # disk. Defaults to no second set, in which case nothing changes.
    if(length(parameters2) > 0 && is.null(nt2)) nt2 <- ni
    if(length(parameters2) == 0) nt2 <- 1
    if(!rtrn.model & !sav.model) stop("There is no way for runNimble to save output. Set either rtrn.model = TRUE or sav.model = TRUE.")
    if(length(parameters2) > 0 && any(parameters2 %in% parameters))
      stop("Parameters appear in both parameters and parameters2: ",
           paste(intersect(parameters2, parameters), collapse = ", "),
           ". The second monitor set is for quantities too large to retain at ",
           "thinning rate nt; monitoring them in both sets defeats that purpose.")
    if(nb < 1 & ((ni - (ni * nb)) / nt) <= 100) stop("Increase iterations (ni), reduce burn-in, or reduce thinning. Too few samples for calculating Rhat.")
    if(nb >= 1 & ((ni - nb) / nt) <= 100) stop("Increase iterations (ni), reduce burn-in, or reduce thinning. Too few samples for calculating Rhat.")
    automate.convergence.checks <- !is.null(check.freq)
    # In manual mode max.samples.saved is the sole stopping criterion
    # (do.gather.check below divides by it); NULL or non-positive is only
    # meaningful with automated convergence checking, where gatherNimble()
    # itself treats NULL as "no cap".
    if(!automate.convergence.checks &&
       (is.null(max.samples.saved) || !is.numeric(max.samples.saved) || max.samples.saved <= 0))
      stop("max.samples.saved must be a positive number when check.freq is NULL ",
           "(manual mode); it is the sole criterion used to decide when sampling stops.")
    if(automate.convergence.checks && !is.null(max.tries) && max.tries < 1)
      stop("max.tries must be at least 1 (or NULL for unlimited retries).")

    # nimble is required (not just imported) because compileNimble()'s C++
    # interop has been unreliable in the past when the package is only
    # available via NAMESPACE imports rather than genuinely attached; every
    # other dependency below is resolved through this package's own NAMESPACE
    # imports (see DESCRIPTION's Imports: and NAMESPACE's importFrom()) and
    # needs no require() here.
    require(nimble)
    if(!is.na(SamplerSourcePath)) { # Included this option for NUTS sampler, which requires HMC
      if(!requireNamespace("nimbleHMC", quietly = TRUE))
        stop("Package 'nimbleHMC' is required when SamplerSourcePath is supplied. ",
             "Install it with install.packages('nimbleHMC').")
      require(nimbleHMC) # attached so a user-supplied SamplerSourcePath script can call its functions unqualified
    }
    if(dir.exists(dump.path)) unlink(dump.path, recursive = TRUE)
    dir.create(dump.path)
    directive.file <- paste0(dump.path, "/runNimbleDirective.txt")
    # Coordination files (directive.file, block<c>Status.txt) are read and
    # written by multiple processes with no locking. write_state() writes to a
    # temp file and renames into place, so a reader can never see a truncated
    # write - the worker script uses the same pattern (finding F26).
    write_state <- function(text, path) {
      tmp <- paste0(path, ".tmp")
      writeLines(text, tmp)
      file.rename(tmp, path)
    }
    write_state("GO", directive.file)

    # --- Failure-safety net -------------------------------------------------
    # Preserves worker logs and stops any still-running workers no matter how
    # this function exits (normal return, stop(), or an interrupt). Without
    # this, an error thrown after workers are launched left them sampling
    # forever against a directive file nobody would ever update again.
    proc <- NULL
    copy_logs <- function() {
      log.dir <- paste0(mod.nam, "_logs")
      dir.create(log.dir, showWarnings = FALSE)
      keep <- c("worker_logs", "worker_joblog.txt", "Check_log.csv", "m.csv",
                list.files(dump.path, pattern = "chn.*_error\\.log$"))
      for(f in keep) {
        src <- file.path(dump.path, f)
        if(file.exists(src)) file.copy(src, log.dir, recursive = TRUE, overwrite = TRUE)
      }
    }
    on.exit({
      try(writeLines("STOP", directive.file), silent = TRUE)
      if(!is.null(proc) && proc$is_alive()) try(proc$kill_tree(), silent = TRUE)
      try(copy_logs(), silent = TRUE)
    }, add = TRUE)

    # A chain that dies mid-run (bad inits, compilation failure, OOM) writes
    # dump.path/chn<c>_error.log via the worker script's tryCatch wrapper
    # before it exits. Call this at the top of every polling loop so a dead
    # worker is detected and reported immediately instead of leaving the
    # parent to poll forever (findings F14, F25, F27).
    check_workers_alive <- function() {
      err.files <- list.files(dump.path, pattern = "chn.*_error\\.log$", full.names = TRUE)
      if(length(err.files) > 0)
        stop("One or more chains failed. See: ", paste(err.files, collapse = "; "))
      if(!proc$is_alive())
        stop("Worker process exited unexpectedly (parallel exit status ",
             proc$get_exit_status(), "). Check ", dump.path,
             "/worker_logs and ", dump.path, "/worker_joblog.txt for details.")
    }
    count_distinct_chains <- function() {
      fl <- list.files(dump.path)
      fl <- fl[str_detect(fl, "mod_chn")]
      if(length(fl) == 0) return(0L)
      length(unique(str_split(fl, "_", simplify = TRUE)[,2]))
    }

    if(automate.convergence.checks) {
      for(cn in 1:nc) {
        status.file <- paste0(dump.path, "/block",cn,"Status.txt")
        write_state("GO", status.file)
      }
      rm(cn, status.file)
      check.log.file <- paste0(dump.path, "/Check_log.csv")
      write.csv(data.frame(Model = character(), Check = numeric(),
                           Time = character(), Status = character()),
                check.log.file, row.names = FALSE)
    }
    save(list = c("model.path", "constants", "data", "inits", "parameters", "ni",
                  "nt", "dump.path", "SamplerSourcePath", "check.freq",
                  "automate.convergence.checks", "directive.file",
                  "parameters2", "nt2"),
         file = paste0(dump.path, "/NimbleObjects.RData"))
    #[Create R script for kicking off nimble run here]. Call it "ModRunScript.R"
    #___________________________________________________________________________#
    writeLines(text = c(
      "require(nimble)",
      "require(manageNimble)",
      "chn <- commandArgs(trailingOnly = TRUE)[[1]]",
      "path.NimbleWorkspace <- commandArgs(trailingOnly = TRUE)[[2]]",
      "load(path.NimbleWorkspace)",
      "error.log <- paste0(dump.path, '/chn', chn, '_error.log')",
      # read_state()/write_state() mirror the parent's own coordination-file
      # handling: write_state() writes to a temp file and renames into place, so
      # a read can never see a truncated write; read_state() retries a read that
      # is momentarily empty or mid-rewrite rather than treating it as a state
      # (finding F26 - the worker previously had no such protection at all).
      "read_state <- function(path) {",
        "x <- try(readLines(path), silent = TRUE)",
        "while(inherits(x, 'try-error') || length(x) != 1) {",
          "Sys.sleep(1)",
          "x <- try(readLines(path), silent = TRUE)",
        "}",
        "x",
      "}",
      "write_state <- function(text, path) {",
        "tmp <- paste0(path, '.tmp')",
        "writeLines(text, tmp)",
        "file.rename(tmp, path)",
      "}",
      "tryCatch({",
      "if(!dir.exists(paste0(dump.path, '/tmp', chn))) dir.create(paste0(dump.path, '/tmp', chn))",
      "source(model.path)",
      "i <- 1",
      "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod.comp <- runNimbleBlock(mod.lst = list(model, constants, data, inits, parameters),",
      "SamplerSourcePath = SamplerSourcePath,",
      "n.iter = ni, n.thin = nt, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path,",
      "monitors2 = parameters2, n.thin2 = nt2)",
      "if(automate.convergence.checks) {",
        "status.file <- paste0(dump.path, '/block',chn,'Status.txt')",
        "status.chain <- read_state(status.file)",
        "i.stop <- check.freq",
      "}",
      "directive <- read_state(directive.file)",
      "while(directive != 'STOP') {",
      "if(automate.convergence.checks) {",
        "if(directive == 'GO' & status.chain == 'GO' & i < i.stop) {",
          "i <- i + 1",
          "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
          "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path,",
          "monitors2 = parameters2)",
          "directive <- read_state(directive.file)",
        "} else if(directive == 'GO' & status.chain == 'GO' & i == i.stop) {",
          "write_state('STOP', status.file)",
          "status.chain <- read_state(status.file)",
        "} else if(directive == 'PAUSE' | directive == 'GO' & status.chain == 'STOP') {",
          "Sys.sleep(10)",
          "directive <- read_state(directive.file)",
          "status.chain <- read_state(status.file)",
        "} else if(status.chain == 'RESUME') {",
          "i.stop <- i.stop + check.freq",
          "write_state('GO', status.file)",
          "status.chain <- read_state(status.file)",
          "directive <- read_state(directive.file)",
        "} else {",
          "stop(paste0('Undefined condition reached during automated convergence process ',",
                      "'on chain ', chn, ': directive=', directive, ', status.chain=', status.chain, '.'))",
        "}",
      "} else {",
        "if(directive == 'GO') {",
          "i <- i + 1",
          "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
          "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path,",
          "monitors2 = parameters2)",
          "directive <- read_state(directive.file)",
        "} else if(directive == 'PAUSE') {",
          "Sys.sleep(10)",
          "directive <- read_state(directive.file)",
        "}",
      "}",
      "}",
      "}, error = function(e) {",
        "writeLines(c(paste0('Time: ', as.character(Sys.time())),",
                    "paste0('Chain: ', chn),",
                    "paste0('Error: ', conditionMessage(e)),",
                    "paste0('Call: ', paste(deparse(conditionCall(e)), collapse = ' '))),",
                  "error.log)",
        "quit(save = 'no', status = 1)",
      "})"
    ),
    con = paste0(dump.path, "/ModRunScript.R"))
    #___________________________________________________________________________#
    # --results/--joblog give a raw stdout/stderr/exit-status record per chain,
    # independent of the tryCatch above - catching crashes below the R level
    # (segfault, OOM kill) that tryCatch cannot.
    proc <- process$new(command = "parallel",
                        args = c("--results", paste0(dump.path, "/worker_logs"),
                                 "--joblog", paste0(dump.path, "/worker_joblog.txt"),
                                 "Rscript", paste0(dump.path, "/ModRunScript.R"),
                                 "{}",
                                 paste0(dump.path, "/NimbleObjects.RData"),
                                 ":::",
                                 1:nc))
    rm(data)
    gc(verbose = FALSE)
    mod.check.result <- FALSE
    run.complete <- FALSE
    nchecks <- 0
    while(count_distinct_chains() < nc) { # Wait until every chain has written at least one block.
      check_workers_alive()
      Sys.sleep(10)
    }
    if(automate.convergence.checks) {
      while(ifelse(is.null(max.tries), !mod.check.result,
               !mod.check.result & nchecks < max.tries)) {
        check_workers_alive()
        status.chains <- character(length = nc)
        for(cn in 1:nc) {
          try.status.check <- try(status.chains[cn] <- readLines(paste0(dump.path, '/block',cn,'Status.txt')), silent = TRUE)
          while(inherits(try.status.check, "try-error")) {
            Sys.sleep(1)
            try.status.check <- try(status.chains[cn] <- readLines(paste0(dump.path, '/block',cn,'Status.txt')), silent = TRUE)
          }
        }
        if(any(status.chains != "STOP")) {  # Wait for chains to finish sampling up to target amount.
          Sys.sleep (60)
          for(cn in 1:nc) {
            try.status.check <- try(status.chains[cn] <- readLines(paste0(dump.path, '/block',cn,'Status.txt')), silent = TRUE)
            while(inherits(try.status.check, "try-error")) {
              Sys.sleep(1)
              try.status.check <- try(status.chains[cn] <- readLines(paste0(dump.path, '/block',cn,'Status.txt')), silent = TRUE)
            }
          }
        } else {

          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
          if(nb > 1) {  # Need to keep sampling until we've passed burnin if burnin is absolute.
            while(!nrow(check.blocks$m) > 0) {
              check_workers_alive()
              for(cn in 1:nc) write_state("RESUME", paste0(dump.path, '/block',cn,'Status.txt'))
              write_state("GO", directive.file)
              Sys.sleep(60)
              check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
            }
          }
          if(length(unique(check.blocks$m[,1])) != nc) stop("Error: Expected ", nc,
              " chains but found ", length(unique(check.blocks$m[,1])), ". Code debug needed somewhere.")
          write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))

          nblks <- check.blocks$nblks
          nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
          ni.now <- ni*nblks

          write_state("PAUSE", directive.file)
          mod.out <- suppressWarnings(
            gatherNimble(read.path = dump.path, directive.file = directive.file,
                         burnin = nb, ni.block = ni, base.thin = nt,
                         max.samples.saved = max.samples.saved)
          )
          mod.check <- checkNimble(mod.out$out, directive.file = directive.file, Rht.required = Rht.required,
                                   neff.required = neff.required, par.ignore = par.ignore, par.dontign = par.dontign,
                                   par.fuzzy.track = par.fuzzy.track, fuzzy.threshold = fuzzy.threshold,
                                   spit.summary = TRUE, mod.nam = mod.nam)
          mod.check.result <- mod.check$result
          prp.fuzzy.not.converged <- mod.check$prp.fuzzy.not.converged
          thin.additional <- mod.out$additional.thin.rate
          nt.now <- ifelse(automate.convergence.checks, nt*thin.additional, nt)
          mcmc.info <- c(nchains = nc, niterations = ni.now,
                         burnin = nb.now, nthin = nt.now)
          sumTab <- mod.check$s
          mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info)
          if(sav.model) R.utils::saveObject(mod, mod.nam)
          if(rtrn.model) assign(mod.nam, mod, envir = .GlobalEnv)
          if(length(par.fuzzy.track) == 0) {
            status <- paste0("Max Rhat = ", max(sumTab$Rhat), " and min neff = ", min(sumTab$n.eff))
          } else {
            status <- paste0("Max Rhat = ", max(sumTab$Rhat), ", min neff = ", min(sumTab$n.eff),
                             ", and proportion fuzzy parameters not converged = ",
                             round(prp.fuzzy.not.converged, digits = 2))
          }
          check.log <- read.csv(check.log.file, colClasses = c("character",
                                                               "numeric",
                                                               "character",
                                                               "character"))
          check.log <- check.log %>% bind_rows(
            data.frame(Model = mod.nam, Check = nchecks + 1, Time = as.character(Sys.time()),
                       Status = status)
          )
          write.csv(check.log, check.log.file, row.names = FALSE)
          nchecks <- nchecks + 1

          if(ifelse(is.null(max.tries), !mod.check.result,
                    !mod.check.result & nchecks < max.tries)) {
            for(cn in 1:nc) write_state("RESUME", paste0(dump.path, '/block',cn,'Status.txt'))
            write_state("GO", directive.file)
            suppressWarnings(rm(mod, mod.check)) # mod.out - need mod.out below, so keep it.
            gc()
          }
        } # Close if(any(status.chains != "STOP")) {} else {}
      } # Close primary while loop, i.e., while(ifelse(is.null(max.tries), ..., ...))
      write_state("STOP", directive.file)
      if(!mod.check.result) {
        warn.message <- paste0("Rhat did not decrease after ", nchecks,
                               " checks. Model abandoned before reaching convergence targets.")
        mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info,
                    warning = warn.message)
        if(sav.model) R.utils::saveObject(mod, mod.nam)
        if(rtrn.model) assign(mod.nam, mod, envir = .GlobalEnv)
      }
      # Wait for workers to actually exit before gathering the second monitor
      # set, so a chain that hasn't yet noticed STOP can't add one more block
      # after the primary gatherNimble() snapshot above (finding F13).
      wait.start <- Sys.time()
      while(proc$is_alive() && as.numeric(Sys.time() - wait.start, units = "secs") < 300) Sys.sleep(2)
      if(proc$is_alive()) warning("Workers still running 5 minutes after STOP was issued; ",
                                  "gathering the second monitor set from current disk contents anyway.")
      # Gather the second monitor set, if any, before the block dumps are
      # removed. Done once at the end rather than at every convergence check:
      # these parameters take no part in convergence and are typically large.
      if(length(parameters2) > 0) {
        gatherNimble2(read.path = dump.path, burnin = nb, ni.block = ni, nt2 = nt2,
                      save.path = paste0(mod.nam, "_monitors2.rds"))
      }
      copy_logs()
      if(delete.blocks) unlink(dump.path, recursive = TRUE)
      gc(verbose = FALSE)
    } else {
      if(nb > 1) {  # Also need to wait until we've passed burnin if burnin is absolute.
        check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        while(!nrow(check.blocks$m) > 0) {
          check_workers_alive()
          Sys.sleep(10)
          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        }
      }
      nblks.previous <- 0 # Will be updated as we go.
      while(!run.complete) {
        check_workers_alive()
        check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        if(length(unique(check.blocks$m[,1])) != nc) stop("Error: Expected ", nc,
            " chains but found ", length(unique(check.blocks$m[,1])), ". Code debug needed somewhere.")
        while(length(unique(check.blocks$m[,1])) < nc) {
          check_workers_alive()
          Sys.sleep(60)
          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        }
        write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))
        nblks <- check.blocks$nblks
        nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
        ni.now <- ni*nblks

        do.gather.check <- (ni.now - nb.now) >= max.samples.saved * nt
        if(do.gather.check) {
          write_state("STOP", directive.file)
          mod.out <- suppressWarnings(
            gatherNimble(read.path = dump.path, directive.file = directive.file,
                         burnin = nb, ni.block = ni, base.thin = nt,
                         max.samples.saved = max.samples.saved)
          )
          if(!is.null(par.ignore)) mod.out.subset <- mcmcOutputSubset(mod.out$out, par.drop = par.ignore)
          mcmc.info <- c(nchains = nc, niterations = ni.now,
                         burnin = nb.now, nthin = nt * mod.out$additional.thin.rate)
          if(!is.null(par.ignore)) {
            sumTab <- summary(mod.out.subset, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
          } else {
            sumTab <- summary(mod.out$out, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
          }
          mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info)
          if(sav.model) R.utils::saveObject(mod, mod.nam)
          if(rtrn.model) assign(mod.nam, mod, envir = .GlobalEnv)
          run.complete <- TRUE
        } else {
          Sys.sleep(300)
        } # Close if(do.gather.check) loop
        nchecks <- nchecks + 1
      } # Close while(!run.complete)
      # Wait for workers to actually exit before gathering the second monitor
      # set (finding F13, see the automated branch above for the rationale).
      wait.start <- Sys.time()
      while(proc$is_alive() && as.numeric(Sys.time() - wait.start, units = "secs") < 300) Sys.sleep(2)
      if(proc$is_alive()) warning("Workers still running 5 minutes after STOP was issued; ",
                                  "gathering the second monitor set from current disk contents anyway.")
      # Gather the second monitor set, if any, before the block dumps are
      # removed. Done once at the end rather than at every convergence check:
      # these parameters take no part in convergence and are typically large.
      if(length(parameters2) > 0) {
        gatherNimble2(read.path = dump.path, burnin = nb, ni.block = ni, nt2 = nt2,
                      save.path = paste0(mod.nam, "_monitors2.rds"))
      }
      copy_logs()
      if(delete.blocks) unlink(dump.path, recursive = TRUE)
      gc(verbose = FALSE)
    }
  }
