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

    require(nimble)
    if(!is.na(SamplerSourcePath)) require(nimbleHMC) # Included this option for NUTS sampler, which requires HMC
    require(processx)
    require(coda)
    require(mcmcOutput)
    if(dir.exists(dump.path)) unlink(dump.path, recursive = TRUE)
    dir.create(dump.path)
    directive.file <- paste0(dump.path, "/runNimbleDirective.txt")
    writeLines("GO", directive.file)
    if(automate.convergence.checks) {
      for(cn in 1:nc) {
        status.file <- paste0(dump.path, "/block",cn,"Status.txt")
        writeLines("GO", status.file)
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
      "if(!dir.exists(paste0(dump.path, '/tmp', chn))) dir.create(paste0(dump.path, '/tmp', chn))",
      "source(model.path)",
      "i <- 1",
      "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod.comp <- runNimbleBlock(mod.lst = list(model, constants, data, inits, parameters, SamplerSourcePath = SamplerSourcePath),",
      "n.iter = ni, n.thin = nt, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path,",
      "monitors2 = parameters2, n.thin2 = nt2)",
      "if(automate.convergence.checks) {",
        "status.file <- paste0(dump.path, '/block',chn,'Status.txt')",
        "status.chain <- readLines(status.file)",
        "i.stop <- check.freq",
      "}",
      "directive <- readLines(directive.file)",
      "while(directive != 'STOP') {",
      "if(automate.convergence.checks) {",
        "if(directive == 'GO' & status.chain == 'GO' & i < i.stop) {",
          "i <- i + 1",
          "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
          "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path,",
          "monitors2 = parameters2)",
          "directive <- readLines(directive.file)",
        "} else if(directive == 'GO' & status.chain == 'GO' & i == i.stop) {",
          "writeLines('STOP', status.file)",
          "status.chain <- readLines(status.file)",
        "} else if(directive == 'PAUSE' | directive == 'GO' & status.chain == 'STOP') {",
          "Sys.sleep(10)",
          "directive <- readLines(directive.file)",
          "status.chain <- readLines(status.file)",
        "} else if(status.chain == 'RESUME') {",
          "i.stop <- i.stop + check.freq",
          "writeLines('GO', status.file)",
          "status.chain <- readLines(status.file)",
          "directive <- readLines(directive.file)",
        "} else {",
          "error.message <- paste0('Undefined condition reached during automated convergence process on chain ', chn, '.')",
          "writeLines(error.message, paste0(dump.path, '/block',cn,'Error.txt'))",
        "}",
      "} else {",
        "if(directive == 'GO') {",
          "i <- i + 1",
          "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
          "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path,",
          "monitors2 = parameters2)",
          "directive <- readLines(directive.file)",
        "} else if(directive == 'PAUSE') {",
          "Sys.sleep(10)",
          "directive <- readLines(directive.file)",
        "}",
      "}",
      "}"
    ),
    con = paste0(dump.path, "/ModRunScript.R"))
    #___________________________________________________________________________#
    proc <<- process$new(command = "parallel",
                        args = c("Rscript", eval(paste0(dump.path, "/ModRunScript.R")),
                                 "{}",
                                 eval(paste0(dump.path, "/NimbleObjects.RData")),
                                 ":::",
                                 1:nc))
    proc
    rm(data)
    gc(verbose = FALSE)
    mod.check.result <- FALSE
    run.complete <- FALSE
    nchecks <- 0
    while(sum(str_detect(list.files(dump.path), "mod_chn")) < nc) {Sys.sleep(10)} # Wait until proc has written at least one file for each chain before going on.
    if(automate.convergence.checks) {
      while(ifelse(is.null(max.tries), !mod.check.result,
               !mod.check.result & nchecks < max.tries)) {
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
              for(cn in 1:nc) writeLines("RESUME", paste0(dump.path, '/block',cn,'Status.txt'))
              writeLines("GO", directive.file)
              Sys.sleep(60)
              check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
            }
          }
          if(length(unique(check.blocks$m[,1])) > nc) stop("Error: Too many sampling blocks created. Code debug needed somewhere.")
          write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))
          
          nblks <- check.blocks$nblks
          nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
          ni.now <- ni*nblks
          
          writeLines("PAUSE", directive.file)
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
          thin.additional <- mod.out$additional.thin.rate
          nt.now <- ifelse(automate.convergence.checks, nt*thin.additional, nt)
          mcmc.info <- c(nchains = nc, niterations = ni.now,
                         burnin = nb.now, nthin = nt.now)
          sumTab <- sumTab.focal <- mod.check$s
          # Exclude par.fuzzy.track as well as par.ignore. checkNimble() now
          # governs fuzzy-tracked parameters by fuzzy.threshold and permits some
          # of them to be unconverged or to have an uncomputable Rhat; leaving
          # them in sumTab.focal would let the "not being sampled" check below
          # halt the run for exactly the parameters the fuzzy tolerance exists to
          # accommodate.
          if(length(c(par.ignore, par.fuzzy.track)) > 0) {
            sumTab.focal <- sumTab %>%
              filter(!(str_split(Parameter, "\\[", simplify = TRUE)[,1]) %in%
                       c(par.ignore, par.fuzzy.track))
          }
          if(any(is.na(sumTab.focal$Rhat)) | any(is.na(sumTab.focal$n.eff))) {
            # proc$kill_tree()
            writeLines("STOP", directive.file)
            write.csv(sumTab.focal, paste0("Model_summary_PID",proc$get_pid(),".csv"))
            stop(paste0("Error: One or more parameters is not being sampled.",
                        " Check data, initial values, etc., and try again.",
                        " See 'Model_summary_PID",proc$get_pid(),
                        ".csv' for parameters missing Rhat or n.eff."))
          }
          if(length(par.fuzzy.track) > 0) {
            Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
            if(!any(names(sumTab) == "Rhat")) {
              # proc$kill_tree()
              writeLines("STOP", directive.file)
              stop("Stopped model run because Rhat not calculated.")
            }
            # Recomputed here only for the status log. Kept identical to the
            # calculation in checkNimble(), including NAs in the numerator, so
            # the logged proportion matches the one that decided convergence.
            Rht.fuzzy <- sumTab %>%
              filter(str_split(Parameter, "\\[", simplify = TRUE)[,1] %in% par.fuzzy.track) %>%
              pull(Rhat)
            prp.fuzzy.not.coverged <-
              (sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = TRUE) +
                 sum(is.na(Rht.fuzzy))) / length(Rht.fuzzy)
          }
          mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info)
          if(sav.model) R.utils::saveObject(mod, mod.nam)
          if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
          if(!mod.check.result & length(par.fuzzy.track) == 0) {
            status <- paste0("Max Rhat = ", max(sumTab.focal$Rhat), " and min neff = ",
                             min(sumTab.focal$n.eff))
          } else if(!mod.check.result & length(par.fuzzy.track) > 0) {
            status <- paste0("Max Rhat = ", max(sumTab.focal$Rhat), ", min neff = ",
                             min(sumTab.focal$n.eff),
                             ", and proportion fuzzy parameters not converged = ",
                             round(prp.fuzzy.not.coverged, digits = 2))
          } else if(mod.check.result & length(par.fuzzy.track) == 0) {
            status <- paste0("Max Rhat = ", max(sumTab.focal$Rhat),
                             " and min neff = ", min(sumTab.focal$n.eff))
          } else {
            status <- paste0("Max Rhat = ", max(sumTab.focal$Rhat),
                             ", min neff = ", min(sumTab.focal$n.eff),
                             ", and proportion fuzzy parameters not converged = ",
                             round(prp.fuzzy.not.coverged, digits = 2))
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
            for(cn in 1:nc) writeLines("RESUME", paste0(dump.path, '/block',cn,'Status.txt'))
            writeLines("GO", directive.file)
            suppressWarnings(rm(mod, mod.check, sumTab, sumTab.focal)) # mod.out - need mod.out below, so keep it.
            gc()
          }
        } # Close if(any(status.chains != "STOP")) {} else {}
      } # Close primary while loop, i.e., while(ifelse(is.null(max.tries), ..., ...))
      writeLines("STOP", directive.file)
      if(!mod.check.result) {
        warn.message <- paste0("Rhat did not decrease after ", nchecks,
                               " checks. Model abandoned before reaching convergence targets.")
        mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info,
                    warning = warn.message)
        if(sav.model) R.utils::saveObject(mod, mod.nam)
        if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
      }
      # Gather the second monitor set, if any, before the block dumps are
      # removed. Done once at the end rather than at every convergence check:
      # these parameters take no part in convergence and are typically large.
      if(length(parameters2) > 0) {
        gatherNimble2(read.path = dump.path, burnin = nb, ni.block = ni,
                      save.path = paste0(mod.nam, "_monitors2.rds"))
      }
      if(delete.blocks) unlink(dump.path, recursive = TRUE)
      gc(verbose = FALSE)
    } else {
      if(nb > 1) {  # Also need to wait until we've passed burnin if burnin is absolute.
        check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        while(!nrow(check.blocks$m) > 0) {
          Sys.sleep(10)
          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        }
      }
      nblks.previous <- 0 # Will be updated as we go.
      while(!run.complete) {
        check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        if(length(unique(check.blocks$m[,1])) > nc) stop("Error: Too many sampling blocks created. Code debug needed somewhere.")
        while(length(unique(check.blocks$m[,1])) < nc) {
          Sys.sleep(60)
          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        }
        write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))
        nblks <- check.blocks$nblks
        nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
        ni.now <- ni*nblks

        do.gather.check <- (ni.now - nb.now) >= max.samples.saved * nt
        if(do.gather.check) {
          writeLines("STOP", directive.file)
          mod.out <- suppressWarnings(
            gatherNimble(read.path = dump.path, directive.file = directive.file,
                         burnin = nb, ni.block = ni, base.thin = nt,
                         max.samples.saved = max.samples.saved)
          )
          if(!is.null(par.ignore)) mod.out.subset <- mcmcOutputSubset(mod.out$out, par.drop = par.ignore)
          mcmc.info <- c(nchains = nc, niterations = ni.now,
                         burnin = nb.now, nthin = nt)
          if(!is.null(par.ignore)) {
            sumTab <- summary(mod.out.subset, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
          } else {
            sumTab <- summary(mod.out$out, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
          }
          mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info)
          if(sav.model) R.utils::saveObject(mod, mod.nam)
          if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
          run.complete <- TRUE
        } else {
          Sys.sleep(300)
        } # Close if(do.gather.check) loop
        nchecks <- nchecks + 1
      } # Close while(!run.complete)
      # Gather the second monitor set, if any, before the block dumps are
      # removed. Done once at the end rather than at every convergence check:
      # these parameters take no part in convergence and are typically large.
      if(length(parameters2) > 0) {
        gatherNimble2(read.path = dump.path, burnin = nb, ni.block = ni,
                      save.path = paste0(mod.nam, "_monitors2.rds"))
      }
      if(delete.blocks) unlink(dump.path, recursive = TRUE)
      gc(verbose = FALSE)
    }
  }
