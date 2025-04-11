runNimble <-
  function(model.path, inits, data, constants, parameters,
           max.samples.saved, par.ignore = c(), par.dontign = c(),
           par.fuzzy.track = c(), fuzzy.threshold = 0.05,
           nc = 2, ni = 2000, nb = 0.5, nt = 1, mod.nam = "mod",
           rtrn.model = F, sav.model = T, check.freq = NULL,
           Rht.required = 1.1, neff.required = 100,
           max.tries = 10, dump.path = "dump",
           SamplerSourcePath = NA, delete.blocks = TRUE) {
    if(!rtrn.model & !sav.model) stop("There is no way for runNimble to save output. Set either rtrn.model = TRUE or sav.model = TRUE.")
    if(nb < 1 & ((ni - (ni * nb)) / nt) <= 100) stop("Increase iterations (ni), reduce burn-in, or reduce thinning. Too few samples for calculating Rhat.")
    if(nb >= 1 & ((ni - nb) / nt) <= 100) stop("Increase iterations (ni), reduce burn-in, or reduce thinning. Too few samples for calculating Rhat.")
    automate.convergence.checks <- !is.null(check.freq)

    require(nimble)
    if(!is.na(SamplerSourcePath)) require(nimbleHMC) # Included this option for NUTS sampler, which requires HMC
    require(processx)
    require(coda)
    require(mcmcOutput)
    if(!dir.exists(dump.path)) dir.create(dump.path)
    directive.file <- paste0(dump.path, "/runNimbleDirective.txt")
    writeLines("GO", directive.file)
    if(automate.convergence.checks) {
      for(cn in 1:nc) {
        status.file <- paste0(dump.path, "/block",cn,"Status.txt")
        writeLines("GO", status.file)
      }
      check.log.file <- paste0(dump.path, "/Check_log.csv")
      write.csv(data.frame(Model = character(), Check = numeric(),
                           Time = character(), Status = character()),
                check.log.file, row.names = FALSE)
    }
    save(list = c("model.path", "constants", "data", "inits", "parameters", "ni",
                  "nt", "dump.path", "SamplerSourcePath", "check.freq",
                  "automate.convergence.checks", "directive.file"),
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
      "i <- 0",
      "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod.comp <- runNimbleBlock(mod.lst = list(model, constants, data, inits, parameters, SamplerSourcePath = SamplerSourcePath),",
      "n.iter = ni, n.thin = nt, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path)",
      "if(automate.convergence.checks) {",
        "status.file <- paste0(dump.path, '/block',cn,'Status.txt')",
        "status.chain <- readLines(status.file)",
        "i.stop <- check.freq",
      "}",
      "directive <- readLines(directive.file)",
      "while(directive != 'STOP') {",
      "if(automate.convergence.checks) {",
        "if(directive == 'GO' & status.chain == 'GO' & i < i.stop) {",
          "i <- i + 1",
          "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
          "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path)",
          "directive <- readLines(directive.file)",
        "} else if(directive == 'GO' & status.chain == 'GO' & i == i.stop) {",
          "writeLines('STOP', status.file)",
          "status.chain <- readLines(status.file)",
        "} else if(directive == 'PAUSE' | directive == 'GO' & status.chain == 'STOP') {",
          "Sys.sleep(10)",
          "directive <- readLines(directive.file)",
          "status.chain <- readLines(status.file)",
        "} else if(status.chain == 'RESUME') {",
          "i.stop <- i + check.freq",
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
          "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path)",
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
    nchecks <- 1
    while(sum(str_detect(list.files(dump.path), "mod_chn")) < nc) {Sys.sleep(10)} # Wait until proc has written at least one file for each chain before going on.
    if(automate.convergence.checks) {
      while(ifelse(is.null(max.tries), !mod.check.result,
               !mod.check.result & nchecks < max.tries)) {
        status.chains <- character(length = nc)
        for(cn in 1:nc) status.chains[cn] <- readLines(paste0(dump.path, '/block',cn,'Status.txt'))
        if(any(status.chains != "STOP")) {  # Wait for chains to finish sampling up to target amount.
          Sys.sleep (60)
          for(cn in 1:nc) status.chains[cn] <- readLines(paste0(dump.path, '/block',cn,'Status.txt'))
        } else {

          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
          if(length(unique(check.blocks$m[,1])) > nc) stop("Error: Too many sampling blocks created. Code debug needed somewhere.")
          write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))
          
          nblks <- check.blocks$nblks
          nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
          ni.now <- ni*nblks
          
          writeLines("PAUSE", paste0(dump.path, "/runNimbleDirective.txt"))
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
          if(length(par.ignore) > 0) {
            sumTab.focal <- sumTab %>%
              filter(!str_detect_any(Parameter, par.ignore))
            if(length(par.dontign) > 0) {
              sumTab.focal <- sumTab.focal %>% bind_rows(
                sumTab %>% filter(str_detect_any(Parameter, par.dontign))
              )
            }
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
            for(p in 1:length(par.fuzzy.track)) {
              pfuz <- par.fuzzy.track[p]
              Rht.fuzzy <- c(Rht.fuzzy,
                             sumTab %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                               pull(Rhat))
            }
            Rht.fuzzy <- Rht.fuzzy[-1]
            prp.fuzzy.not.coverged <- (sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = T) /
                                         length(Rht.fuzzy))
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
            data.frame(Model = mod.nam, Check = nchecks, Time = as.character(Sys.time()),
                       Status = status)
          )
          write.csv(check.log, check.log.file, row.names = FALSE)
          nchecks <- nchecks + 1
          
          if(ifelse(is.null(max.tries), !mod.check.result,
                    !mod.check.result & nchecks < max.tries)) {
            for(cn in 1:nc) writeLines("RESUME", paste0(dump.path, '/block',cn,'Status.txt'))
            writeLines("GO", directive.file)
            suppressWarnings(rm(mod, mod.out, mod.check, sumTab, sumTab.focal))
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
      if(delete.blocks) unlink(dump.path, recursive = TRUE)
      gc(verbose = FALSE)
    }
  }
