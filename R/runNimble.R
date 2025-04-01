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
      check.log.file <- paste0(dump.path, "/Check_log.csv")
      write.csv(data.frame(Model = character(), Check = numeric(), Time = character(),
                           Min_waited = numeric(), Status = character()),
                check.log.file, row.names = FALSE)
    }
    save(list = c("model.path", "constants", "data", "inits", "parameters", "ni",
                  "nt", "dump.path", "SamplerSourcePath"),
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
      "n.iter = ni, n.thin = nt, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path)",
      "GO <- readLines(paste0(dump.path, '/runNimbleDirective.txt')) == 'GO'",
      "while(GO) {",
      "i <- i + 1",
      "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod.comp <- runNimbleBlock(comp.mcmc = mod.comp, n.iter = ni, tmp.path = paste0(dump.path, '/tmp', chn), dump.file.path = dump.file.path)",
      "GO <- readLines(paste0(dump.path, '/runNimbleDirective.txt')) == 'GO'",
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
    if(nb > 1) {  # Also need to wait until we've passed burnin if burnin is absolute.
      check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
      while(!nrow(check.blocks$m) > 0) {
        Sys.sleep(10)
        check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
      }
    }
    nblks.previous <- 0 # Will be updated as we go.
    while(if(automate.convergence.checks) {
        ifelse(is.null(max.tries), !mod.check.result,
               !mod.check.result & nchecks < max.tries)
      } else {
        !run.complete
      }) {
      if(automate.convergence.checks) {
        Sys.sleep(check.freq)
        status <- "No new convergence information."
      }
      check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
      if(length(unique(check.blocks$m[,1])) > nc) stop("Error: Too many sampling blocks created. Code debug needed somewhere.")
      if(automate.convergence.checks) {
        while(length(unique(check.blocks$m[,1])) < nc |
              (check.blocks$nblks * ni) <= 100) {
          Sys.sleep(check.freq)
          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        }
      } else {
        while(length(unique(check.blocks$m[,1])) < nc) {
          Sys.sleep(60)
          check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
        }
      }
      write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))
      if(check.blocks$nblks > nblks.previous) {
        nblks.previous <- nblks <- check.blocks$nblks
        nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
        nt.now <- ifelse(automate.convergence.checks, nt*thin.additional, nt)
        ni.now <- ni*nblks
        
        do.gather.check <- automate.convergence.checks | (!automate.convergence.checks &
             ni.now >= max.samples.saved * nt)
        if(do.gather.check) {
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
            writeLines("STOP", paste0(dump.path, "runNimbleDirective.txt"))
            write.csv(sumTab.focal, paste0("Model_summary_PID",proc$get_pid(),".csv"))
            stop(paste0("Error: One or more parameters is not being sampled.",
                        " Check data, initial values, etc., and try again.",
                        " See 'Model_summary_PID",proc$get_pid(),
                        ".csv' for parameters missing Rhat or n.eff."))
          }
          if(automate.convergence.checks & length(par.fuzzy.track) > 0) {
            Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
            if(!any(names(sumTab) == "Rhat")) {
              # proc$kill_tree()
              writeLines("STOP", paste0(dump.path, "runNimbleDirective.txt"))
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
          if(!automate.convergence.checks) run.complete <- TRUE
        } # Close if(do.gather.check) loop
        
        if(automate.convergence.checks) {
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
        }
      }
      if(automate.convergence.checks) {
        check.log <- read.csv(check.log.file, colClasses = c("character", "numeric", "character",
                                                             "numeric", "character"))
        check.log <- check.log %>% bind_rows(
          data.frame(Model = mod.nam, Check = nchecks, Time = as.character(Sys.time()),
                     Min_waited = (nchecks * check.freq / 60), Status = status)
        )
        write.csv(check.log, check.log.file, row.names = FALSE)
      }
      nchecks <- nchecks + 1
    }
    # proc$kill_tree()
    writeLines("STOP", directive.file)
    if(automate.convergence.checks & !mod.check.result) {
      warn.message <- paste0("Rhat did not decrease after ", nchecks,
                            " checks. Model abandoned before reaching convergence targets.")
      mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info,
                  warning = warn.message)
      if(sav.model) R.utils::saveObject(mod, mod.nam)
      if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
    } # Close primary while loop, i.e., while(if(automate.convergence.checks) {...} else {...})
    if(delete.blocks) unlink(dump.path, recursive = TRUE)
    gc(verbose = FALSE)
  }
