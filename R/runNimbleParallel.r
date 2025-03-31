runNimbleParallel <-
  function(model.path, inits, data, constants, parameters,
           par.ignore = c(), par.dontign = c(),
           par.fuzzy.track = c(), fuzzy.threshold = 0.05,
           nc = 2, ni = 2000, nb = 0.5, nt = 1, mod.nam = "mod",
           max.samples.saved = 10000, rtrn.model = F, sav.model = T,
           Rht.required = 1.1, neff.required = 100,
           check.freq = 10, max.tries = NULL, dump.path = "dump",
           SamplerSourcePath = NA) {
    if(!rtrn.model & !sav.model) stop("There is no way for RunNimbleParallel to save output. Set either rtrn.model = TRUE or sav.model = TRUE.")
    if(nb < 1 & ((ni - (ni * nb)) / nt) <= 100) stop("Increase iterations (ni), reduce burn-in, or reduce thinning. Too few samples for calculating Rhat.")
    if(nb >= 1 & ((ni - nb) / nt) <= 100) stop("Increase iterations (ni), reduce burn-in, or reduce thinning. Too few samples for calculating Rhat.")

    require(nimble)
    if(!is.na(SamplerSourcePath)) require(nimbleHMC) # Included this option for NUTS sampler, which required HMC
    require(processx)
    require(coda)
    require(mcmcOutput)
    if(!dir.exists(dump.path)) dir.create(dump.path)
    check.log.file <- paste0(dump.path, "/Check_log.csv")
    write.csv(data.frame(Model = character(), Check = numeric(), Time = character(),
                         Min_waited = numeric(), Status = character()),
              check.log.file, row.names = FALSE)
    save(list = c("model.path", "constants", "data", "inits", "parameters", "ni", "nt", "SamplerSourcePath"),
         file = paste0(dump.path, "/NimbleObjects.RData"))
    #[Create R script for kicking off nimble run here]. Call it "ModRunScript.R"
    #___________________________________________________________________________#
    writeLines(text = c(
      "require(nimble)",
      "require(FunctionsBCR)",

      "chn <- commandArgs(trailingOnly = TRUE)[[1]]",
      "dump.path <- commandArgs(trailingOnly = TRUE)[[2]]",
      "path.NimbleWorkspace <- commandArgs(trailingOnly = TRUE)[[3]]",
      
      "load(path.NimbleWorkspace)",
      "source(model.path)",
      "i <- 1",
      "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod.comp <- runNimble(mod.lst = list(model, constants, data, inits, parameters, SamplerSourcePath = SamplerSourcePath),",
      "n.iter = ni, n.thin = nt, dump.file.path = dump.file.path)",
      "repeat{",
      "i <- i + 1",
      "dump.file.path <- paste0(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod.comp <- runNimble(comp.mcmc = mod.comp, n.iter = ni, dump.file.path = dump.file.path)",
      "}"
    ),
    con = paste0(dump.path, "/ModRunScript.R"))
    #___________________________________________________________________________#
    proc <<- process$new(command = "parallel",
                        args = c("Rscript", eval(paste0(dump.path, "/ModRunScript.R")),
                                 "{}",
                                 eval(dump.path),
                                 eval(paste0(dump.path, "/NimbleObjects.RData")),
                                 ":::",
                                 1:nc))
    proc
    rm(data)
    gc(verbose = FALSE)
    mod.check.result <- FALSE
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
    while(ifelse(is.null(max.tries), !mod.check.result, !mod.check.result & nchecks < max.tries)) {
      Sys.sleep(check.freq)
      status <- "No new convergence information."
      check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
      if(length(unique(check.blocks$m[,1])) > nc) stop("Error: Too many sampling blocks created. Code debug needed somewhere.")
      while(length(unique(check.blocks$m[,1])) < nc |
            (check.blocks$nblks * ni) <= 100) {
        Sys.sleep(check.freq)
        check.blocks <- countNimbleBlocks(read.path = dump.path, burnin = nb, ni.block = ni)
      }
      write.csv(check.blocks$m, paste0(dump.path, "/m.csv"))
      if(check.blocks$nblks > nblks.previous) {
        nblks.previous <- check.blocks$nblks
        
        mod.out <- suppressWarnings(
          gatherNimble(read.path = dump.path, burnin = nb, ni.block = ni,
                       base.thin = nt, max.samples.saved = max.samples.saved)
        )
        mod.check <- checkNimble(mod.out$out, Rht.required = Rht.required, neff.required = neff.required,
                                 par.ignore = par.ignore, par.dontign = par.dontign,
                                 par.fuzzy.track = par.fuzzy.track, fuzzy.threshold = fuzzy.threshold,
                                 spit.summary = TRUE, mod.nam = mod.nam)
        mod.check.result <- mod.check$result
        nblks <- mod.out$nblks
        thin.additional <- mod.out$additional.thin.rate
        mcmc.info <- c(nchains = nc, niterations = ni*nblks,
                       burnin = ifelse(nb<1, nb*ni*nblks, nb),
                       nthin = nt*thin.additional)
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
          proc$kill_tree()
          write.csv(sumTab.focal, paste0("Model_summary_PID",proc$get_pid(),".csv"))
          stop(paste0("Error: One or more parameters is not being sampled.",
                      " Check data, initial values, etc., and try again.",
                      " See 'Model_summary_PID",proc$get_pid(),
                      ".csv' for parameters missing Rhat or n.eff."))
        }
        if(length(par.fuzzy.track) > 0) {
          Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
          if(!any(names(sumTab) == "Rhat")) {
            proc$kill_tree()
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
      }
      check.log <- read.csv(check.log.file, colClasses = c("character", "numeric", "character",
                                                           "numeric", "character"))
      check.log <- check.log %>% bind_rows(
        data.frame(Model = mod.nam, Check = nchecks, Time = as.character(Sys.time()),
                   Min_waited = (nchecks * check.freq / 60), Status = status)
      )
      write.csv(check.log, check.log.file, row.names = FALSE)
      nchecks <- nchecks + 1
    }
    proc$kill_tree()
    if(!mod.check.result) {
      warn.message <- paste0("Rhat did not decrease after ", nchecks,
                            " checks. Model abandoned before reaching convergence targets.")
      mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info,
                  warning = warn.message)
      if(sav.model) R.utils::saveObject(mod, mod.nam)
      if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
    }
    unlink(dump.path, recursive = TRUE)
    gc(verbose = FALSE)
  }
