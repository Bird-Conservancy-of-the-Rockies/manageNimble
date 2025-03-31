checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
                        par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                        spit.summary = FALSE, mod.nam = "mod") {
  require(mcmcOutput)
  
  if(!is.null(par.ignore) & is.null(par.fuzzy.track)) mcmcOutput <- mcmcOutputSubset(mcmcOutput,
                                                                                     par.drop = par.ignore)
  if(!is.null(par.ignore) & !is.null(par.fuzzy.track)) mcmcOutput <- mcmcOutputSubset(mcmcOutput,
                                                                                      par.keep = c(par.dontign,
                                                                                                   par.fuzzy.track),
                                                                                      par.drop = par.ignore)
  
  s <- summary(mcmcOutput, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
  if(!any(names(s) == "Rhat")) {
    proc$kill_tree()
    stop("Rhat not calculated. Troubleshoot mcmcOuput.")
  }
  s <- s %>%
    as_tibble() %>%
    mutate(Parameter = row.names(s)) %>%
    dplyr::select(Parameter, mean:f)
  if(length(par.ignore) > 0) {
    if(length(par.dontign) == 0) {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore))
    } else {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore) &
                            !str_detect_any(s$Parameter, par.dontign))
    }
    if(length(ind.ignore) > 0) {
      s.focal <- s %>% slice(-ind.ignore)
    } else {
      s.focal <- s
    }
    if(any(is.na(s.focal$Rhat))) {
      write.csv(s.focal %>% filter(is.na(Rhat) | Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
      stop(paste0("Parameters missing Rhat. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
    }
    if(any(s.focal$Rhat %in% c(Inf, -Inf))) { # This doesn't seem to be working. Getting -Inf for Rhat with function continuing to run.
      write.csv(s.focal %>% filter(Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
      stop(paste0("Parameters with Inf or -Inf Rhats. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
    }
    result <- max(s.focal$Rhat) <= Rht.required & min(s.focal$n.eff) >= neff.required
  } else {
    result <- max(s$Rhat) <= Rht.required & min(s$n.eff) >= neff.required
  }
  if(length(par.fuzzy.track) > 0) {
    Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
    if(!any(names(s) == "Rhat")) {
      proc$kill_tree()
      write.csv(s, str_c(species, "_sum_at_fail.csv"), row.names = FALSE)
      stop("Stopped model run because Rhat not calculated.")
    }
    for(p in 1:length(par.fuzzy.track)) {
      pfuz <- par.fuzzy.track[p]
      Rht.fuzzy <- c(Rht.fuzzy,
                     s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                       pull(Rhat))
    }
    Rht.fuzzy <- Rht.fuzzy[-1]
    if(any(is.na(Rht.fuzzy))) {
      write.csv(s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                  filter(is.na(Rhat) | Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy parameters missing Rhat. Check Bad_pars_fuzzy_", mod.nam,
                  ".csv and possibly try alternative initial values or check data."))
    }
    if(any(s.focal$Rhat %in% c(Inf, -Inf))) {
      write.csv(s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                  filter(Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy parameters with Inf or -Inf Rhats. Check Bad_pars_fuzzy_", mod.nam,
                  ".csv and possibly try alternative initial values or check data."))
    }
    if(((sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = TRUE) + sum(is.na(Rht.fuzzy))) / length(Rht.fuzzy)) >
        (length(Rht.fuzzy) * fuzzy.threshold)) result <- FALSE
  }
  if(spit.summary) {
    return(mget(c("result", "s")))
  } else {
    return(mget(c("result")))
  }
}
