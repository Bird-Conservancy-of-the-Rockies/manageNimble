checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
                        par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                        spit.summary = FALSE, mod.nam = "mod", directive.file = "") {
  require(mcmcOutput)
  
  if(!is.null(par.ignore) & is.null(par.fuzzy.track)) mcmcOutput <- mcmcOutputSubset(mcmcOutput,
                                                                                     par.drop = par.ignore)
  if(!is.null(par.ignore) & !is.null(par.fuzzy.track)) mcmcOutput <- mcmcOutputSubset(mcmcOutput,
                                                                                      par.keep = c(par.dontign,
                                                                                                   par.fuzzy.track),
                                                                                      par.drop = par.ignore)
  
  s <- summary(mcmcOutput, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
  if(!any(names(s) == "Rhat")) {
    # proc$kill_tree()
    if(directive.file != "") writeLines("STOP", directive.file)
    stop("Rhat not calculated. Troubleshoot mcmcOuput.")
  }
  s <- s %>%
    as_tibble() %>%
    mutate(Parameter = row.names(s)) %>%
    dplyr::select(Parameter, mean:f)
  # Parameters governed by the strict Rhat / n.eff criterion. Both par.ignore and
  # par.fuzzy.track are excluded here: par.ignore because it should not be
  # checked at all, par.fuzzy.track because those parameters are governed instead
  # by the fuzzy.threshold tolerance below.
  #
  # Previously the fuzzy-tracked parameters stayed in the strict test, which made
  # the fuzzy block unreachable: if the strict test passed then no fuzzy
  # parameter could be unconverged, and if it failed then result was already
  # FALSE. This branch also used to define s.focal only when par.ignore was
  # non-empty, so the fuzzy block below failed with "object 's.focal' not found"
  # whenever par.fuzzy.track was used without par.ignore.
  par.base    <- str_split(s$Parameter, "\\[", simplify = TRUE)[,1]
  ind.exclude <- which(par.base %in% c(par.ignore, par.fuzzy.track))
  s.focal     <- if(length(ind.exclude) > 0) s %>% slice(-ind.exclude) else s

  if(nrow(s.focal) == 0)
    stop("No parameters remain for the strict convergence check after removing ",
         "par.ignore and par.fuzzy.track. At least one parameter must be checked ",
         "strictly, or the run can never be declared converged.")

  if(any(is.na(s.focal$Rhat))) {
    write.csv(s.focal %>% filter(is.na(Rhat) | Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
    stop(paste0("Parameters missing Rhat. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
  }
  if(any(s.focal$Rhat %in% c(Inf, -Inf))) { # This doesn't seem to be working. Getting -Inf for Rhat with function continuing to run.
    write.csv(s.focal %>% filter(Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
    stop(paste0("Parameters with Inf or -Inf Rhats. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
  }
  result <- max(s.focal$Rhat) <= Rht.required & min(s.focal$n.eff) >= neff.required
  if(length(par.fuzzy.track) > 0) {
    # Rows of s belonging to any tracked family. Built once, before the checks,
    # so that the diagnostic CSVs below cover every tracked family. Previously
    # these filtered on `pfuz`, which after the loop held only the LAST element
    # of par.fuzzy.track, so diagnostics silently omitted the other families.
    # Matching on the base name (split at "[") is the same convention used for
    # par.ignore above.
    ind.fuzzy <- which(str_split(s$Parameter, "\\[", simplify = TRUE)[,1] %in%
                         par.fuzzy.track)
    s.fuzzy   <- s %>% slice(ind.fuzzy)
    Rht.fuzzy <- s.fuzzy$Rhat

    # NAs are counted toward the unconverged proportion rather than being fatal.
    # The proportion below explicitly adds sum(is.na(Rht.fuzzy)) to its
    # numerator, which was unreachable while this was a stop(). Counting them is
    # also the point of fuzzy tracking: a data-poor species can legitimately
    # yield an uncomputable Rhat, and that should not halt a run that is
    # otherwise converged. The diagnostic file is still written so the affected
    # parameters remain visible. Change warning() back to stop() to restore the
    # previous fatal behaviour.
    if(any(is.na(Rht.fuzzy))) {
      write.csv(s.fuzzy %>% filter(is.na(Rhat) | Rhat %in% c(Inf, -Inf)),
                paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      warning(paste0(sum(is.na(Rht.fuzzy)), " of ", length(Rht.fuzzy),
                     " fuzzy-tracked parameters are missing Rhat and are being ",
                     "counted as unconverged. See Bad_pars_fuzzy_", mod.nam, ".csv."))
    }
    # Was `s.focal$Rhat`, which is a hard error when par.ignore is empty: s.focal
    # is only created inside the `if(length(par.ignore) > 0)` branch above, so
    # using par.fuzzy.track without par.ignore failed with "object 's.focal' not
    # found". It also tested the wrong parameters when it did run.
    if(any(Rht.fuzzy %in% c(Inf, -Inf))) {
      write.csv(s.fuzzy %>% filter(Rhat %in% c(Inf, -Inf)),
                paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy parameters with Inf or -Inf Rhats. Check Bad_pars_fuzzy_", mod.nam,
                  ".csv and possibly try alternative initial values or check data."))
    }

    # Proportion of tracked parameters that have not converged.
    # NOTE: the comparison below was previously against
    # `length(Rht.fuzzy) * fuzzy.threshold`, a count, while the left side is a
    # proportion. For any family with more than 1 / fuzzy.threshold members the
    # right side exceeded 1 and the condition could never be TRUE, so the fuzzy
    # check never fired. See also the note on par.fuzzy.track and s.focal below.
    prp.not.converged <-
      (sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = TRUE) +
         sum(is.na(Rht.fuzzy))) / length(Rht.fuzzy)
    if(prp.not.converged > fuzzy.threshold) result <- FALSE
  }
  if(spit.summary) {
    return(mget(c("result", "s")))
  } else {
    return(mget(c("result")))
  }
}
