checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
                        par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                        spit.summary = FALSE, mod.nam = "mod", directive.file = "") {
  # par.dontign is a rescue list for parameters that par.ignore's substring
  # matching would otherwise sweep up unintentionally (see man/checkNimble.Rd).
  # It must be applied whenever par.ignore is, not only when par.fuzzy.track is
  # also supplied - previously it was silently a no-op in the (common)
  # par.ignore-only case.
  if(length(par.ignore) > 0) {
    mcmcOutput <- mcmcOutputSubset(mcmcOutput,
                                   par.keep = c(par.dontign, par.fuzzy.track),
                                   par.drop = par.ignore)
  }

  s <- summary(mcmcOutput, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
  if(!any(names(s) == "Rhat")) {
    if(directive.file != "") writeLines("STOP", directive.file)
    stop("Rhat not calculated. Troubleshoot mcmcOuput.")
  }
  s <- s %>%
    as_tibble() %>%
    mutate(Parameter = row.names(s)) %>%
    dplyr::select(Parameter, mean:f)

  # Parameters governed by the strict Rhat / n.eff criterion. par.ignore is
  # excluded because it should not be checked at all; par.fuzzy.track is
  # excluded because those parameters are governed instead by fuzzy.threshold
  # below.
  par.base    <- str_split(s$Parameter, "\\[", simplify = TRUE)[,1]
  ind.exclude <- which(par.base %in% c(par.ignore, par.fuzzy.track))
  s.focal     <- if(length(ind.exclude) > 0) s %>% slice(-ind.exclude) else s

  if(nrow(s.focal) == 0)
    stop("No parameters remain for the strict convergence check after removing ",
         "par.ignore and par.fuzzy.track. At least one parameter must be checked ",
         "strictly, or the run can never be declared converged.")

  # NA in either Rhat or n.eff is fatal, for both focal and fuzzy-tracked
  # parameters. Without this, max()/min() silently propagate the NA into
  # `result`, which crashes downstream `if(!mod.check.result)` checks in
  # runNimble() ("missing value where TRUE/FALSE needed") - and before that
  # crash, a parameter that was never actually sampled would not have been
  # caught at all.
  if(any(is.na(s.focal$Rhat)) | any(is.na(s.focal$n.eff))) {
    write.csv(s.focal %>% filter(is.na(Rhat) | is.na(n.eff) | Rhat %in% c(Inf, -Inf)),
              paste0("Bad_pars_", mod.nam, ".csv"))
    stop(paste0("Parameters missing Rhat or n.eff. Check Bad_pars_", mod.nam,
                ".csv and possibly try alternative initial values or check data."))
  }
  if(any(s.focal$Rhat %in% c(Inf, -Inf))) {
    write.csv(s.focal %>% filter(Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
    stop(paste0("Parameters with Inf or -Inf Rhats. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
  }
  result <- max(s.focal$Rhat) <= Rht.required & min(s.focal$n.eff) >= neff.required

  prp.fuzzy.not.converged <- NA_real_
  if(length(par.fuzzy.track) > 0) {
    # Rows of s belonging to any tracked family, matched on the base name (same
    # convention as par.ignore above) so every family is covered, not just one.
    ind.fuzzy <- which(str_split(s$Parameter, "\\[", simplify = TRUE)[,1] %in%
                         par.fuzzy.track)
    s.fuzzy   <- s %>% slice(ind.fuzzy)

    if(nrow(s.fuzzy) == 0)
      stop("par.fuzzy.track matched no parameters (", paste(par.fuzzy.track, collapse = ", "),
           "). Check for a typo, or that this model actually monitors these parameters.")

    Rht.fuzzy  <- s.fuzzy$Rhat
    neff.fuzzy <- s.fuzzy$n.eff

    if(any(is.na(Rht.fuzzy)) | any(is.na(neff.fuzzy))) {
      write.csv(s.fuzzy %>% filter(is.na(Rhat) | is.na(n.eff) | Rhat %in% c(Inf, -Inf)),
                paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy-tracked parameters missing Rhat or n.eff. Check Bad_pars_fuzzy_",
                  mod.nam, ".csv and possibly try alternative initial values or check data."))
    }
    if(any(Rht.fuzzy %in% c(Inf, -Inf))) {
      write.csv(s.fuzzy %>% filter(Rhat %in% c(Inf, -Inf)),
                paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy parameters with Inf or -Inf Rhats. Check Bad_pars_fuzzy_", mod.nam,
                  ".csv and possibly try alternative initial values or check data."))
    }

    # Proportion of fuzzy-tracked parameters failing EITHER criterion: a high
    # Rhat or a below-target n.eff. fuzzy.threshold governs both, with the same
    # tolerance. NA is excluded from this calculation because it is fatal above.
    prp.fuzzy.not.converged <-
      sum(Rht.fuzzy > Rht.required | neff.fuzzy < neff.required) / length(Rht.fuzzy)
    if(prp.fuzzy.not.converged > fuzzy.threshold) result <- FALSE
  }

  out <- list(result = result, prp.fuzzy.not.converged = prp.fuzzy.not.converged)
  if(spit.summary) {
    # mcmcOutputSubset()'s par.keep/par.drop matching is intentionally
    # unanchored (see its own tests/docs) - so a par.ignore'd parameter can be
    # rescued back into `s` as a side effect of an unrelated par.keep entry
    # (e.g. "sd.pnt_occupancy" rescued by "pnt_occupancy" appearing in
    # par.fuzzy.track). Re-apply the same exact base-name drop already used for
    # s.focal above, so par.ignore actually means "not shown", while
    # par.dontign still rescues anything a caller explicitly wants kept.
    drop.out <- which(par.base %in% par.ignore & !(par.base %in% par.dontign))
    out$s <- if(length(drop.out) > 0) s %>% slice(-drop.out) else s
  }
  return(out)
}
