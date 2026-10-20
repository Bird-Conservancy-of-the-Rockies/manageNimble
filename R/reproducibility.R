# Internal helpers behind runNimble()'s `inits` forms and `seed`. Not exported.
# All of the parent-side work happens here, before any worker launches, so it
# can be unit-tested without NIMBLE or GNU parallel.

# The full generator kind (generator, normal, sample) for seeded runs. Set
# explicitly so a caller's non-default normal.kind/sample.kind can't leak into
# the streams: element 1 of .Random.seed encodes all three, and workers adopt
# whatever it encodes.
RNG_KIND <- c(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")

# Snapshot and restore the session's RNG: the three kinds RNGkind() reports,
# plus .Random.seed in the global environment (NULL if the session has never
# drawn a random number). Used with on.exit() around anything here that
# changes the RNG, so runNimble() leaves the caller's random numbers exactly as
# they were. The seed is read before RNGkind() is called, so the snapshot
# holds the seed exactly as the session had it. Restoring a NULL seed removes
# .Random.seed rather than leaving the one created in between, so an unseeded
# session stays unseeded.
saveRNGState <- function() {
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  list(kind = RNGkind(), seed = seed)
}

restoreRNGState <- function(state) {
  # A session using sample.kind = "Rounding" gets R's standard warning about
  # it on every RNGkind() call; restoring that session's own setting
  # shouldn't repeat it.
  suppressWarnings(RNGkind(state$kind[1], state$kind[2], state$kind[3]))
  if(is.null(state$seed)) {
    if(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", state$seed, envir = globalenv())
  }
}

# One L'Ecuyer-CMRG stream per chain, derived from `seed`. Streams from
# nextRNGStream() are guaranteed non-overlapping, unlike set.seed(seed + chn),
# whose seeds also collide across runs (seed 1/chain 2 == seed 2/chain 1).
# Leaves the session's RNG kind and state exactly as it found them.
deriveChainStreams <- function(seed, nc) {
  state <- saveRNGState()
  on.exit(restoreRNGState(state), add = TRUE)
  RNGkind(RNG_KIND[1], RNG_KIND[2], RNG_KIND[3])
  set.seed(seed)
  s <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  streams <- vector("list", nc)
  for(cn in seq_len(nc)) {
    s <- nextRNGStream(s)
    streams[[cn]] <- s
  }
  streams
}

# One chain's inits list: a named list (possibly empty), with no list-valued
# elements. nimbleModel() silently ignores inits names it doesn't recognise,
# list-valued ones included, so a nested list here is almost certainly a
# per-chain list passed in the wrong shape - e.g. list(chain1 = ..., chain2 =
# ...) - which would otherwise run with every inits value quietly discarded.
checkChainInits <- function(x, what) {
  if(!is.list(x))
    stop(what, " must be a named list of initial values, not ",
         class(x)[1], ".", call. = FALSE)
  if(length(x) == 0) return(invisible(NULL))
  nms <- names(x)
  if(is.null(nms) || any(is.na(nms) | nms == ""))
    stop(what, " must be a named list of initial values: every element ",
         "needs a name.", call. = FALSE)
  nested <- vapply(x, is.list, logical(1))
  if(any(nested))
    stop(what, " has list-valued element(s): ", paste(nms[nested], collapse = ", "),
         ". Initial values must be numeric vectors/arrays, and NIMBLE would ",
         "silently ignore these. For different initial values per chain, pass ",
         "an unnamed list with one named list per chain, e.g. ",
         "list(inits.chain1, inits.chain2), or a function - see ?runNimble.",
         call. = FALSE)
  invisible(NULL)
}

# How runNimble() reads `inits` - a rule on the `names` attribute, which NIMBLE
# already requires of a single inits list, not a guess from element contents:
#   function                                -> "function" (called once per chain)
#   list with no names, length > 0         -> "per.chain" (one inits list per chain)
#   list fully named, list(), or non-list  -> "single" (passed through unchanged)
#   list with some names empty             -> error
# Non-list values are passed through as "single" rather than rejected, so no
# call that NIMBLE itself accepts starts failing here.
classifyInits <- function(inits, nc) {
  if(is.function(inits)) return("function")
  if(!is.list(inits) || length(inits) == 0) return("single")
  nms <- names(inits)
  if(is.null(nms) || all(is.na(nms) | nms == "")) {
    if(length(inits) != nc)
      stop("inits is an unnamed list of length ", length(inits), ", so it is read ",
           "as one inits list per chain, but nc = ", nc, ". Supply exactly nc ",
           "per-chain lists, or name every element of a single inits list.",
           call. = FALSE)
    for(cn in seq_len(nc)) checkChainInits(inits[[cn]], paste0("inits[[", cn, "]]"))
    return("per.chain")
  }
  if(any(is.na(nms) | nms == ""))
    stop("inits has some named and some unnamed elements. Name every element ",
         "(a single inits list), or none of them (one inits list per chain).",
         call. = FALSE)
  checkChainInits(inits, "inits")
  "single"
}

# Resolves `inits` into what the workers receive:
#   inits     - saved to NimbleObjects.RData as `inits` (the original object for
#               "single"; a list of nc inits lists otherwise)
#   per.chain - saved as `inits.per.chain`; each worker then keeps only
#               inits[[chn]]
#   used      - the per-chain inits actually used, recorded as mod$inits;
#               NULL for "single"
# A function is called here, in the parent, once per chain: f(chn) if it takes
# an argument, f() otherwise. With `streams`, chain chn's call runs in
# nextRNGSubStream(streams[[chn]]), so generated inits are reproducible and
# independent of the stream the chain samples from; without, it uses the
# session RNG. Either way the session's RNG state is restored afterwards.
prepareChainInits <- function(inits, nc, streams = NULL) {
  type <- classifyInits(inits, nc)
  if(type == "single") return(list(inits = inits, per.chain = FALSE, used = NULL))
  if(type == "per.chain") return(list(inits = inits, per.chain = TRUE, used = inits))

  state <- saveRNGState()
  on.exit(restoreRNGState(state), add = TRUE)
  takes.chain <- length(formals(inits)) > 0
  out <- vector("list", nc)
  for(cn in seq_len(nc)) {
    if(!is.null(streams)) {
      RNGkind(RNG_KIND[1], RNG_KIND[2], RNG_KIND[3])
      assign(".Random.seed", nextRNGSubStream(streams[[cn]]), envir = globalenv())
    }
    val <- tryCatch(if(takes.chain) inits(cn) else inits(),
                    error = function(e)
                      stop("The inits function failed for chain ", cn, ": ",
                           conditionMessage(e), call. = FALSE))
    checkChainInits(val, paste0("The inits function's return value for chain ", cn))
    out[[cn]] <- val
  }
  list(inits = out, per.chain = TRUE, used = out)
}

# Smallest block count satisfying manual mode's stopping rule in runNimble()
# (do.gather.check), computed with the identical expression so the two can't
# disagree. Used with a seed to gather exactly this many blocks, whatever
# overshoot the 300-s polling allowed.
minGatherBlocks <- function(ni, nb, nt, max.samples.saved) {
  n <- 1
  repeat {
    nb.now <- ifelse(nb < 1, nb * ni * n, nb)
    if((ni * n - nb.now) >= max.samples.saved * nt) return(n)
    n <- n + 1
  }
}
