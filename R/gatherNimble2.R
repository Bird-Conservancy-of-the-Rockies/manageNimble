gatherNimble2 <- function(read.path, burnin, ni.block, nt2 = 1, base.thin,
                          retained.iters, save.path = NULL) {
  # Gather the second monitor set (samp2) from the block dumps written by
  # runNimbleBlock(). Companion to gatherNimble(), which gathers the primary set.
  #
  # The kept rows are derived directly from gatherNimble()'s own selection
  # (retained.iters - the exact absolute NIMBLE iteration numbers it retained
  # for the primary set, see gatherNimble()'s own comments) rather than
  # sampled independently: a row is kept here iff its absolute iteration is in
  # retained.iters. Since nt2 must be a positive integer multiple of base.thin
  # (validated below), every nt2-spaced iteration is guaranteed to also be an
  # base.thin-spaced iteration - i.e. a candidate member of the primary set's
  # full (pre-cap) retained-iteration sequence - so this intersection is exact,
  # not a resampling. This guarantees two properties together: the ratio of
  # kept second-set rows to kept primary-set rows is always nt2/base.thin
  # regardless of how many retry blocks a run needed (retained.iters already
  # reflects whatever burn-in trim and max.samples.saved-driven subsample
  # gatherNimble() actually applied), and every kept row's (chn, iter) is an
  # explicit, verifiable join key back to the corresponding primary-set row -
  # see the returned iter.key.
  #
  # Chains are stacked into a single [ndraw, nparam] matrix. The chain structure
  # is not preserved, because the second set is intended for post hoc
  # calculations over the posterior (posterior predictive checks, derived
  # quantities) rather than for convergence diagnostics; the returned iter.key
  # carries chain identity back for any caller that needs it.
  #
  # Returns NULL with a warning if no block dump contains a second set, which is
  # the expected result when runNimble() was called without parameters2.
  if(nt2 <= 0 || nt2 %% 1 != 0 || nt2 %% base.thin != 0)
    stop("gatherNimble2: nt2 (", nt2, ") must be a positive integer multiple of ",
         "the primary thinning rate base.thin (", base.thin, "). Every second-set ",
         "draw's absolute iteration must be guaranteed to also be a candidate ",
         "member of the primary set's retained-iteration sequence, or the two ",
         "monitor sets cannot be kept in verifiable correspondence.")

  cNB  <- countNimbleBlocks(read.path, burnin, ni.block)
  m    <- cNB$m
  if(nrow(m) == 0) {
    warning("gatherNimble2: no blocks remain after burn-in.")
    return(NULL)
  }
  chns <- unique(m[, "chn"])

  # runNimbleBlock() writes each block via a temp file + atomic rename, so any
  # file visible under its final mod_chn<c>_<b>.RData name is guaranteed
  # complete - no polling for a still-open file descriptor is needed here.
  gathr <- lapply(chns, FUN = function(s) {
    blks <- unique(m[m[, "chn"] == s, "blk"])
    lst <- lapply(blks, FUN = function(b) {
      fl <- paste0(read.path, "/", rownames(m)[m[, "chn"] == s & m[, "blk"] == b])
      samp2 <- NULL          # Stays NULL if the dump predates monitors2 support
      load(file = fl)        # Brings in samp and, if present, samp2
      if(is.null(samp2)) {
        rm(samp); gc(verbose = FALSE)
        return(NULL)
      }
      # Row k (1-based) of this block's samp2 was recorded at absolute NIMBLE
      # iteration (b-1)*ni.block + k*nt2 - thinning runs continuously across
      # block boundaries (runNimbleBlock() reuses the compiled sampler with
      # reset = FALSE on every continuation block), so this is exact, the same
      # formula pattern gatherNimble() uses for the primary set at base.thin.
      abs.iter <- (b - 1) * ni.block + seq_len(nrow(samp2)) * nt2
      keep <- which(abs.iter %in% retained.iters)
      out <- if(length(keep) == 0) NULL else
        list(mat = samp2[keep, , drop = FALSE], iter = abs.iter[keep])
      # load() also brings in `samp` (the primary set, potentially large), which
      # this function never uses. Drop it and force reclaim rather than letting
      # it accumulate across blocks, same rationale as gatherNimble()'s own loop.
      rm(samp, samp2); gc(verbose = FALSE)
      out
    })
    lst <- lst[!sapply(lst, is.null)]
    if(length(lst) == 0) return(NULL)
    # chn carried explicitly (not inferred positionally) so a chain dropped
    # below for contributing nothing can't silently shift another chain's id.
    list(chn = s, mat = do.call(rbind, lapply(lst, `[[`, "mat")),
         iter = unlist(lapply(lst, `[[`, "iter")))
  })

  gathr <- gathr[!sapply(gathr, is.null)]
  if(length(gathr) == 0) {
    warning("gatherNimble2: no second monitor set found in any block dump. ",
            "Was runNimble() called with parameters2?")
    return(NULL)
  }

  # A chain contributing fewer draws than the others usually means one of its
  # block dumps predates parameters2, or lacks a second set for some other
  # reason - every chain draws from the identical retained.iters target set, so
  # under normal operation (every block's samp2 present) all chains contribute
  # exactly the same count. Since chain identity isn't kept in `out`, this would
  # otherwise silently weight chains unequally in any downstream calculation.
  row.counts <- vapply(gathr, function(x) nrow(x$mat), integer(1))
  if(length(unique(row.counts)) > 1)
    warning("gatherNimble2: chains contributed unequal numbers of draws (",
            paste(row.counts, collapse = ", "), "). A block dump is likely ",
            "missing its second monitor set.")

  out <- do.call(rbind, lapply(gathr, `[[`, "mat"))
  iter.key <- data.frame(chn  = rep(vapply(gathr, `[[`, numeric(1), "chn"), row.counts),
                         iter = unlist(lapply(gathr, `[[`, "iter")))
  gc(verbose = FALSE)
  result <- list(out = out, iter.key = iter.key)
  if(!is.null(save.path)) saveRDS(result, save.path)
  return(result)
}
