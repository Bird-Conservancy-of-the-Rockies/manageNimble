gatherNimble2 <- function(read.path, burnin, ni.block, nt2 = 1, save.path = NULL,
                          max.rows = NULL) {
  # Gather the second monitor set (samp2) from the block dumps written by
  # runNimbleBlock(). Companion to gatherNimble(), which gathers the primary set.
  #
  # The second set exists for quantities too large to retain at the primary
  # thinning rate, so it is handled more simply than gatherNimble(): blocks
  # falling entirely within burn-in are dropped (via countNimbleBlocks), any
  # residual burn-in within the first retained block is trimmed per chain (at
  # nt2's granularity), and max.rows (if supplied) systematically subsamples
  # the concatenated result - no thinning is otherwise applied, since thin2
  # already set the rate at sampling time.
  #
  # Chains are stacked into a single [ndraw, nparam] matrix. The chain structure
  # is not preserved, because the second set is intended for post hoc
  # calculations over the posterior (posterior predictive checks, derived
  # quantities) rather than for convergence diagnostics.
  #
  # Returns NULL with a warning if no block dump contains a second set, which is
  # the expected result when runNimble() was called without parameters2.
  cNB  <- countNimbleBlocks(read.path, burnin, ni.block)
  m    <- cNB$m
  if(nrow(m) == 0) {
    warning("gatherNimble2: no blocks remain after burn-in.")
    return(NULL)
  }
  chns <- unique(m[, "chn"])

  burnin.needed <- (cNB$burnin - cNB$burnin.realized) / nt2
  if(burnin.needed < 0)
    stop("gatherNimble2: additional burnin needed is negative. countNimbleBlocks ",
         "burned extra samples and needs to be checked.")

  # runNimbleBlock() writes each block via a temp file + atomic rename, so any
  # file visible under its final mod_chn<c>_<b>.RData name is guaranteed
  # complete - no polling for a still-open file descriptor is needed here.
  gathr <- lapply(chns, FUN = function(s) {
    blks <- unique(m[m[, "chn"] == s, "blk"])
    lst <- lapply(blks, FUN = function(b) {
      fl <- paste0(read.path, "/", rownames(m)[m[, "chn"] == s & m[, "blk"] == b])
      samp2 <- NULL          # Stays NULL if the dump predates monitors2 support
      load(file = fl)        # Brings in samp and, if present, samp2
      out <- samp2
      # load() also brings in `samp` (the primary set, potentially large), which
      # this function never uses. Drop it and force reclaim rather than letting
      # it accumulate across blocks, same rationale as gatherNimble()'s own loop.
      rm(samp); gc(verbose = FALSE)
      out
    })
    lst <- lst[!sapply(lst, is.null)]
    if(length(lst) == 0) return(NULL)
    mat <- do.call(rbind, lst)
    if(burnin.needed > 0) {
      if(burnin.needed >= nrow(mat))
        stop("gatherNimble2: burnin.needed (", burnin.needed, ") meets or exceeds ",
             "the ", nrow(mat), " draws available for chain ", s,
             ". Check nt2 against burnin and ni.block.")
      mat <- mat[-seq_len(burnin.needed), , drop = FALSE]
    }
    mat
  })

  gathr <- gathr[!sapply(gathr, is.null)]
  if(length(gathr) == 0) {
    warning("gatherNimble2: no second monitor set found in any block dump. ",
            "Was runNimble() called with parameters2?")
    return(NULL)
  }

  # A chain contributing fewer draws than the others usually means one of its
  # block dumps predates parameters2, or lacks a second set for some other
  # reason. Since chain identity isn't kept in the output, this would otherwise
  # silently weight chains unequally in any downstream calculation.
  row.counts <- sapply(gathr, nrow)
  if(length(unique(row.counts)) > 1)
    warning("gatherNimble2: chains contributed unequal numbers of draws (",
            paste(row.counts, collapse = ", "), "). A block dump is likely ",
            "missing its second monitor set.")

  out <- do.call(rbind, gathr)
  if(!is.null(max.rows) && max.rows < nrow(out)) {
    ind.sav <- unique(round(seq(1, nrow(out), length.out = max.rows)))
    out <- out[ind.sav, , drop = FALSE]
  }
  gc(verbose = FALSE)
  if(!is.null(save.path)) saveRDS(out, save.path)
  return(out)
}
