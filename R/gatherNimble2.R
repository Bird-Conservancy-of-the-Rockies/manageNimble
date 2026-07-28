gatherNimble2 <- function(read.path, burnin, ni.block, save.path = NULL) {
  # Gather the second monitor set (samp2) from the block dumps written by
  # runNimbleBlock(). Companion to gatherNimble(), which gathers the primary set.
  #
  # The second set exists for quantities too large to retain at the primary
  # thinning rate, so it is handled much more simply than gatherNimble(): blocks
  # falling entirely within burn-in are dropped (via countNimbleBlocks), the
  # remainder are concatenated, and no further thinning or sample-count capping
  # is applied - the thinning was already imposed by thin2 at sampling time.
  #
  # Chains are stacked into a single [ndraw, nparam] matrix. The chain structure
  # is not preserved, because the second set is intended for post hoc
  # calculations over the posterior (posterior predictive checks, derived
  # quantities) rather than for convergence diagnostics.
  #
  # Returns NULL with a warning if no block dump contains a second set, which is
  # the expected result when runNimble() was called without parameters2.
  require(stringr)

  cNB  <- countNimbleBlocks(read.path, burnin, ni.block)
  m    <- cNB$m
  if(nrow(m) == 0) {
    warning("gatherNimble2: no blocks remain after burn-in.")
    return(NULL)
  }
  chns <- unique(m[, "chn"])

  Sys.sleep(5) # To provide time for lingering files to finish writing.

  gathr <- lapply(chns, FUN = function(s) {
    blks <- unique(m[m[, "chn"] == s, "blk"])
    lst <- lapply(blks, FUN = function(b) {
      fl <- paste0(read.path, "/", rownames(m)[m[, "chn"] == s & m[, "blk"] == b])
      x <- suppressWarnings(system2(command = "lsof", args = fl, stdout = TRUE))
      while(length(x) > 0) {
        Sys.sleep(5)
        x <- suppressWarnings(system2(command = "lsof", args = fl, stdout = TRUE))
      }
      samp2 <- NULL          # Stays NULL if the dump predates monitors2 support
      load(file = fl)        # Brings in samp and, if present, samp2
      samp2
    })
    lst <- lst[!sapply(lst, is.null)]
    if(length(lst) == 0) return(NULL)
    do.call(rbind, lst)
  })

  gathr <- gathr[!sapply(gathr, is.null)]
  if(length(gathr) == 0) {
    warning("gatherNimble2: no second monitor set found in any block dump. ",
            "Was runNimble() called with parameters2?")
    return(NULL)
  }

  out <- do.call(rbind, gathr)
  gc(verbose = FALSE)
  if(!is.null(save.path)) saveRDS(out, save.path)
  return(out)
}
