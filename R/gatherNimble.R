gatherNimble <- function(read.path, burnin, ni.block, base.thin, max.samples.saved,
                         directive.file = "") {
  # coda/mcmcOutput/stringr functions used below are resolved through this
  # package's own NAMESPACE imports, not the caller's search path, so no
  # require() is needed here (or anywhere else in this package) for a hard
  # dependency - see DESCRIPTION's Imports: and NAMESPACE's importFrom().
  cNB <- countNimbleBlocks(read.path, burnin, ni.block)
  m <- cNB$m
  if(nrow(m) == 0) {
    if(directive.file != "") writeLines("STOP", directive.file)
    stop("gatherNimble: no blocks remain after burn-in. Either burnin exceeds ",
         "everything sampled so far, or no blocks have been written yet.")
  }
  nblks <- cNB$nblks
  chns <- unique(m[,"chn"])
  burnin <- cNB$burnin
  burnin.realized <- cNB$burnin.realized
  # burnin - burnin.realized is the residual burn-in falling inside the first
  # retained block (in [0, ni.block) by construction in countNimbleBlocks()),
  # which is generally not an exact multiple of base.thin. floor() gives the
  # number of already-thinned draws whose iteration index is strictly before
  # the true burn-in cutoff (each successive draw advances the iteration index
  # by base.thin, so that many complete base.thin-steps fit inside the
  # residual); that is exactly the count to drop, with every remaining draw
  # landing at or past the cutoff.
  #
  # Flooring here, rather than leaving burnin.needed fractional, is essential:
  # row.start below is seeded from it and nrow(samp) (always a whole number,
  # checked at line ~64) is added to it once per block, so an unflored
  # fractional burnin.needed makes row.start, and therefore every element of
  # block.rows, fractional for the rest of the chain. block.rows %in% ind.sav
  # then compares fractional values against ind.sav's integers and never
  # matches, so `keep` is empty for every block, `pieces` ends up entirely
  # empty, do.call(rbind, list()) returns NULL, and as.mcmc(NULL) errors -
  # deterministically, not as a rare edge case: this fires whenever the
  # residual (in iterations) isn't an exact multiple of base.thin, which
  # includes the very first convergence check whenever burn-in hasn't been
  # fully satisfied by whole blocks yet.
  burnin.needed <- floor((burnin - burnin.realized) / base.thin)
  if(burnin.needed < 0) {
    if(directive.file != "") writeLines("STOP", directive.file)
    stop("Additional burnin needed is negative. countNimbleBlocks burned extra samples and needs to be checked.")
  }

  # Every block runs the same n.iter at the same thinning rate (runNimbleBlock()
  # is always called with n.iter = ni.block, and n.thin = base.thin is fixed at
  # the first call for the life of the chain), so every retained block has
  # exactly ni.block / base.thin rows. That means the post-burn-in chain length,
  # and therefore which rows max.samples.saved will keep, can be computed from
  # countNimbleBlocks()'s block accounting alone - before loading a single file.
  blks.per.chain <- lapply(chns, function(s) sort(unique(m[m[,"chn"] == s, "blk"])))
  n.blks.each <- vapply(blks.per.chain, length, integer(1))
  if(length(unique(n.blks.each)) > 1)
    stop("gatherNimble: chains have different numbers of surviving blocks after ",
         "countNimbleBlocks(); this should not happen and needs to be checked.")
  rows.per.block <- ni.block / base.thin
  rows.per.chain <- n.blks.each[1] * rows.per.block
  if(burnin.needed >= rows.per.chain)
    stop("gatherNimble: burnin.needed (", burnin.needed, ") meets or exceeds the ",
         rows.per.chain, " rows retained per chain. Check nb against ni and nt.")
  chain.length.now <- rows.per.chain - burnin.needed

  if(is.null(max.samples.saved)) max.samples.saved <- chain.length.now
  if(max.samples.saved < chain.length.now) {
    ind.sav <- unique(round(seq(1, chain.length.now, length.out = max.samples.saved)))
    additional.thin.rate <- chain.length.now / length(ind.sav)
  } else {
    ind.sav <- seq_len(chain.length.now)
    additional.thin.rate <- 1
  }

  # ind.sav is a set of LOCAL, 1-based post-burnin draw indices (r = 1 is the
  # first retained draw of any chain), shared identically across every chain -
  # chains are already trimmed to the same block count above, so this doesn't
  # vary by chain. Converting it to ABSOLUTE NIMBLE iteration numbers gives an
  # explicit, verifiable key for exactly which draws survived (burn-in trim and
  # any max.samples.saved-driven subsample together) - needed so gatherNimble2()
  # can derive the second monitor set's kept rows from this same selection,
  # rather than sampling independently (see gatherNimble2()'s own comments).
  #
  # Local row r's absolute iteration: burnin.realized accounts for whole
  # dropped blocks; (burnin.needed + r) is r's position counting from the start
  # of the first retained block (burnin.needed local rows before r=1 were
  # trimmed there); multiplying by base.thin converts that draw count back to
  # an iteration number, since each successive retained draw advances the
  # iteration index by base.thin. Cross-checked by execution against every
  # existing test in test-gatherNimble.R, including the fractional-burnin.needed
  # cases and a max.samples.saved cap together, using the fixture's own
  # true-iteration data column as ground truth.
  retained.iters <- burnin.realized + (burnin.needed + ind.sav) * base.thin

  # Per-row (chn, iter) key aligned with as.matrix(out)'s row order: chain-major
  # (chns in order), and within each chain, ascending iteration order - the same
  # order the block-loading loop below produces, since ind.sav is sorted
  # ascending and blocks/rows are processed in ascending order throughout.
  iter.key <- data.frame(chn = rep(chns, each = length(retained.iters)),
                         iter = rep(retained.iters, times = length(chns)))

  # Load one block at a time and immediately keep only the rows this chain
  # actually needs (the residual burn-in trim on the first retained block, and
  # the systematic subsample selected above), discarding the rest before the
  # next block is loaded. Peak memory therefore never holds more than one raw
  # block plus the shrinking accumulator, rather than every retained block for
  # every chain at once (finding F31).
  gathr <- lapply(seq_along(chns), FUN = function(ci) {
    s <- chns[ci]
    blks <- blks.per.chain[[ci]]
    row.start <- 1 - burnin.needed  # first block's local row 1 starts here
    pieces <- lapply(blks, FUN = function(b) {
      load(file = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b]))
      if(nrow(samp) != rows.per.block)
        stop("gatherNimble: block ", b, " of chain ", s, " has ", nrow(samp),
             " rows, expected ", rows.per.block, " (ni.block / base.thin). ",
             "A block dump may be from a different configuration than assumed.")
      block.rows <- row.start:(row.start + nrow(samp) - 1)
      row.start <<- row.start + nrow(samp)
      keep <- which(block.rows %in% ind.sav)
      out <- if(length(keep) == 0) NULL else samp[keep, , drop = FALSE]
      # R's generational GC does not necessarily reclaim `samp` (one full block)
      # right away, so without forcing it here garbage from earlier blocks in
      # this loop can pile up and peak memory ends up close to holding every
      # block at once anyway - exactly what this rewrite exists to avoid.
      rm(samp); gc(verbose = FALSE)
      out
    })
    pieces <- pieces[!sapply(pieces, is.null)]
    as.mcmc(do.call(rbind, pieces))
  })
  gc(verbose = FALSE)

  out <- mcmcOutput(as.mcmc.list(gathr))
  return(mget(c("out", "nblks", "additional.thin.rate", "retained.iters", "iter.key")))
}
