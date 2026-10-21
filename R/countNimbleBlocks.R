countNimbleBlocks <- function(read.path, burnin, ni.block, max.nblks = NULL) {
  ##  This could be made more flexible, but for now use a stereotyped approach:
  fl <- list.files(read.path)
  # Anchored to end in ".RData" so a still-being-written mod_chn<c>_<b>.RData.tmp
  # (runNimbleBlock() writes each block via temp file + atomic rename) is never
  # counted as a complete dump - a readdir caught mid-rename can list both names
  # at once, which would otherwise look like the same block number twice. A
  # malformed ".RData" name (e.g. a truncated write) still matches, reaching the
  # NA-parsing check below rather than being silently dropped.
  fl <- fl[which(str_detect(fl, "^mod_chn.*\\.RData$"))]

  m <- cbind(chn = str_split(fl, "_", simplify = TRUE)[,2] %>%
               str_sub(4, -1) %>% as.integer,
             blk = str_split(str_split(fl, "_", simplify = TRUE)[,3],
                             "\\.", simplify = TRUE)[,1] %>% as.integer)
  rownames(m) <- fl
  if(any(is.na(m)))
    stop("countNimbleBlocks: could not parse a chain or block number from ",
         "file name(s): ", paste(fl[apply(is.na(m), 1, any)], collapse = ", "),
         ". Expected the pattern mod_chn<chain>_<block>.RData.")
  m <- m[order(m[,"chn"],m[,"blk"]),,drop = FALSE]
  for(chn in unique(m[,"chn"])) { # Block numbers are sometimes gapped (a worker died mid-run,
                                  # or a dump was lost). Renumbering here would silently pair
                                  # unrelated iteration ranges across chains, so refuse instead.
    ind.blks <- unname(m[which(m[,"chn"] == chn), "blk"])
    if(!identical(ind.blks, seq_len(length(ind.blks))))
      stop("countNimbleBlocks: chain ", chn, " is missing block dump(s) ",
           paste(setdiff(seq_len(max(ind.blks)), ind.blks), collapse = ", "),
           ". Refusing to renumber - block number encodes iteration range, and ",
           "renumbering would silently align mismatched iterations across chains. ",
           "Investigate why the block(s) are missing (see worker logs) before resuming.")
  }
  nblks <- min(tapply(m[,2], m[,1], max))
  # max.nblks only ever lowers the count, to ignore blocks beyond a fixed
  # target - runNimble()'s seeded manual mode uses it so that polling overshoot
  # can't change what gets gathered.
  if(!is.null(max.nblks)) nblks <- min(nblks, as.integer(max.nblks))
  m <- m[which(m[,2] <= nblks),,drop = FALSE] # Make chains same length (chop off excess blocks on longer chains)

  # Drop blocks that only contain burnin
  if(burnin < 1) burnin <- ni.block * nblks * burnin
  burnin.block <- floor(burnin / ni.block)
  m <- m[which(m[,"blk"] > burnin.block),,drop = FALSE]
  burnin.realized <- burnin.block * ni.block

  return(mget(c("m", "nblks", "burnin", "burnin.realized")))
}
