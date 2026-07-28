countNimbleBlocks <- function(read.path, burnin, ni.block) {
  ##  This could be made more flexible, but for now use a stereotyped approach:
  fl <- list.files(read.path)
  fl <- fl[which(str_detect(fl, "mod_chn"))]

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
  m <- m[which(m[,2] <= nblks),,drop = FALSE] # Make chains same length (chop off excess blocks on longer chains)

  # Drop blocks that only contain burnin
  if(burnin < 1) burnin <- ni.block * nblks * burnin
  burnin.block <- floor(burnin / ni.block)
  m <- m[which(m[,"blk"] > burnin.block),,drop = FALSE]
  burnin.realized <- burnin.block * ni.block

  return(mget(c("m", "nblks", "burnin", "burnin.realized")))
}
