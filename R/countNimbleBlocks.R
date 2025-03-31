countNimbleBlocks <- function(read.path, burnin, ni.block) {
  require(stringr)
  
  ##  This could be made more flexible, but for now use a stereotyped approach:
  fl <- list.files(read.path)
  fl <- fl[which(str_detect(fl, "mod_chn"))]
  
  m <- cbind(chn = str_split(fl, "_", simplify = TRUE)[,2] %>%
               str_sub(4, -1) %>% as.integer,
             blk = str_split(str_split(fl, "_", simplify = TRUE)[,3],
                             "\\.", simplify = TRUE)[,1] %>% as.integer)
  rownames(m) <- fl
  m <- m[order(m[,"chn"],m[,"blk"]),,drop = FALSE]
  for(chn in 1:max(m[,"chn"])) { # It seems that in some cases block numbers are getting skipped,
                                  # so adding this to make sure blocks are numbered sequentially.
    ind.blks <- m[which(m[,"chn"] == chn), "blk"]
    if(any(ind.blks != 1:length(ind.blks)))
      m[which(m[,"chn"] == chn), "blk"] <- 1:length(ind.blks)
  }
  nblks <- min(tapply(m[,2], m[,1], max))
  m <- m[which(m[,2] <= nblks),] # Make chains same length (chop off excess blocks on longer chains)
  
  # Drop blocks that only contain burnin
  if(burnin < 1) burnin <- ni.block * nblks * burnin
  burnin.block <- floor(burnin / ni.block)
  m <- m[which(m[,"blk"] > burnin.block),,drop = FALSE]
  burnin.realized <- burnin.block * ni.block
  
  return(mget(c("m", "nblks", "burnin", "burnin.realized")))
}
