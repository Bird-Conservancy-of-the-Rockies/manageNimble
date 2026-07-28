mcmcOutputSubset <- function(mcmcOutput, par.keep = c(), par.drop = c()) {
  ## here, "mcmcOutput" is the incoming mcmcOutput object
  ## "mcmcOutput2" is the outgoing (subsetted) mcmcOutput object

  # mcmcOutput/coda functions used below are resolved through this package's
  # own NAMESPACE imports; no library()/require() is needed here.
  nc <- attr(mcmcOutput, "nChains")
  chn.iter <- nrow(mcmcOutput)/nc

  if(nc == 1) {
    # 1:(nc-1) is 1:0 = c(1,0) for nc == 1, which builds a nonsensical
    # descending-range index matrix. Handle the single-chain case directly.
    m <- matrix(c(1, chn.iter), nrow = 2)
  } else {
    m <- chn.iter * matrix(c(0, rep(1:(nc-1), each = 2), nc), nrow = 2) + c(1,0)
  }

  ## need to replace actual periods with escaped periods:
  dot.escape <- function (string) {
    gsub(x = string, pattern = "\\.", replacement = "\\\\.") 
  }
  
  if(is.null(par.keep) & is.null(par.drop))
    stop("One or both of 'par.keep' and 'par.drop' need to be specified.")
  
  if(any(par.keep %in% par.drop))
    warning("One or more element of 'par.keep' is in 'par.drop'. All parameters in 'par.keep' will be summarized even if they appear in 'par.drop'.")
  
  if(!is.null(par.keep)) {
    if(length(par.keep) > 0) {
      want.cols <- sapply(par.keep, FUN = function (s) {
        base <- dot.escape(s)
        mult <- paste0(base, "\\[")  ## maybe a vector or array
        singl <- paste0(base, "$")   ## maybe a scalar
        c(grep(x = colnames(mcmcOutput), pattern = mult, value = TRUE),
          grep(x = colnames(mcmcOutput), pattern = singl, value = TRUE))
      })
      want.cols <- unlist(want.cols)
      if(length(want.cols) == 0)
        stop("mcmcOutputSubset: par.keep matched no columns. par.keep = ",
             paste(par.keep, collapse = ", "), "; available columns: ",
             paste(colnames(mcmcOutput), collapse = ", "))
    } else {
      want.cols <- character(0)
    }
  }

  if(!is.null(par.drop)) {
    noWant <- sapply(par.drop, FUN = function (s) {
      base <- dot.escape(s)
      mult <- paste0(base, "\\[")  ## maybe a vector or array
      singl <- paste0(base, "$")   ## maybe a scalar
      c(grep(x = colnames(mcmcOutput), pattern = mult, value = TRUE),
        grep(x = colnames(mcmcOutput), pattern = singl, value = TRUE))
    })
    noWant <- unlist(noWant)
    
    if(!is.null(par.keep)) {
      want.cols <- c(want.cols, colnames(mcmcOutput)[!colnames(mcmcOutput) %in% noWant])
    } else {
      want.cols <- colnames(mcmcOutput)[!colnames(mcmcOutput) %in% noWant]
    }
  }
  
  want.cols <- unique(want.cols) # Shouldn't be an issue, but just in case
  if(length(want.cols) == 0)
    stop("mcmcOutputSubset: no columns remain after applying par.keep/par.drop.")

  full <- mcmcOutput[,]  # materialize once; every chain slice below reads this
  mcmcOutput.lst <- lapply(as.data.frame(m), FUN = function(v) {
    as.mcmc(full[v[1]:v[2], want.cols, drop = FALSE])
  })

  names(mcmcOutput.lst) <- NULL
  
  mcmcOutput2 <- mcmcOutput(as.mcmc.list(mcmcOutput.lst))
  stopifnot(all(dimnames(mcmcOutput2)[[2]] == want.cols))
  
  return(mcmcOutput2)
}
