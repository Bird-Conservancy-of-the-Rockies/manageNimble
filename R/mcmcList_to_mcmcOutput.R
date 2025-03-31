mcmcList_to_mcmcOutput <- function(mcmcList) {
  require(coda)
  require(mcmcOutput)
  require(dplyr)
  
  mcmcList <- coda::as.mcmc.list(lapply(mcmcList, coda::as.mcmc))
  mcmcList <- mcmcOutput(mcmcList)
  sumTab <- summary(mcmcList, MCEpc = FALSE, Rhat = TRUE, n.eff = TRUE,
                    f = FALSE, overlap0 = FALSE, verbose = FALSE)
  sumTab <- sumTab %>%
    as_tibble() %>%
    mutate(Parameter = row.names(sumTab)) %>%
    dplyr::select(Parameter, mean:n.eff)
  mcmcList <- list(mcmcOutput = mcmcList, summary = sumTab)
  return(mcmcList)
}
