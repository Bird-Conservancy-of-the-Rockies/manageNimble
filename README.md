# manageNimble

## Package contents
```runNimble``` Primary function that allows parallel processing and automated convergences and sampling checks for Bayesian model fitting with Nimble<br>
```runNimbleBlock``` Base function called by `runNimble` to initiate or continue sampling and save posterior samples as a block<br>
```countNimbleBlocks``` Base function called by `runNimble` used to count saved NIMBLE blocks and calculate sampling progress<br>
```gatherNimble``` Base function called by `runNimble` used to gather and organize NIMBLE sample blocks into an mcmcOutput object<br>
```checkNimble``` Base function called by `runNimble` used to summarize posterior parameter distributions, calculate Rhat and n_effective, and assess whether convergence and MCMC sampling targets have been met<br>
```mcmcOutputSubset``` Base function that subsets mcmcOutput. Used by `checkNimble` prior to summarizing MCMC samples stored as mcmcOutput objects when not all parameters sampled need to be summarized.