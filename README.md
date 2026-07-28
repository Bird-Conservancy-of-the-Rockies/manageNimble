# manageNimble

## Package contents
```runNimble``` Primary (and only exported entry point besides the three below) function that allows parallel processing and automated convergence and sampling checks for Bayesian model fitting with Nimble. Supports an optional second monitor set (`parameters2`/`nt2`) for quantities too large to retain at the primary thinning rate - see `?runNimble`.<br>
```runNimbleBlock``` Exported. Base function called by `runNimble`'s generated worker script to initiate or continue sampling and save posterior samples as a block.<br>
```checkNimble``` Exported. Base function called by `runNimble` used to summarize posterior parameter distributions, calculate Rhat and n_effective, and assess whether convergence and MCMC sampling targets have been met.<br>
```mcmcOutputSubset``` Exported. Utility function that subsets an mcmcOutput object. Used by `checkNimble` prior to summarizing MCMC samples when not all parameters sampled need to be summarized.<br>
```countNimbleBlocks``` Internal (not exported). Used by `runNimble` to count saved NIMBLE blocks and calculate sampling progress.<br>
```gatherNimble``` Internal (not exported). Used by `runNimble` to gather and organize NIMBLE sample blocks into an mcmcOutput object.<br>
```gatherNimble2``` Internal (not exported). Companion to `gatherNimble`; gathers the second monitor set, if any, into a plain matrix with chains stacked.

## Worker logs
`runNimble()`'s parallel workers write their console output and exit status under `dump.path/worker_logs` and `dump.path/worker_joblog.txt`, and any R-level error to `dump.path/chn<n>_error.log`. These are copied to `<mod.nam>_logs/` before `dump.path` is (optionally) deleted, whether the run succeeds or fails - see `?runNimble` for details.

See `NEWS.md` for what changed in the most recent review pass, including which changes affect results for existing callers.
