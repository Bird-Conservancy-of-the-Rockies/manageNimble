# manageNimble

## Package contents
```runNimble``` Primary (and only exported entry point besides the three below) function that allows parallel processing and automated convergence and sampling checks for Bayesian model fitting with Nimble. Supports an optional second monitor set (`parameters2`/`nt2`) for quantities too large to retain at the primary thinning rate - see `?runNimble`.<br>
```runNimbleBlock``` Exported. Base function called by `runNimble`'s generated worker script to initiate or continue sampling and save posterior samples as a block.<br>
```checkNimble``` Exported. Base function called by `runNimble` used to summarize posterior parameter distributions, calculate Rhat and n_effective, and assess whether convergence and MCMC sampling targets have been met.<br>
```mcmcOutputSubset``` Exported. Utility function that subsets an mcmcOutput object. Used by `checkNimble` prior to summarizing MCMC samples when not all parameters sampled need to be summarized.<br>
```countNimbleBlocks``` Internal (not exported). Used by `runNimble` to count saved NIMBLE blocks and calculate sampling progress.<br>
```gatherNimble``` Internal (not exported). Used by `runNimble` to gather and organize NIMBLE sample blocks into an mcmcOutput object.<br>
```gatherNimble2``` Internal (not exported). Companion to `gatherNimble`; gathers the second monitor set, if any, into a plain matrix with chains stacked.

## How it works

Worth reading before modifying anything - the design is not obvious from any single file.

### `ni` is a block size, not a run length

`runNimble()` does not run `ni` iterations and stop. `ni` is the size of one **block**. Chains sample a block at a time and the run continues until the convergence criteria are met or `max.tries` blocks have been sampled, so total sampling is `ni * nblks` with `nblks` unknown in advance. Under the defaults that ceiling is `ni * 10`; callers raising `max.tries` should size `ni` accordingly.

### Launching workers

`runNimble()` saves the model objects to `dump.path/NimbleObjects.RData`, generates an R script at `dump.path/ModRunScript.R` (built as a character vector and written with `writeLines`), then launches `nc` copies of it through GNU `parallel` via `processx::process$new()`. Each worker `source()`s the model file - so anything the model code defines at the top level, such as a `nimbleFunction` used inside `nimbleCode`, is available in the worker.

Each worker calls `runNimbleBlock()` to build, compile and sample one block, writes the samples to `dump.path/mod_chn<c>_<b>.RData`, then loops for further blocks. The model is built and C++-compiled **once per chain, concurrently**, inside each worker - so peak memory during startup is roughly `nc` times a single model's footprint, and each chain generates its own compiled code under `dump.path/tmp<c>`.

Blocks are written to a `.tmp` name and moved into place with `file.rename()`. Because that rename is atomic, a file visible under its final `mod_chn<c>_<b>.RData` name is always complete, which is what lets the readers below list the directory without any locking or polling.

### Coordination

The parent and its workers communicate entirely through small text files in `dump.path`:

- `runNimbleDirective.txt` holds `GO`, `PAUSE` or `STOP` and is read by every worker.
- `block<c>Status.txt` holds per-chain `GO`, `STOP` or `RESUME`.

Both sides write these atomically (temp file plus rename) and retry a malformed read, since there is no locking. The parent polls with `Sys.sleep`, and on every polling iteration also checks worker liveness and the per-chain error logs, so a dead worker surfaces immediately rather than leaving the parent waiting indefinitely.

### The convergence loop

`countNimbleBlocks()` inventories the dump files, truncates all chains to the length of the shortest, and drops blocks falling entirely within burn-in. A gap in a chain's block sequence is an error, not something to work around - it would mean pairing unrelated iteration ranges across chains.

`gatherNimble()` computes which rows survive burn-in, thinning and the `max.samples.saved` cap **analytically, from the block accounting, before opening any file**. Each block is then loaded, immediately reduced to the rows it contributes, and discarded. Peak memory therefore tracks the cap rather than the total volume sampled - which is the point of the whole block-and-dump design, and is easy to undo by accident when editing this function.

`checkNimble()` summarizes Rhat and n.eff and returns a pass/fail verdict. Parameters named in `par.ignore` are excluded from the check entirely; parameters named in `par.fuzzy.track` are excluded from the strict test and instead required to have no more than `fuzzy.threshold` of their members unconverged. Note that `par.ignore` and friends match by unanchored substring - `par.dontign` is the designed rescue for the over-matching that causes.

The parent then loops: gather, check, and if not converged write `RESUME` and let the chains add another block, up to `max.tries`. The model object is re-saved after every check, so partial output can be inspected mid-run without stopping anything.

With `check.freq = NULL` the automated checking is skipped entirely and the run instead stops once enough samples exist to fill `max.samples.saved` (which is a **per-chain** cap - the saved object holds `nc * max.samples.saved` draws).

### The second monitor set

`parameters2`/`nt2` map onto NIMBLE's `monitors2`/`thin2`, for quantities too large to retain at the primary thinning rate. They take no part in convergence checking. `gatherNimble2()` collects them once at the end, after the parent has waited for the workers to actually exit - otherwise a worker that had not yet noticed `STOP` could add a block after the primary gather, giving the two monitor sets different effective burn-in. The result is written to `<mod.nam>_monitors2.rds` as a plain `[ndraw, nparam]` matrix with chains stacked.

### Cleanup

An `on.exit()` handler sets the directive to `STOP`, kills any surviving worker and copies the logs out, on success, error or interrupt alike - all best-effort, so a failure during cleanup cannot mask the original error.

## Worker logs
`runNimble()`'s parallel workers write their console output and exit status under `dump.path/worker_logs` and `dump.path/worker_joblog.txt`, and any R-level error to `dump.path/chn<n>_error.log`. These are copied to `<mod.nam>_logs/` before `dump.path` is (optionally) deleted, whether the run succeeds or fails - see `?runNimble` for details.

See `NEWS.md` for what changed in the most recent review pass, including which changes affect results for existing callers.
