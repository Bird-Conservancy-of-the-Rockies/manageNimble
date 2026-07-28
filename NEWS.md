# manageNimble 0.0

This release implements every finding from a full code review of `R/` (34
findings, commit `68b7708` and earlier). Full technical detail — how each was
found, whether by execution, static analysis, or reading — was captured during
the review and is summarized below by finding ID. Findings marked
**CHANGES RESULTS** affect what a caller using the relevant argument gets back;
everything else is a crash/reliability/packaging fix with no effect on any run
that was working correctly before.

Backward compatibility: unless marked **CHANGES RESULTS**, no existing call to
`runNimble()`, `checkNimble()`, `mcmcOutputSubset()`, or `runNimbleBlock()`
changes its output. `parameters2`/`nt2` remain fully optional and inert at their
defaults.

## Silent wrongness (highest priority — could make a converged verdict wrong)

- **F1 — CHANGES RESULTS.** `SamplerSourcePath` was passed as an element of
  `list(...)`, not as an argument to `runNimbleBlock()`, so it was silently
  ignored: neither `nimbleHMC` nor the custom sampler script it named was ever
  loaded, and every such run fitted with NIMBLE's default samplers regardless
  of what was requested. Fixed by passing it as a real argument. **Any
  analysis that supplied `SamplerSourcePath` will now actually use the
  requested samplers** — results from before and after this fix are not
  comparable.

- **F2 — CHANGES RESULTS.** `par.dontign` (the documented rescue for
  `par.ignore`'s substring over-matching — see `par.ignore = "p"` also
  matching `beta.p[1]`/`mu.p`) was only applied when `par.fuzzy.track` was
  *also* supplied. In the far more common `par.ignore`-only case it was
  silently a no-op. Fixed: it now always applies whenever `par.ignore` is
  supplied. **Any caller using `par.dontign` without `par.fuzzy.track`
  will now have more parameters checked than before** — convergence becomes
  harder, correctly, for those calls.

- **F3 — CHANGES RESULTS.** `countNimbleBlocks()` renumbered gapped block
  sequences (e.g. chain 1 missing block 4) instead of reporting the gap,
  silently pairing unrelated iteration ranges across chains and computing
  Rhat over data that didn't cover the same span of the run. It now **stops**
  with the chain and missing block number instead. A run that previously
  limped along on misaligned data will now halt with a diagnosable error —
  investigate the worker logs (see below) to find out why the block went
  missing before resuming.

- **F4/F5 — CHANGES RESULTS (per explicit decision).** Two related fixes to
  fuzzy-tracked parameter handling:
  - NA `Rhat` *or* NA `n.eff`, for any checked parameter (focal or
    fuzzy-tracked), is now always fatal. Fuzzy-tracked NA Rhat had briefly
    been a warning in the immediately preceding commit; that is reverted.
    NA `n.eff` was never actually guarded at all before this — arithmetic on
    `NA` silently produced a crash further downstream (`runNimble()`'s
    `if(!mod.check.result)`) rather than a clear diagnosis.
  - Fuzzy-tracked parameters are now governed by `neff.required` as well as
    `Rht.required`: a low `n.eff` counts toward the fuzzy unconverged
    proportion the same way a high `Rhat` does. Previously they escaped the
    `n.eff` requirement entirely.
  - `checkNimble()` also stops with a clear message if `par.fuzzy.track`
    matches no parameters (a typo, or a family absent from a particular
    model), rather than crashing on a `0/0` `NaN` comparison at the first
    convergence check.

- **F11 — CHANGES RESULTS (minor).** The fuzzy Rhat comparison used to round
  to 1 decimal place before comparing against `Rht.required`, so e.g. 1.14
  rounded down to pass at `Rht.required = 1.1` even though the same value
  fails the (unrounded) strict test. Both now use the same precision;
  fuzzy tolerance is very slightly tighter than before.

- **F6 (documentation only, no code change).** `max.samples.saved` is, and
  remains, a **per-chain** cap (`nrow(out) = nc * max.samples.saved`), not
  "across chains" as the old docs said. Fixing the code would change the size
  of every saved model object and the manual-mode stopping point for every
  existing caller; fixing the documentation does not. `?runNimble` now states
  this correctly, and warns that `neff.required > nc * max.samples.saved`
  makes convergence unreachable.

- **F7.** `mcmc.info$nthin` in manual mode ignored the extra thinning
  `gatherNimble()` applies when `max.samples.saved` forces subsampling,
  understating the realized thinning in the saved metadata. Fixed
  (`nthin = nt * additional.thin.rate`, matching automated mode). The samples
  themselves were never affected — this is metadata only.

- **F8.** `gatherNimble2()` dropped whole burn-in blocks but never trimmed the
  residual burn-in within the first retained block, leaking burn-in draws
  into the second monitor set whenever `nt2 < ni`. Fixed — it now trims the
  same way `gatherNimble()` does, at `nt2`'s granularity.

- **F9.** A chain contributing fewer draws to the second monitor set than the
  others (e.g. one block dump predates `parameters2`) was silently accepted,
  quietly unbalancing the stacked matrix. Now warns with the per-chain counts.

- **F10.** `rtrn.model = TRUE` did `assign("mod", mod.nam, ...)` — assigning
  the *string* `"mod"`, not the model object, discarding the actual result.
  With `sav.model = FALSE` (explicitly permitted) this meant the entire run's
  output was lost. Fixed: `assign(mod.nam, mod, ...)`.

## Latent crashes (undefined variables, unexercised branches)

- **F14.** The worker's "undefined condition reached" branch referenced `cn`,
  the *parent's* loop variable, which doesn't exist in the worker process —
  so instead of reporting the problem, the worker died with
  `object 'cn' not found` and no diagnostic reached anyone. Folded into the
  new error-logging mechanism (see below): the branch now raises `stop()` with
  the worker's own `chn` and the current directive/status values, caught by
  the enclosing `tryCatch` and logged.
- **F15.** `countNimbleBlocks()` dropped `drop = FALSE` on one subscript,
  collapsing a matrix to a plain vector and crashing whenever `nc = 1`. Fixed.
- **F16.** `par.fuzzy.track` matching no parameters crashed on a `0/0` `NaN`
  comparison. Now a clear `stop()` (see F4/F5 above).
- **F17.** `mcmcOutputSubset()` couldn't handle a single-chain (`nc = 1`)
  object at all (`1:(nc-1)` degenerates to `1:0`). Fixed.
- **F18.** `par.keep` matching no columns crashed with an opaque
  `subscript out of bounds`. Now names the unmatched values.
- **F19.** `gatherNimble()` crashed opaquely (`invalid 'type' (list)`) when no
  blocks remained after burn-in. Now a clear `stop()`.
- **F20.** `max.samples.saved = NULL` in manual mode (where it is the sole
  stopping criterion) crashed with `argument is of length zero`, *after*
  workers were already launched. Now validated up front, before anything is
  spawned.
- **F23.** An unparseable block/chain number in a filename produced `NA`,
  which crashed the (now-removed) renumbering loop. `countNimbleBlocks()` now
  validates parsed numbers explicitly and names the offending file. Also
  fixed for free by F22's atomic writes: a `.tmp` file is never visible under
  the `mod_chn` pattern in the first place.
- **F24.** `max.tries = 0` skipped the convergence loop entirely, then crashed
  referencing variables that were never assigned. `max.tries` is now
  validated to be at least 1 (or `NULL`) up front.
- **F28.** `proc` was assigned with `<<-`, writing into the calling user's
  global environment — clobbering any object named `proc` there, and letting
  a second concurrent `runNimble()` call overwrite the first's handle. Now a
  properly scoped local variable.

## Concurrency and worker error visibility

The parent and its parallel workers coordinate entirely through plain text
files with no locking. Before this pass, a dead or errored worker was
invisible: its console output was discarded, nothing checked whether it was
still alive, and the parent could poll forever waiting for a chain that was
never coming back.

- **F21.** `dplyr`, `stringr`, `coda`, `mcmcOutput`, and `R.utils` were listed
  under `Suggests:` with no `require()` calls (for most of them) and no
  `NAMESPACE` imports, so `runNimble()` depended on the calling session
  happening to have them attached — and failed *after* launching `nc`
  workers if not, leaving them sampling forever. Fixed properly: all five are
  now `Imports:` with explicit `NAMESPACE` `importFrom()` declarations, so
  package code resolves them through its own imports regardless of what the
  caller's session has attached. `nimbleHMC` remains `Suggests:` (it's
  genuinely optional) with a `requireNamespace()` availability check before
  the conditional `require()` that a user-supplied `SamplerSourcePath` script
  needs.
- **F22.** `gatherNimble()`/`gatherNimble2()` polled `system2("lsof", ...)` to
  avoid reading a block file still being written — Linux-only, and on
  Windows a hard error rather than the assumed no-op. Removed entirely:
  `runNimbleBlock()` now writes each block via a temp file + atomic rename,
  so a file visible under its final name is always complete. This also
  removes a TOCTOU race the polling never actually closed.
- **F25/F27.** Nothing tracked whether workers were still alive, and their
  console output was discarded (`processx` defaults to `NULL` = discard).
  Now: workers' stdout/stderr/exit status are captured per chain via GNU
  `parallel`'s `--results`/`--joblog` under `dump.path/worker_logs` and
  `dump.path/worker_joblog.txt`; the generated worker script wraps its entire
  body in `tryCatch()`, logging any R-level error to
  `dump.path/chn<n>_error.log` with a timestamp before exiting; and the
  parent checks for these on every polling iteration, stopping immediately
  with a pointer to the relevant file(s) instead of polling indefinitely.
- **F26.** Coordination-file reads/writes (the directive and per-chain status
  files) had no protection against torn reads in the worker (a read hitting a
  file mid-rewrite could see zero or multiple lines and crash on the next
  comparison); the parent had a `try()`/retry loop, the worker had none. Both
  sides now write atomically (temp file + rename) and the worker retries a
  malformed read the same way the parent already did.
- **F29.** The initial "has every chain written a block yet" wait counted raw
  matching filenames, satisfiable by `nc` files from a single fast chain. Now
  counts distinct chains.
- **F12.** Only "too many chains" was ever detected; a chain that silently
  wrote nothing was invisible. Now both directions are checked.
- **F13.** `gatherNimble2()` re-inventoried the dump directory independently
  at the end of the run, so a worker that hadn't yet noticed `STOP` could add
  one more block after the primary `gatherNimble()` snapshot, giving the two
  monitor sets different effective burn-in. The parent now waits (up to 5
  minutes, with a warning if exceeded) for workers to actually exit before
  gathering the second set.
- **F30.** No `on.exit()` handler existed anywhere in `runNimble()`, so an
  error thrown after workers were launched left them running indefinitely,
  holding cores. Now: on any exit (success, error, or interrupt), the
  directive file is set to `STOP`, any still-alive worker process is killed,
  and logs are copied out — all best-effort, so a failure in cleanup itself
  can't mask the original error.
- **Log retention.** Worker logs (`worker_logs/`, `worker_joblog.txt`, any
  `chn<n>_error.log`, `Check_log.csv`, `m.csv`) are now copied to
  `<mod.nam>_logs/` before `dump.path` is (optionally) deleted, regardless of
  whether the run succeeded or failed.

## Resources — memory

The package exists to avoid overloading RAM by dumping samples to disk in
blocks; two of its own gather functions defeated that at the gather step.

- **F31.** `gatherNimble()` loaded every retained block for every chain,
  `rbind`ed them all, and only then applied `max.samples.saved` — so peak
  memory scaled with the *total* amount of data sampled so far, not with the
  requested cap (measured at 6–24x the final object size before this fix).
  Redesigned: the retained row indices are now computed analytically from
  `countNimbleBlocks()`'s block accounting *before* any file is loaded (every
  block has a fixed, known row count), so each block is loaded, immediately
  reduced to the rows it actually contributes, and discarded before the next
  one loads. Verified to produce byte-identical output to the previous
  implementation. Measured: with a fixed cap, peak memory now stays flat as
  total retained data grows (tested from 61 MB to 244 MB retained with no
  increase in peak) — it no longer scales with total retained data at all,
  only with the cap.
- **F32.** `mcmcOutputSubset()` re-materialized the full input object inside
  a per-chain loop (once per chain instead of once total). Hoisted out.
- **F33.** `gatherNimble2()` had no size cap at all. Added `max.rows`,
  applied the same way as `gatherNimble()`'s `max.samples.saved`.

## Packaging

- `DESCRIPTION`: `mcmcOutput`, `nimble`, `processx`, `dplyr`, `stringr`,
  `coda`, `R.utils` moved to `Imports:` (were `Suggests:` with no declared
  runtime dependency mechanism); `nimbleHMC` added to `Suggests:`; `Date:`
  typo fixed (day/month were transposed). `License:` intentionally left
  untouched — that's a decision for the maintainer, not a bug fix.
- `NAMESPACE`: replaced the blanket `exportPattern("^[[:alpha:]]+")` with
  explicit exports of the documented public API (`runNimble`, `checkNimble`,
  `mcmcOutputSubset`, `runNimbleBlock`) plus `importFrom()` for every function
  actually used from Imports. `countNimbleBlocks`, `gatherNimble`,
  `gatherNimble2` are now internal (accessible via `manageNimble:::`, as the
  test suite does) rather than exported — **if any existing analysis script
  calls one of these three directly rather than through `runNimble()`, that
  call will now fail to find the function by its bare name and needs the
  `manageNimble:::` prefix, or the script should be pointed at `runNimble()`
  instead.**
- `man/`: added `gatherNimble2.Rd`; fixed every code/documentation mismatch
  `R CMD check` reported (missing `mod.nam`, stale `nt` default of 10 vs the
  actual default of 1, missing `parameters2`/`nt2`, missing
  `monitors2`/`n.thin2`, missing `directive.file` in `\usage`); documented the
  corrected `par.dontign`/fuzzy-tracking behaviour and the new worker-logging
  outputs; fixed a stale `runNimbleBlock.Rd` return-value description that
  didn't match the actual return value even before this review.
- Minor: `mcmcOutputSubset.R`'s `gsub()` call relied on partial argument
  matching (`replace` for `replacement`); spelled out. Added
  `utils::globalVariables()` for column names referenced via `dplyr`
  non-standard evaluation, which `R CMD check`'s codetools pass otherwise
  reports as undefined globals.

## Known limitation, not fixed here

`man/runNimbleBlock.Rd`, `mcmcOutputSubset.Rd`, and elsewhere: `par.ignore`
(and by extension `par.drop`/`par.keep` in `mcmcOutputSubset()`) match by
unanchored substring, so `par.ignore = "p"` also matches `beta.p[1]` and
`mu.p`. This is documented, intentional behaviour, not a bug — `par.dontign`
(now working correctly, see F2) is the designed rescue for it. Anchoring the
match with `^`/`$` instead would be a real behaviour change for every existing
caller and was deliberately left for a separate decision.

Also not addressed: `mcmcOutputSubset()`'s column reordering when `par.keep`
is supplied (kept columns move to the front, changing summary-table row
order). Not identified as a correctness issue, just a quirk worth knowing
about — see the corresponding test in `test-mcmcOutputSubset.R`.

## Testing

`tests/testthat/` — see `tests/README.md` for what's covered, what needs a
compiled NIMBLE model (`MANAGENIMBLE_TEST_NIMBLE=true`), and what remains
untestable off the Linux server (the actual `parallel` process launch, GNU
`parallel`'s `--results` directory layout, genuine concurrent multi-process
races).
