# manageNimble test suite

Run the whole suite:

```bash
Rscript -e "devtools::load_all('.'); testthat::test_dir('tests/testthat', package='manageNimble', load_package='none')"
```

Or, against a real install (closer to what `R CMD check` does):

```bash
R CMD build .
R CMD check --no-manual manageNimble_*.tar.gz
```

## What each file covers

| File | Covers | Needs |
|---|---|---|
| `helper-fixtures.R` | Shared fixtures: block-dump builders (`make_empty_dump_dir`, `make_dump_dir`), a synthetic `mcmcOutput` with controlled Rhat/n.eff per parameter (`make_mcmcOutput`), and the `manageNimble:::` aliasing for internal functions. Nothing here needs a Linux server or GNU `parallel`. | Nothing extra |
| `test-countNimbleBlocks.R` | Chain truncation, burn-in accounting, block-gap detection (stops rather than renumbers), NA/unparseable filenames. | Nothing extra |
| `test-mcmcOutputSubset.R` | `par.keep`/`par.drop` interaction (including the "union, not restriction" semantics), single-chain handling, unmatched-name errors. | Nothing extra |
| `test-checkNimble.R` | The full convergence-verdict matrix: `par.ignore`/`par.dontign`/`par.fuzzy.track` in every combination, NA-fatal paths, the fuzzy proportion calculation. This is where almost every review finding lived, so it's the highest-value file if you're auditing a future change here. | Nothing extra |
| `test-gatherNimble.R` | Exact retained iteration ranges (every dump cell encodes its true global iteration index, so misalignment is directly detectable), residual burn-in, `max.samples.saved` capping, and that it works with zero reliance on `lsof`. | Nothing extra |
| `test-gatherNimble2.R` | Same, for the second monitor set: concatenation, residual burn-in at `nt2`'s granularity, `max.rows` capping, unequal-chain-contribution warnings. | Nothing extra |
| `test-runNimbleBlock.R` | Compiles a real (tiny) NIMBLE model and runs actual MCMC: `monitors2`/`n.thin2` dimensions and thinning, that a continuation block doesn't duplicate the second monitor set, `SamplerSourcePath` actually being sourced. | A working C++ toolchain (same one NIMBLE itself needs to compile). Takes a few minutes — set `MANAGENIMBLE_TEST_NIMBLE=false` to skip. Skipped automatically under `R CMD check` (`skip_on_cran()`). |
| `test-worker-script.R` | Statically extracts the character vector `runNimble()` hands to `writeLines()` to build the worker script, and checks: it parses; the multi-line `runNimbleBlock()` calls survive the `writeLines` split; every symbol it references resolves to something real; and (the one dynamic test in this file) that the `tryCatch`/error-log mechanism actually catches a forced error and writes a correct `chn<n>_error.log` when run as a real Rscript subprocess. | Nothing extra for the static checks; the dynamic error-log test spawns a real `Rscript` subprocess (a few seconds). |

## What's still untestable without the Linux server

Two things make the orchestration layer Linux-only: `process$new(command = "parallel", ...)` in `runNimble()`, and (before this review) `system2("lsof", ...)` in the gather functions — the latter is now gone entirely (see `NEWS.md`, finding F22). What remains genuinely untestable here:

- The actual `parallel` launch and multi-process coordination end to end — the `GO`/`PAUSE`/`STOP`/`RESUME` protocol driving real, separate `Rscript` processes rather than the static/single-process checks above.
- GNU `parallel`'s exact `--results`/`--joblog` directory layout for this specific version — the `tryCatch`/`chn<n>_error.log` mechanism is verified independently (see `test-worker-script.R`) and doesn't depend on this working as expected, but the raw stdout/stderr capture itself is unverified here.
- Genuine concurrent file races (two real OS processes reading/writing the same coordination file at the same instant) — the atomic-write and retry-on-read logic is exercised for correctness of the *mechanism*, not under real concurrent load.

If you get access to the Linux server this package is meant to run on, that's the order worth testing in: launch a real 2-chain run with `check.freq` set low, deliberately kill one worker mid-run and confirm the parent reports it (rather than hanging), then a full run to convergence.

## Conventions if you're adding to this suite

- Fixtures that write real block dumps (not just empty files) encode the true global iteration index into the data itself (e.g. `samp[, "iter"]`), so a test can assert exactly which iterations survived burn-in/thinning rather than just checking row counts. Prefer this pattern for anything touching burn-in or thinning.
- `checkNimble()` writes diagnostic CSVs (`Bad_pars_*.csv`) into the working directory on its failure paths. Any test that trips one of those paths should run inside `in_tempdir()` (see `helper-fixtures.R`), not the package directory.
- `countNimbleBlocks`, `gatherNimble`, `gatherNimble2` are internal (not exported) — call them by bare name as usual; `helper-fixtures.R` aliases them from `manageNimble:::` when the suite is running against a real install rather than `devtools::load_all()`.
