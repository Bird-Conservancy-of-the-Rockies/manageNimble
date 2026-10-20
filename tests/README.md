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
| `helper-fixtures.R` | Shared fixtures: block-dump builders (`make_empty_dump_dir`, `make_dump_dir`), a synthetic `mcmcOutput` with controlled Rhat/n.eff per parameter (`make_mcmcOutput`), the generated worker script (`worker_script_lines`), and the `manageNimble:::` aliasing for internal functions. Nothing here needs a Linux server or GNU `parallel`. | Nothing extra |
| `test-countNimbleBlocks.R` | Chain truncation, burn-in accounting, block-gap detection (stops rather than renumbers), NA/unparseable filenames, the `max.nblks` cap. | Nothing extra |
| `test-mcmcOutputSubset.R` | `par.keep`/`par.drop` interaction (including the "union, not restriction" semantics), single-chain handling, unmatched-name errors. | Nothing extra |
| `test-checkNimble.R` | The full convergence-verdict matrix: `par.ignore`/`par.dontign`/`par.fuzzy.track` in every combination, NA-fatal paths, the fuzzy proportion calculation. The convergence verdict is the most consequential logic in the package, so this is the highest-value file if you're auditing a change here. | Nothing extra |
| `test-gatherNimble.R` | Exact retained iteration ranges (every dump cell encodes its true global iteration index, so misalignment is directly detectable), residual burn-in, `max.samples.saved` capping, the `retained.iters`/`iter.key` correspondence key, the `max.nblks` cap, and that it needs no external tools such as `lsof` on the PATH. | Nothing extra |
| `test-gatherNimble2.R` | The second monitor set: kept rows derived from the primary set's `retained.iters` (exact `nt2/nt` ratio and a real `iter.key` join, including under a `max.samples.saved` cap), burn-in, single chain, `nt2` validation, unequal-chain-contribution warnings, the `max.nblks` cap. | Nothing extra |
| `test-reproducibility-helpers.R` | The parent-side helpers behind `inits` and `seed`: every `inits` form and edge case, the list-valued-inits guard, generator calling conventions, per-chain stream derivation, exact restoration of the session's RNG state, `minGatherBlocks()`. | Nothing extra |
| `test-runNimble-nt2-validation.R` / `test-runNimble-inits-seed-validation.R` | Genuine `runNimble()` calls that fail in its up-front validation (`nt2`; `seed` and `inits`), relying on lazy evaluation so nothing past validation is ever reached - plus checks that valid input gets past it. | Nothing extra (the second needs `nimble` installed, since `runNimble()` attaches it before preparing inits) |
| `test-runNimble-monitors2.R` | `consolidate.monitors2` and `mod$iter.key`, mirroring `runNimble()`'s attach-and-save logic against real block dumps, including a real join of the two monitor sets after a save/load round trip. | A working C++ toolchain; set `MANAGENIMBLE_TEST_NIMBLE=false` to skip. |
| `test-runNimble-reproducibility.R` | `seed` and per-chain `inits` end to end: runs the real generated worker script in separate `Rscript` processes (2 chains, 2 blocks each, three runs) and checks identical draws and inits for the same seed, different draws for a different seed or chain, and per-chain inits reaching the right chain. Under `devtools::load_all()` each worker first loads the source tree, so it tests the code in this checkout rather than whatever version is installed. | A working C++ toolchain; compiles 6 models at once (a few minutes). Set `MANAGENIMBLE_TEST_NIMBLE=false` to skip. |
| `test-runNimbleBlock.R` | Compiles a real (tiny) NIMBLE model and runs actual MCMC: `monitors2`/`n.thin2` dimensions and thinning, that a continuation block doesn't duplicate the second monitor set, `SamplerSourcePath` actually being sourced. | A working C++ toolchain (same one NIMBLE itself needs to compile). Takes a few minutes — set `MANAGENIMBLE_TEST_NIMBLE=false` to skip. Skipped automatically under `R CMD check` (`skip_on_cran()`). |
| `test-worker-script.R` | Statically extracts the character vector `runNimble()` hands to `writeLines()` to build the worker script, and checks: it parses; the multi-line `runNimbleBlock()` calls survive the `writeLines` split; every symbol it references resolves to something real; and (the one dynamic test in this file) that the `tryCatch`/error-log mechanism actually catches a forced error and writes a correct `chn<n>_error.log` when run as a real Rscript subprocess. | Nothing extra for the static checks; the dynamic error-log test spawns a real `Rscript` subprocess (a few seconds). |

## What's untestable without the Linux server

`process$new(command = "parallel", ...)` in `runNimble()` makes the orchestration layer Linux-only. What can't be tested here:

- The actual `parallel` launch and the parent's side of the coordination loop inside `runNimble()` (polling, issuing `RESUME`/`PAUSE`/`STOP`, the seeded manual-mode block cap) — these are checked statically. The worker's side is exercised for real by `test-runNimble-reproducibility.R`, which drives real `Rscript` workers through `GO` -> `STOP` via the status and directive files.
- GNU `parallel`'s exact `--results`/`--joblog` directory layout for this specific version — the `tryCatch`/`chn<n>_error.log` mechanism is verified independently (see `test-worker-script.R`) and doesn't depend on this working as expected, but the raw stdout/stderr capture itself is unverified here.
- Genuine concurrent file races (two real OS processes reading/writing the same coordination file at the same instant) — the atomic-write and retry-on-read logic is exercised for correctness of the *mechanism*, not under real concurrent load.

On the Linux server, `smoke-test/runNimble_smoke_test.R` runs `runNimble()` end to end on a toy model (see its header for how to run it against a separate library). Beyond that, worth checking by hand: deliberately kill one worker mid-run and confirm the parent reports it rather than hanging.

## Conventions if you're adding to this suite

- Fixtures that write real block dumps (not just empty files) encode the true global iteration index into the data itself (e.g. `samp[, "iter"]`), so a test can assert exactly which iterations survived burn-in/thinning rather than just checking row counts. Prefer this pattern for anything touching burn-in or thinning.
- `checkNimble()` writes diagnostic CSVs (`Bad_pars_*.csv`) into the working directory on its failure paths. Any test that trips one of those paths should run inside `in_tempdir()` (see `helper-fixtures.R`), not the package directory.
- `countNimbleBlocks`, `gatherNimble`, `gatherNimble2` and the helpers in `R/reproducibility.R` are internal (not exported) — call them by bare name as usual; `helper-fixtures.R` aliases them from `manageNimble:::` when the suite is running against a real install rather than `devtools::load_all()`.
- Tests that change the session's RNG kind or `.Random.seed` must restore both (see `local_session_rng()` in `test-reproducibility-helpers.R`), or later tests that rely on `set.seed()` can silently behave differently.
