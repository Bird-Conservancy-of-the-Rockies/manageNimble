# Plan: chain-specific inits and reproducible seeding in `runNimble()`

Status: **implemented 2026-10-20, uncommitted on `main`** (version 0.4.0). Full
suite 739 pass / 0 fail; `R CMD check` 0 errors, the 2 pre-existing warnings,
0 notes. Not verified by execution: `runNimble()`'s own parent loop (needs GNU
`parallel`), so E1's wait and E2's cap are checked statically - first thing to
exercise on the Linux server.
Target version: 0.4.0 (non-breaking, one flagged exception - see Decision 1).
Version on `main` before this work: 0.3.1.

This file is the authority for this work. Re-read it before writing any code
for it; update sections in place when decisions change.

## Goal

1. **Chain-specific initial values.** `inits` accepts, besides today's single
   named list, (a) a per-chain list of `nc` inits lists, or (b) a function
   returning an inits list.
2. **Reproducible sampling.** New `seed = NULL` argument. `NULL` keeps today's
   behaviour exactly (no seeding). A supplied seed makes every chain's draws,
   and the saved model object, reproducible.

Hard constraint: any existing call with a single inits list and no seed must
behave exactly as it does now.

## Current behaviour (confirmed by reading the code, 2026-10-20)

- `runNimble()` saves one `inits` object to `dump.path/NimbleObjects.RData`.
- The generated worker script (`ModRunScript.R`) loads that file and passes
  `inits` unchanged: `runNimbleBlock(mod.lst = list(model, constants, data,
  inits, parameters), ...)` -> `nimbleModel(inits = mod.lst[[4]])`. All chains
  start identically.
- Nothing touches the RNG. Workers are independent `Rscript` processes, so
  chains differ but runs are not reproducible.
- `forceFinishNimble()` no longer exists (removed in `e97c6b6`). No
  resume/restart path: `runNimble()` wipes `dump.path` on entry; a dead worker
  ends the run with an error.

## Verified facts the design rests on

- **NIMBLE's compiled samplers draw from the worker process's R RNG and honour
  `RNGkind("L'Ecuyer-CMRG")`.** Verified 2026-10-20 with four separate
  worker-like `Rscript` processes (toy model; RW + conjugate samplers): same
  seed + same chain -> identical draws, for the first block AND a
  `reset = FALSE` continuation block; different chain or seed -> different;
  block 2 != block 1 (stream continues, not reseeded).
  **`nimbleHMC` NUTS verified too** (2026-10-20, same method, `buildDerivs =
  TRUE`, NUTS on `mu`/`sigma`): same seed -> identical first and continuation
  blocks; different seed -> different.
- **`nimbleModel()` silently ignores unknown inits names, even list-valued
  ones** (verified 2026-10-20).

## Decision 1 - representing and detecting per-chain inits

Applied up front in `runNimble()`, by internal helper `classifyInits()`:

| `inits` is...                                                                 | Treated as     |
|---|---|
| anything that is not a list and not a function (incl. `NULL`), or a list of length 0, or a list whose names are all non-empty | **single** (today's path, passed through unchanged) |
| a list with **no** names (`names()` NULL or all `""`), length > 0            | **per-chain**: must have length `nc`, each element a named list (or empty list) |
| a function                                                                    | **generator**  |
| a list with some names empty and some not                                     | **error**      |

Non-list, non-function values are passed through as "single" (not rejected),
so no existing call that NIMBLE itself accepts starts failing here.

**Guard (decided A: error).** A single (fully named) inits list - and each
per-chain / generated chain's list - containing any list-valued element is an
error. Catches `list(chain1 = ..., chain2 = ...)`, which NIMBLE would accept
and silently ignore. This is the one deliberate exception to "single-list path
unchanged" (a call that silently ran with its inits discarded now errors);
call it out in NEWS.md.

Generator signature: one or more formal arguments -> called `f(chn)` with the
integer chain number; none -> `f()`. Must return a named list (validated).

## Decision 2 - where the inits function is evaluated

**Decided (B): in the parent, before any worker launches**, once per chain
(deviates from the original request, which had it called in the worker).
Reasons: the common idiom `function() list(z = zst, ...)` refers to session
objects that don't exist in a worker; the real call doubles as validation,
before launch; the realized values can be recorded.

The full per-chain list (generated, or supplied per-chain) is saved to
`NimbleObjects.RData` as `inits`, with `inits.per.chain = TRUE`; every worker
loads it and keeps only its own element. Workers never see a function.

With a seed, chain `c`'s call runs inside `nextRNGSubStream(stream[[c]])`.
With `seed = NULL`, calls use the session RNG. Either way the session's RNG
kind and `.Random.seed` are restored afterwards (no side effects on the
caller's session). With `seed = NULL`, inits are reproducible only if the
user `set.seed()`s first; sampling is not.

## Decision 3 - seed derivation and worker seeding

Parent (only when `seed` is non-NULL), helper `deriveChainStreams()`:
1. Save the session's `RNGkind()` and `.Random.seed`; restore both on exit.
2. `RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")` (all three set
   explicitly, so a caller's non-default `normal.kind`/`sample.kind` can't leak
   into the streams); `set.seed(seed)`;
   `stream[[1]] <- nextRNGStream(.Random.seed)`,
   `stream[[c]] <- nextRNGStream(stream[[c - 1]])`.
3. Chain `c`'s inits (if generated) use `nextRNGSubStream(stream[[c]])`; its
   sampling uses `stream[[c]]`.
4. Save `rng.streams` (list of `nc` integer vectors, or `NULL`) to
   `NimbleObjects.RData`.

Worker (generated script): first statements inside the `tryCatch`, before
`source(model.path)` and the first `runNimbleBlock()`:

    if(!is.null(rng.streams)) {
      RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
      assign('.Random.seed', rng.streams[[as.integer(chn)]], envir = globalenv())
    }
    if(inits.per.chain) inits <- inits[[as.integer(chn)]]

Seeded once per worker; continuation blocks continue the stream. The literal
`list(model, constants, data, inits, parameters)` call is unchanged.

Trade-off vs `set.seed(seed + chn)`: simpler and keeps Mersenne-Twister, but
streams collide across runs (seed 1/chain 2 == seed 2/chain 1) and MT seeds
carry no non-overlap guarantee. `nextRNGStream()` streams are guaranteed
non-overlapping (2^127 apart) - base R's standard for parallel work. Cost:
seeded runs use a different generator kind from unseeded ones.

`seed` validation: `NULL`, or a single, non-NA, finite whole number within
integer range (numeric or integer accepted; coerced to integer).

`parallel` (base R) must be declared: DESCRIPTION `Imports:` and NAMESPACE
`importFrom(parallel, nextRNGStream, nextRNGSubStream)`.

## Decision 4 - what is recorded in the output

Appended after the existing elements of every `mod <- list(...)` in
`runNimble()` (three sites; the `monitors2` re-save reuses `mod`):
- `seed` supplied: `mod$rng <- list(seed = <int>, kind = c("L'Ecuyer-CMRG",
  "Inversion", "Rejection"), streams = <nc x 7 integer matrix; row c = chain
  c's sampling start state>)`.
- **Decided (C):** `mod$inits` = the per-chain inits actually used (list of
  `nc`), whenever `inits` was per-chain or a function. Provenance more than
  reproducibility (a seed already regenerates them), but survives the function
  or its inputs changing, is the only record when `seed = NULL`, and lets one
  chain be re-run by hand. Cost: size, if inits hold large latent arrays.
- Single inits + `seed = NULL`: **no new elements** - existing calls' saved
  structure is unchanged.

## Decision 5 (E) - making the gathered block count deterministic

Root cause of the remaining non-reproducibility: how many blocks exist when
the parent gathers can depend on timing, and that count feeds proportional
burn-in, `chain.length.now` and the `max.samples.saved` subsample.

**E1 - automated mode, absolute burn-in warm-up (applies always).** Today the
`nb > 1` warm-up loop writes `RESUME`, sleeps 60 s and re-counts without
waiting for chains to finish; repeated `RESUME`s can give chains different
targets, and the first check gathers mid-run. Fix: after writing `RESUME`,
wait (polling, with `check_workers_alive()`) until every chain reports
`STOP`, then re-count. Each `RESUME` then adds exactly `check.freq` blocks to
every chain. A race fix, not a behaviour change in substance - its outcome is
one the current code can already produce - so it applies with or without a
seed. Parent gains a `read_state()` (retrying read, mirroring the worker's)
and a `wait_for_all_stop()` helper.

**E2 - manual mode (applies only when `seed` is supplied).** Compute up front
`nblks.cap` = the smallest block count satisfying the existing stopping rule
`(ni * n - nb.now(n)) >= max.samples.saved * nt` (helper `minGatherBlocks()`,
using the same expression as the loop). In the polling loop, once the shortest
chain has at least `nblks.cap` blocks, treat `nblks` as `nblks.cap`; gather
with `max.nblks = nblks.cap` so overshoot blocks are ignored; `mcmc.info`
uses `nblks.cap`. Without a seed, `nblks.cap` is `NULL` and the loop runs
exactly as today (applying it unconditionally would change existing results).

Plumbing: new optional `max.nblks = NULL` argument on `countNimbleBlocks()`
(caps `nblks` before trimming and burn-in), passed through by `gatherNimble()`
and `gatherNimble2()`. `NULL` = today's behaviour.

## Reproducibility scope (to document in `?runNimble` / README)

With a seed: identical saved model object given identical model/data/
constants/inits, `nc`, `ni`, `nb`, `nt`, `nt2`, `parameters`/`parameters2`,
`SamplerSourcePath`, `check.freq`, `max.tries`, `max.samples.saved`,
`Rht.required`, `neff.required`, `par.*`, and the same R / NIMBLE / compiler /
platform. Both modes (after E1/E2). Not guaranteed across machines or package
versions. Without a seed: unchanged (not reproducible).

## Backward compatibility (explicit)

Single inits + `seed = NULL`: identical values reach `nimbleModel()`; no RNG
calls added; saved model object unchanged. Internal-only differences:
`ModRunScript.R` gains the lines above; `NimbleObjects.RData` gains
`inits.per.chain` (FALSE) and `rng.streams` (NULL). Two behaviour changes,
both flagged in NEWS.md: the list-valued-inits guard (Decision 1), and E1's
synchronised warm-up (automated mode with `check.freq <= floor(nb / ni)`).

## Files to change

- `R/reproducibility.R` (new, internal, not exported): `classifyInits()`,
  `deriveChainStreams()`, `prepareChainInits()`, `minGatherBlocks()`, and an
  RNG save/restore helper.
- `R/runNimble.R`: `seed` argument (last in signature); validation + parent
  prep after existing checks, after `require(nimble)`, before `dump.path` is
  wiped (keeps `test-runNimble-nt2-validation.R`'s lazy-evaluation tests
  valid); `save(list = ...)` adds `inits.per.chain`, `rng.streams`; worker
  lines; provenance elements in all `mod <- list(...)` sites; E1; E2.
- `R/countNimbleBlocks.R`, `R/gatherNimble.R`, `R/gatherNimble2.R`:
  `max.nblks = NULL` pass-through.
- `R/runNimbleBlock.R`: no change.
- `DESCRIPTION` (0.4.0, date, `parallel` in Imports), `NAMESPACE`.
- `man/runNimble.Rd` (`inits`, `seed`, details on reproducibility, `\value`
  `rng`/`inits`); `man/countNimbleBlocks.Rd`, `man/gatherNimble.Rd`,
  `man/gatherNimble2.Rd` (`max.nblks`).
- `README.md`: "How it works" subsection on inits/seeding/reproducibility;
  **(D)** remove `forceFinishNimble()` from Acknowledgement and Package
  contents; correct "The second monitor set" (`list(out, iter.key)`, derived
  from primary set's retained iterations, joinable via `mod$iter.key`).
- `NEWS.md`: 0.4.0 entry.
- `tests/README.md`: new test files.
- `.Rbuildignore`: excludes this plan file (done).

## Tests

- `test-reproducibility-helpers.R` (new, pure R): every `inits` form and edge
  case; generator called once per chain, with/without chain argument; return
  validation; guard; stream derivation deterministic/distinct/seed-dependent;
  session RNG kind and `.Random.seed` restored exactly (incl. when no
  `.Random.seed` existed); seeded generator inits reproducible;
  `minGatherBlocks()` agrees with the loop's stopping expression.
- `test-runNimble-inits-seed-validation.R` (new): genuine `runNimble()` calls
  failing in up-front validation (bad seeds, bad inits forms, generator
  returning a non-list) - and no `dump.path` created.
- `test-countNimbleBlocks.R` / `test-gatherNimble.R` / `test-gatherNimble2.R`:
  `max.nblks` caps blocks; `NULL` unchanged.
- `test-worker-script.R` (update): saved-symbol list incl. new objects; seed
  install and inits selection inside `tryCatch` and before
  `source(model.path)` / first `runNimbleBlock()`; E1 structure (warm-up waits
  for all `STOP`, no bare `Sys.sleep(60)` re-count); E2 structure.
- Reproducibility integration test (new, gated by `MANAGENIMBLE_TEST_NIMBLE`):
  runs the **real generated `ModRunScript.R`** in `processx` subprocesses
  (2 chains, automated mode, `check.freq = 2` -> 2 blocks incl. a continuation
  block), with `NimbleObjects.RData` built through the real helpers. The test
  model file saves the `inits` each worker actually selected. Same seed ->
  identical draws and inits; different seed -> different; chain 1 != chain 2;
  per-chain inits reach the right chain.
