### What the package does

Wrapper functions around NIMBLE for fitting large Bayesian models, per `DESCRIPTION`:
parallel processing, automated convergence checking (Rhat, n.eff), and dumping
posterior samples to disk to avoid overloading RAM.

Used by at least the NFWF_RMR analysis (`03-Analysis.R`), which is about to fit a
Bayesian hierarchical community model with roughly 8 million graph nodes, ~5,000
monitored parameters, and an expected runtime measured in days to weeks. Correctness
matters more than elegance here — a wrong convergence verdict costs weeks.

### Architecture

Understanding this first will save time; it is not obvious from any single file.

- **`runNimble()`** is the entry point. It does *not* run `ni` iterations and stop.
  `ni` is the size of one **block**. The function writes the model objects to
  `dump.path/NimbleObjects.RData`, generates an R script
  (`dump.path/ModRunScript.R`) as a character vector passed to `writeLines`, and
  launches `nc` copies of it via GNU `parallel` through `processx::process$new()`.
- Each worker `source()`s the model file, calls **`runNimbleBlock()`** to build,
  compile and sample one block, saves the samples to
  `dump.path/mod_chn<c>_<i>.RData`, then loops for further blocks.
- **Coordination is via files.** `dump.path/runNimbleDirective.txt` holds
  `GO`/`PAUSE`/`STOP` and is read by every worker. `dump.path/block<c>Status.txt`
  holds per-chain `GO`/`STOP`/`RESUME`. The parent polls these with `readLines` in
  `while` loops with `Sys.sleep`.
- **`countNimbleBlocks()`** inventories the dump files, truncates chains to equal
  length, and drops blocks falling entirely within burn-in.
- **`gatherNimble()`** loads every retained block for every chain, `rbind`s them,
  applies residual burn-in and thinning, and caps the total at
  `max.samples.saved`. It polls `lsof` to avoid reading a file still being written.
- **`checkNimble()`** summarises Rhat and n.eff and returns a pass/fail. The parent
  loops: gather, check, and if not converged write `RESUME` and let the chains add
  another block, up to `max.tries`.
- **`mcmcOutputSubset()`** subsets an `mcmcOutput` object by parameter name.

### Uncommitted / newly committed changes — review these too

Four files were changed very recently and **have never been parsed or executed**.
`git show HEAD` (or `git diff` if not yet committed) isolates them. They should be in
scope for the review, not assumed correct.

- `R/checkNimble.R` — five bug fixes, one of which deliberately changes behaviour
  (see below)
- `R/runNimble.R` — new `parameters2`/`nt2` arguments; two knock-on fixes to fuzzy
  parameter handling
- `R/runNimbleBlock.R` — new `monitors2`/`n.thin2` arguments
- `R/gatherNimble2.R` — new file

**What the new feature does.** NIMBLE's `configureMCMC()` supports a second monitor
set (`monitors2`/`thin2`) saved at its own thinning rate. The package never exposed
it. It's needed for latent state arrays too large to retain at the primary thinning
rate. `nt2` defaults to `ni`, i.e. one retained draw per block per chain. All
defaults are "no second set", so existing calls should be unaffected — **verify that
claim**.

**Bugs already fixed in `checkNimble()`** (don't re-report; do check the fixes):

1. `s.focal` was referenced inside the `par.fuzzy.track` block but only defined
   inside `if(length(par.ignore) > 0)`. Using `par.fuzzy.track` without `par.ignore`
   was a hard error.
2. The fuzzy tolerance compared a proportion against a count
   (`length(Rht.fuzzy) * fuzzy.threshold`), so it could never fire.
3. Diagnostic CSVs filtered on `pfuz`, which after its loop held only the last
   tracked family.
4. `str_c(species, "_sum_at_fail.csv")` referenced an undefined `species`, in an
   unreachable branch. Branch removed.
5. **Behaviour change:** fuzzy-tracked parameters are now excluded from the strict
   Rhat/n.eff test, so `fuzzy.threshold` actually governs them. Previously they were
   subject to the strict test *as well*, which made the fuzzy block provably dead
   code. This makes convergence easier to achieve than before.

**Judgment call worth a second opinion:** NA `Rhat` among fuzzy-tracked parameters
was fatal; it is now a warning, counted toward the unconverged proportion. Rationale:
the proportion formula explicitly adds `sum(is.na(Rht.fuzzy))` to its numerator,
which was unreachable while any NA was fatal. Argue the other side if you think it's
wrong — silently tolerating unsampled parameters is a real risk.

### Known-suspicious areas not yet chased

Starting points, not an exhaustive list.

- **`mcmcOutputSubset()` interaction at the top of `checkNimble()`.** When both
  `par.ignore` and `par.fuzzy.track` are supplied, the call is
  `par.keep = c(par.dontign, par.fuzzy.track)`, `par.drop = par.ignore`. If
  `par.dontign` is empty this appears to keep *only* the fuzzy parameters and discard
  everything else. Is that the intent? What is `par.dontign` actually for? It is
  undocumented and used in exactly one place.
- **`# This doesn't seem to be working. Getting -Inf for Rhat with function
  continuing to run.`** — author's own comment, still unresolved, in `checkNimble()`.
  Root-cause it.
- **`countNimbleBlocks()` renumbers blocks** with the comment `It seems that in some
  cases block numbers are getting skipped`. That is a symptom being papered over.
  Why are blocks skipped? Renumbering could silently misalign chains or burn-in.
- **Burn-in semantics.** `nb < 1` is proportional, `nb >= 1` absolute.
  `countNimbleBlocks()` recomputes `burnin <- ni.block * nblks * burnin`, so as the
  run grows the proportional burn-in grows with it and blocks previously retained get
  discarded. Defensible, but confirm it's intended and that `burnin.realized` /
  `burnin.needed` stay consistent.
- **`gatherNimble()` loads everything into memory.** Every retained block for every
  chain is `rbind`ed before thinning. `DESCRIPTION` claims the package dumps to disk
  "to avoid overloading RAM". Does the gather step defeat that at scale?
- **File-based IPC without locking.** `readLines`/`writeLines` on the directive and
  status files from concurrent processes. A worker can read a partially written file.
  There is a `try()`/retry loop around status reads in one place but not everywhere.
- **`lsof` polling** in `gatherNimble()` (and now `gatherNimble2()`) to detect files
  still being written. Linux-only, and it's worth checking whether it's actually
  reliable for this purpose or just usually fast enough.
- **`proc <<- process$new(...)`** assigns into the global environment with `<<-`.
  Side effect on the user's workspace; also means a second concurrent `runNimble()`
  call would clobber it.
- **Error handling in workers.** If a worker dies — bad inits, compilation failure,
  OOM — does the parent notice, or does it poll forever? The
  `while(sum(str_detect(list.files(dump.path), "mod_chn")) < nc)` loop at the start
  has no timeout.
- **`delete.blocks = TRUE`** unlinks `dump.path` at the end. Anything not gathered
  before that point is lost. Check every exit path, including error paths.
- **Packaging.** `DESCRIPTION` lists nimble, dplyr, stringr, coda, mcmcOutput,
  processx under `Suggests:`, but the code uses `%>%`, `filter()`, `str_split()` etc.
  directly and calls `require()` inside functions. These are hard dependencies and
  belong in `Imports:` with `importFrom` in `NAMESPACE`.
  `NAMESPACE` is `exportPattern("^[[:alpha:]]+")`, exporting everything including
  internals. `man/` has no `.Rd` for `gatherNimble2` and the existing `.Rd` files
  don't document the new arguments, so `R CMD check` will flag mismatches.
- **No tests.** There is no `tests/` directory. Worth discussing what could be tested
  cheaply — `countNimbleBlocks()` and `mcmcOutputSubset()` are pure functions over
  small inputs and are good candidates.

### Constraints

- Runs on a Linux server; `parallel` and `lsof` are assumed present. Development is
  on Windows, where the parallel path can't be exercised.
- Backward compatibility matters — other analyses depend on this package. Flag any
  proposed change that would alter results or convergence behaviour for an existing
  caller, so it can be taken deliberately rather than by accident.

### What can be tested on Windows

Only two things make the package Linux-only: `process$new(command = "parallel", ...)`
in `runNimble()`, and `system2("lsof", ...)` in the gather functions. Everything else
is portable, so do as much as possible locally.

**Run these regardless:**

- `parse()` on every file in `R/`.
- `R CMD check` / `devtools::check()`. This runs fine on Windows and its `codetools`
  pass reports undefined global variables — which is exactly the class of bug already
  found by hand (`species`, `s.focal`). Expect it also to flag the missing
  `gatherNimble2.Rd`, undocumented new arguments, and the `Suggests:`-should-be-
  `Imports:` problem.

**Unit-testable on Windows with no server:**

- `countNimbleBlocks()` — pure function over `list.files()`. Create a temp directory
  containing empty files named `mod_chn1_1.RData`, `mod_chn1_2.RData`,
  `mod_chn2_1.RData` and so on, and assert the chain truncation, block renumbering
  and burn-in dropping. The renumbering path (the `block numbers are getting skipped`
  comment) is worth targeted tests with deliberately gapped filenames.
- `checkNimble()` — takes an `mcmcOutput` and returns pass/fail. Construct small
  synthetic `mcmcOutput` objects with known Rhat and n.eff patterns and assert the
  verdict. **This is where the known bugs were**, so it is the highest-value target.
  Cover at minimum: `par.ignore` empty vs not, `par.fuzzy.track` empty vs not, both
  together, fuzzy parameters above and below `fuzzy.threshold`, and NA Rhat among
  fuzzy parameters.
- `mcmcOutputSubset()` — pure function; test the `par.keep`/`par.drop` interaction
  directly, which should settle the open question about `par.dontign`.
- `runNimbleBlock()` — no `parallel`, no `lsof`. Runs on Windows. Fit a toy NIMBLE
  model (a few nodes, a handful of iterations) and check that the dump file is
  written and, with `monitors2` supplied, contains both `samp` and `samp2` with the
  expected dimensions and thinning. **This exercises essentially all of the new
  feature.**
- `gatherNimble2()` — hand-write two or three fake block dumps containing `samp` and
  `samp2` and check the concatenation and burn-in handling. Note that
  `system2("lsof", ...)` on Windows warns and returns nothing, so the wait loop is
  skipped; confirm that rather than assume it.
- **The generated worker script.** `runNimble()` builds `ModRunScript.R` as a
  character vector passed to `writeLines`. You can't execute it on Windows, but you
  can generate it and `parse()` the result. This validates the multi-line
  `runNimbleBlock(...)` calls that were recently split across `writeLines` elements —
  a real risk, since a misplaced comma there fails only at runtime inside a spawned
  process, where the error is easy to miss.

**Needs the Linux server:**

- `runNimble()` end to end: the `parallel` launch, the file-based `GO`/`PAUSE`/`STOP`
  coordination, the convergence loop, worker failure handling.

Note that workers are separate `Rscript` processes that call `require(manageNimble)`,
so `devtools::load_all()` will not reach them — the server test needs the package
actually installed. Install it to a private library and point `R_LIBS_USER` at it for
that session; `processx` passes the parent environment to the workers, so they will
find it, and the shared library stays untouched for anything else running.
