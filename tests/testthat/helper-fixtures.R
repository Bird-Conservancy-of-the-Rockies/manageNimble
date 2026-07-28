# Shared fixtures for the manageNimble test suite.
#
# CONVENTION USED THROUGHOUT THIS SUITE
# -------------------------------------
# Tests fall into three groups, distinguishable by their test_that() label:
#
#   * Plain behavioural tests - assert what the function is supposed to do.
#     Should stay green forever.
#
#   * "FIXED (Fn)" - regression tests for a specific code review finding.
#     During the review these were written as "CHARACTERISES BUG n" tests
#     asserting the *current, wrong* behaviour, so the suite passed against
#     unfixed code; each was flipped to assert the corrected behaviour once
#     its fix landed. See NEWS.md for the full write-up of every finding.
#
#   * Tests documenting an intentional design characteristic that was
#     reviewed and deliberately NOT changed (e.g. par.ignore's unanchored
#     substring matching) - these exist so a future change to that behaviour
#     is a deliberate decision, not an accident.
#
# Nothing here requires the Linux server, GNU parallel, or lsof.

suppressPackageStartupMessages({
  library(coda)
  library(dplyr)
  library(stringr)
})

# countNimbleBlocks/gatherNimble/gatherNimble2 are internal (not exported) as
# of the NAMESPACE tightening in this review - runNimble() is the documented
# entry point; those three are its implementation details (see NEWS.md). Under
# devtools::load_all() (the normal dev loop) they're visible by bare name
# already; under a real install + library(manageNimble) (what R CMD check's
# test_check() does) only ::: reaches them. This alias keeps every test in the
# suite working under both, without qualifying dozens of call sites.
if (!exists("countNimbleBlocks", mode = "function")) {
  countNimbleBlocks <- manageNimble:::countNimbleBlocks
  gatherNimble <- manageNimble:::gatherNimble
  gatherNimble2 <- manageNimble:::gatherNimble2
}

# ---------------------------------------------------------------------------
# Working-directory helper. checkNimble() writes diagnostic CSVs into getwd()
# on its failure paths, so tests that trip those paths must not run in the
# package directory.
# ---------------------------------------------------------------------------
in_tempdir <- function(code) {
  d <- new_dir("wd")
  old <- setwd(d)
  on.exit(setwd(old), add = TRUE)
  force(code)
}

# A fresh, empty directory. Uses tempfile() rather than sample(), because
# make_mcmcOutput() calls set.seed() and would otherwise make sampled directory
# names repeat across test files, silently mixing fixtures together.
new_dir <- function(prefix = "d") {
  d <- tempfile(pattern = prefix)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# Source text of a function body, whitespace-normalised, for the static checks
# on runNimble(). deparse() wraps long lines, so patterns must be matched
# against a squished single string.
fn_text <- function(f) {
  gsub("[[:space:]]+", " ", paste(deparse(body(f), width.cutoff = 500L), collapse = " "))
}

# ---------------------------------------------------------------------------
# Block-dump fixtures
# ---------------------------------------------------------------------------

# Empty files named mod_chn<c>_<b>.RData. Enough for countNimbleBlocks(),
# which only ever looks at file names.
# A chain may be given integer(0) to mean "this chain wrote nothing". Note the
# explicit length check: paste0("mod_chn2_", integer(0), ".RData") does NOT
# return character(0), it returns "mod_chn2_.RData", because paste treats
# zero-length arguments as "".
make_block_names <- function(blocks_per_chain) {
  unlist(lapply(seq_along(blocks_per_chain), function(cn) {
    blks <- blocks_per_chain[[cn]]
    if (length(blks) == 0) return(character(0))
    paste0("mod_chn", cn, "_", blks, ".RData")
  }))
}

make_empty_dump_dir <- function(blocks_per_chain) {
  d <- new_dir("dmp")
  for (f in make_block_names(blocks_per_chain)) file.create(file.path(d, f))
  d
}

# Real dumps containing a `samp` matrix (and optionally `samp2`) whose cells
# record the true global iteration index, so that any mis-alignment between
# chains or mis-applied burn-in is directly detectable in the gathered output.
make_dump_dir <- function(blocks_per_chain, ni.block, nt = 1,
                          with_samp2 = FALSE, nt2 = ni.block) {
  d <- new_dir("dmp")
  for (cn in seq_along(blocks_per_chain)) {
    for (b in blocks_per_chain[[cn]]) {
      iters <- seq((b - 1) * ni.block + nt, b * ni.block, by = nt)
      samp <- cbind(iter = iters,
                    chain = rep(cn, length(iters)),
                    alpha = iters + cn / 10)
      f <- file.path(d, paste0("mod_chn", cn, "_", b, ".RData"))
      if (with_samp2) {
        it2 <- seq((b - 1) * ni.block + nt2, b * ni.block, by = nt2)
        samp2 <- cbind(iter2 = it2, chain2 = rep(cn, length(it2)))
        save(samp, samp2, file = f)
      } else {
        save(samp, file = f)
      }
    }
  }
  d
}

# ---------------------------------------------------------------------------
# Synthetic mcmcOutput objects with per-parameter control of Rhat / n.eff
# ---------------------------------------------------------------------------
#
# kind:
#   "good"    - iid, chains agree            -> Rhat ~ 1,   n.eff large
#   "bad"     - chain-specific mean offset   -> Rhat >> 1.1
#   "const"   - zero variance                -> Rhat NA,    n.eff NA
#   "sticky"  - heavily autocorrelated  -> Rhat ~ 1,   n.eff << 100
#
# "sticky" gives every chain the SAME underlying AR(1) path plus a little
# chain-specific noise. Between-chain variance is therefore near zero (Rhat ~ 1)
# while within-chain autocorrelation stays high (n.eff far below 100). An
# independent AR(1) per chain would not do: with n.eff of a few dozen, Rhat
# fluctuates around 1.1 by chance and the test would be flaky.
make_mcmcOutput <- function(spec, n = 1000, nchain = 3, seed = 1) {
  set.seed(seed)
  ar1 <- function(n, rho = 0.995) {
    x <- numeric(n)
    for (i in 2:n) x[i] <- rho * x[i - 1] + rnorm(1, 0, sqrt(1 - rho^2))
    x
  }
  shared <- lapply(names(spec), function(p) if (spec[[p]] == "sticky") ar1(n) else NULL)
  names(shared) <- names(spec)

  lst <- lapply(seq_len(nchain), function(k) {
    cols <- lapply(names(spec), function(p) {
      switch(spec[[p]],
             good   = rnorm(n),
             bad    = rnorm(n) + k * 50,
             const  = rep(1, n),
             sticky = shared[[p]] + rnorm(n, 0, 0.01),
             stop("unknown kind: ", spec[[p]]))
    })
    m <- do.call(cbind, cols)
    colnames(m) <- names(spec)
    coda::as.mcmc(m)
  })
  mcmcOutput::mcmcOutput(coda::as.mcmc.list(lst))
}

# The character vector runNimble() hands to writeLines() to create
# dump.path/ModRunScript.R, extracted statically from the function body so the
# generated worker script can be inspected without launching any processes.
worker_script_lines <- function() {
  found <- NULL
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("writeLines")) &&
          !is.null(names(x)) && "text" %in% names(x)) {
        found <<- x[["text"]]
      }
      for (i in seq_along(x)) if (!is.null(x[[i]])) try(walk(x[[i]]), silent = TRUE)
    }
  }
  walk(body(manageNimble::runNimble))
  stopifnot(!is.null(found))
  eval(found)
}
