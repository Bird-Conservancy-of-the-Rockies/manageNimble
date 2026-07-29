# forceFinishNimble() rescues a runNimble() dump.path (live or abandoned) by
# gathering whatever blocks exist under a manually-specified burn-in and
# checking convergence exactly as runNimble() would. Unlike gatherNimble()'s
# own tests (which encode the true iteration index to pin down exact rows
# retained - see make_dump_dir() in helper-fixtures.R), these fixtures use
# make_dump_dir_iid(): real random draws, so a genuine gatherNimble() +
# checkNimble() pass produces a finite, meaningful Rhat/n.eff rather than the
# Inf/NA that an encoded iteration ramp would give checkNimble() on its
# "chain" column. These tests are about forceFinishNimble()'s own plumbing
# (ni/nt recovery, commit/preview branching, the saved object's structure,
# the live-process gate) - gatherNimble()'s own retained-iteration accounting
# is already covered by test-gatherNimble.R and is not re-verified here.

test_that("errors clearly when dump.path does not exist", {
  expect_error(
    forceFinishNimble("/no/such/dump/dir", burnin = 0.5, mod.nam = "m",
                      check.live = FALSE),
    "does not exist"
  )
})

test_that("errors clearly when dump.path is not a runNimble() dump (no NimbleObjects.RData)", {
  d <- new_dir("notadump")
  expect_error(
    forceFinishNimble(d, burnin = 0.5, mod.nam = "m", check.live = FALSE),
    "does not look like a runNimble\\(\\) dump directory"
  )
})

test_that("ni/nt are read from NimbleObjects.RData, not assumed - a mismatch is caught", {
  # Fixture's blocks genuinely have 100 rows each; NimbleObjects.RData claims ni = 999.
  # If forceFinishNimble() silently trusted a caller-supplied ni instead of the dump's
  # own metadata, this mismatch would go undetected.
  d <- make_dump_dir_iid(list(1:2, 1:2), ni.block = 100, nt = 1, kind = "good")
  add_nimble_objects(d, ni = 999, nt = 1)
  expect_error(
    forceFinishNimble(d, burnin = 0.5, mod.nam = "m", plot.result = FALSE,
                      quiet = TRUE, check.live = FALSE),
    "rows, expected|different configuration"
  )
})

test_that("burnin exceeding what's on disk is a clear error, not a cryptic one", {
  d <- make_dump_dir_iid(list(1:2, 1:2), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)
  expect_error(
    forceFinishNimble(d, burnin = 5000, mod.nam = "m", plot.result = FALSE,
                      quiet = TRUE, check.live = FALSE),
    "exceeds the iterations accumulated"
  )
})

test_that("preview mode (commit = FALSE) reports convergence and writes nothing", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)
  in_tempdir({
    out <- forceFinishNimble(d, burnin = 0.5, mod.nam = "preview_good",
                             plot.result = FALSE, quiet = TRUE, check.live = FALSE)
    expect_true(out$result)
    expect_true(is.data.frame(out$summary))
    expect_false(file.exists("preview_good"))
  })
})

test_that("preview mode correctly reports non-convergence for badly-mixed chains", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "bad")
  add_nimble_objects(d, ni = 100, nt = 1)
  out <- forceFinishNimble(d, burnin = 0.5, mod.nam = "preview_bad",
                           plot.result = FALSE, quiet = TRUE, check.live = FALSE)
  expect_false(out$result)
})

test_that("par.ignore reaches checkNimble()", {
  # This fixture monitors a single parameter, "theta", and ignores it - so
  # nothing at all survives mcmcOutputSubset()'s par.drop step inside
  # checkNimble(), which raises its OWN "no columns remain" error before
  # checkNimble() ever builds a summary table to apply its separate
  # "No parameters remain for the strict convergence check" guard against.
  # (That second guard fires only when par.ignore leaves the mcmcOutput
  # object non-empty but the *fuzzy-adjusted* focal set empty - see
  # test-checkNimble.R's "excluding every parameter is refused..." test,
  # which uses two parameters for exactly this reason.) Either message is a
  # clear, non-cryptic error - this test only confirms par.ignore genuinely
  # reaches checkNimble() through forceFinishNimble(), which it does.
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good", par = "theta")
  add_nimble_objects(d, ni = 100, nt = 1)
  expect_error(
    forceFinishNimble(d, burnin = 0.5, mod.nam = "m", par.ignore = "theta",
                      plot.result = FALSE, quiet = TRUE, check.live = FALSE),
    "no columns remain"
  )
})

test_that("commit = TRUE saves the same structure runNimble() itself would produce", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)
  in_tempdir({
    mod <- forceFinishNimble(d, burnin = 0.5, mod.nam = "commit_good", commit = TRUE,
                             plot.result = FALSE, quiet = TRUE, check.live = FALSE)
    expect_true(file.exists("commit_good"))
    expect_setequal(names(mod), c("mcmcOutput", "summary", "mcmc.info"))
    reloaded <- R.utils::loadObject("commit_good")
    expect_equal(names(reloaded), names(mod))
    expect_equal(unname(mod$mcmc.info["nchains"]), 2)
  })
})

test_that("commit = TRUE records a warning when forced convergence still fails", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "bad")
  add_nimble_objects(d, ni = 100, nt = 1)
  in_tempdir({
    mod <- forceFinishNimble(d, burnin = 0.5, mod.nam = "commit_bad", commit = TRUE,
                             plot.result = FALSE, quiet = TRUE, check.live = FALSE)
    expect_false(is.null(mod$warning))
    expect_match(mod$warning, "forceFinishNimble")
  })
})

test_that("sav.model = FALSE returns the object but writes nothing to disk", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)
  in_tempdir({
    mod <- forceFinishNimble(d, burnin = 0.5, mod.nam = "nosave", commit = TRUE,
                             sav.model = FALSE, plot.result = FALSE, quiet = TRUE,
                             check.live = FALSE)
    expect_false(file.exists("nosave"))
    expect_true(is.list(mod))
    expect_setequal(names(mod), c("mcmcOutput", "summary", "mcmc.info"))
  })
})

test_that("delete.blocks = TRUE removes dump.path after a successful commit", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)
  in_tempdir({
    forceFinishNimble(d, burnin = 0.5, mod.nam = "del", commit = TRUE,
                      sav.model = FALSE, delete.blocks = TRUE,
                      plot.result = FALSE, quiet = TRUE, check.live = FALSE)
  })
  expect_false(dir.exists(d))
})

test_that("plot.result = TRUE does not error", {
  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({ grDevices::dev.off(); unlink(pdf_file) }, add = TRUE)
  expect_no_error(
    forceFinishNimble(d, burnin = 0.5, mod.nam = "mplot", plot.result = TRUE,
                      quiet = TRUE, check.live = FALSE)
  )
})

# ---------------------------------------------------------------------------
# FIXED: liveNimbleWorkers()'s pgrep pattern is now regex-escaped
# ---------------------------------------------------------------------------
# pgrep -f matches its pattern as an extended regex, not a literal substring.
# Before this fix, dump.path was passed to it unescaped, so a path containing
# "(", "+", "[", or other metacharacters - all plausible, unexceptional
# directory names - could fail to match a genuinely live process at all,
# silently defeating the commit-safety gate exactly when it matters most.
# pgrep itself isn't available on this platform to exercise end-to-end (see
# the skip_if() guards below and in the live-process section), so these tests
# verify the escaping helper directly, plus a round-trip through R's own
# regex engine standing in for pgrep's (both are POSIX-extended-flavoured for
# these basic metacharacters).

test_that("regex.escape() escapes every extended-regex metacharacter", {
  metachars <- c(".", "^", "$", "*", "+", "?", "(", ")", "[", "]", "{", "}", "|", "\\")
  for (m in metachars) {
    escaped <- regex.escape(m)
    expect_true(grepl(escaped, m, perl = TRUE),
               info = paste("metacharacter:", m))
  }
  # ordinary characters are passed through unchanged
  expect_equal(regex.escape("dump_species1"), "dump_species1")
})

test_that("regex.escape()'d dump.path matches only its own literal command line, not a false positive", {
  # Mirrors the exact paths that broke the unescaped version, verified against
  # a fake command line the way pgrep would see it (path baked into a worker's
  # command-line arguments).
  paths <- c("dump", "dump_species1", "./dump", "dump(v2)", "output/run.1/dump",
             "dump+final", "C:/analysis/dump[bird]", "dump{2}", "a|b")
  for (p in paths) {
    cmdline <- paste0("Rscript ModRunScript.R 1 ", p, "/NimbleObjects.RData")
    expect_true(grepl(regex.escape(p), cmdline, perl = TRUE), info = paste("path:", p))
  }

  # the unescaped metacharacters really would have broken this - confirms the
  # test fixture is genuinely discriminating, not just checking a tautology
  broken <- "dump(v2)"
  cmdline <- paste0("Rscript ModRunScript.R 1 ", broken, "/NimbleObjects.RData")
  expect_false(grepl(broken, cmdline, perl = TRUE, fixed = FALSE))

  # and no false positives: the escaped pattern must not match a DIFFERENT
  # string that merely contains the same characters, unescaped-regex-style
  expect_false(grepl(regex.escape("dump(v2)"), "dumpXv2)", perl = TRUE))
})

# ---------------------------------------------------------------------------
# Live-process safety gate
# ---------------------------------------------------------------------------
# forceFinishNimble() recognizes a still-running run by pgrep-matching
# dump.path against every process's command line (see liveNimbleWorkers() in
# R/forceFinishNimble.R) - it has no knowledge of runNimble()'s actual
# `parallel`/Rscript workers, so a plain background process invoked with
# dump.path somewhere on its command line is a faithful enough stand-in here.

test_that("commit = TRUE is refused while a process referencing dump.path is alive", {
  skip_if(nchar(Sys.which("pgrep")) == 0 || nchar(Sys.which("bash")) == 0,
          "pgrep/bash not available")

  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)

  px <- processx::process$new("bash", args = c("-c", "sleep 30", d))
  on.exit(try(px$kill(), silent = TRUE), add = TRUE)
  Sys.sleep(0.3) # let it actually start before pgrep looks for it

  in_tempdir({
    expect_error(
      forceFinishNimble(d, burnin = 0.5, mod.nam = "m_live", commit = TRUE,
                        plot.result = FALSE, quiet = TRUE),
      "Refusing to commit"
    )
    # check.live = FALSE bypasses the gate even though the process is still alive.
    mod <- forceFinishNimble(d, burnin = 0.5, mod.nam = "m_live2", commit = TRUE,
                             sav.model = FALSE, plot.result = FALSE, quiet = TRUE,
                             check.live = FALSE)
    expect_true(is.list(mod))
  })
})

test_that("commit = TRUE succeeds once the live process is gone", {
  skip_if(nchar(Sys.which("pgrep")) == 0 || nchar(Sys.which("bash")) == 0,
          "pgrep/bash not available")

  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)

  px <- processx::process$new("bash", args = c("-c", "sleep 30", d))
  Sys.sleep(0.3)
  px$kill()
  px$wait()
  Sys.sleep(0.2) # let the OS finish reaping it

  in_tempdir({
    expect_no_error(
      forceFinishNimble(d, burnin = 0.5, mod.nam = "m_dead", commit = TRUE,
                        sav.model = FALSE, plot.result = FALSE, quiet = TRUE)
    )
  })
})

test_that("commit = TRUE detects a live process even when dump.path contains regex metacharacters", {
  # End-to-end regression test for the regex.escape() fix, exercising the real
  # pgrep call (unlike the unit tests above, which stand in for it with R's
  # own regex engine). Only runs where pgrep/bash are actually available.
  skip_if(nchar(Sys.which("pgrep")) == 0 || nchar(Sys.which("bash")) == 0,
          "pgrep/bash not available")

  d0 <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  d <- paste0(d0, "(v2)")  # inject a regex metacharacter into the path itself
  file.rename(d0, d)
  add_nimble_objects(d, ni = 100, nt = 1)

  px <- processx::process$new("bash", args = c("-c", "sleep 30", d))
  on.exit(try(px$kill(), silent = TRUE), add = TRUE)
  Sys.sleep(0.3)

  in_tempdir({
    expect_error(
      forceFinishNimble(d, burnin = 0.5, mod.nam = "m_live_paren", commit = TRUE,
                        plot.result = FALSE, quiet = TRUE),
      "Refusing to commit"
    )
  })
})

test_that("preview (commit = FALSE) never requires the live-process check to pass", {
  skip_if(nchar(Sys.which("pgrep")) == 0 || nchar(Sys.which("bash")) == 0,
          "pgrep/bash not available")

  d <- make_dump_dir_iid(list(1:4, 1:4), ni.block = 100, kind = "good")
  add_nimble_objects(d, ni = 100, nt = 1)

  px <- processx::process$new("bash", args = c("-c", "sleep 30", d))
  on.exit(try(px$kill(), silent = TRUE), add = TRUE)
  Sys.sleep(0.3)

  expect_no_error(
    forceFinishNimble(d, burnin = 0.5, mod.nam = "m_preview_live",
                      plot.result = FALSE, quiet = TRUE) # commit defaults to FALSE
  )
})
