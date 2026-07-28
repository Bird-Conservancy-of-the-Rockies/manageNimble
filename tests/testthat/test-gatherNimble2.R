# gatherNimble2() no longer shells out to lsof at all (finding F22 - blocks are
# now written via temp-file + atomic rename). These tests don't need lsof on PATH.
#
# nt2 passed to gatherNimble2() below always matches the nt2 the fixture dumps
# were written at (make_dump_dir's own nt2 argument), exactly as runNimble()
# threads its own nt2 through - see test-worker-script.R.

test_that("second monitor set is concatenated across blocks and chains", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 100)

  expect_true(is.matrix(o))
  expect_equal(dim(o), c(4L, 2L))                   # blocks 3,4 x 2 chains
  expect_equal(sort(o[, "iter2"]), c(300, 300, 400, 400))
  expect_equal(sort(o[, "chain2"]), c(1, 1, 2, 2))
})

test_that("chains are stacked, not kept separate (documented behaviour)", {
  d <- make_dump_dir(list(1:4, 1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 50)

  o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 50)

  expect_length(dim(o), 2)                          # a plain matrix, not [draw, chain, param]
  expect_equal(nrow(o), 3 * 2 * 2)                  # 3 chains x 2 blocks x 2 draws
  expect_null(attr(o, "nChains"))
})

test_that("save.path writes an rds that round-trips", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)
  p <- file.path(tempdir(), "monitors2.rds")
  on.exit(unlink(p), add = TRUE)

  o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 100, save.path = p)

  expect_true(file.exists(p))
  expect_equal(readRDS(p), o)
})

test_that("dumps with no samp2 return NULL with a warning", {
  # The expected result when runNimble() ran without parameters2.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = FALSE)

  expect_warning(o <- gatherNimble2(d, burnin = 0.5, ni.block = 100),
                 "no second monitor set found")
  expect_null(o)
})

test_that("no blocks left after burn-in returns NULL with a warning", {
  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  expect_warning(o <- gatherNimble2(d, burnin = 5000, ni.block = 100, nt2 = 100),
                 "no blocks remain after burn-in")
  expect_null(o)
})

test_that("max.rows systematically subsamples the concatenated result", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 25)

  o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 25, max.rows = 4)

  expect_equal(nrow(o), 4)
})

# ---------------------------------------------------------------------------
# Regression tests for fixed bugs (see NEWS.md for the full write-up)
# ---------------------------------------------------------------------------

test_that("FIXED (F8): residual burn-in is now applied to the second set too", {
  # gatherNimble2 used to drop whole burn-in blocks but, unlike gatherNimble,
  # never dropped the remaining burn-in draws inside the first retained block.
  # With nt2 < ni.block this leaked burn-in draws into the second monitor set.
  # burnin = 0.5 * 100 * 5 = 250, realized 200, so 50 iterations / nt2 = 25 ->
  # 2 draws dropped from the front of each chain's concatenated blocks. Block 3
  # (the first retained block) samples samp2 at iterations 225, 250, 275, 300 -
  # dropping the first 2 (225, 250 - both within/at the burn-in cutoff) leaves
  # 275 as the first surviving draw.
  d <- make_dump_dir(list(1:5, 1:5), ni.block = 100, with_samp2 = TRUE, nt2 = 25)

  o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 25)
  kept <- sort(unique(o[, "iter2"]))

  expect_true(all(kept > 250))
  expect_equal(kept[1], 275)

  # gatherNimble, on the same directory at its own (finer) thinning rate,
  # independently agrees that nothing at or before iteration 250 survives
  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = NULL)
  expect_equal(min(g$out[, "iter"]), 251)
})

test_that("FIXED (F9): chains contributing unequal numbers of draws now warns", {
  # If one block dump lacks samp2 - an older dump, or a block written before
  # parameters2 was added - that chain contributes fewer draws and the stacked
  # matrix is unbalanced across chains. Previously silent.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  f <- file.path(d, "mod_chn2_4.RData")
  e <- new.env(); load(f, envir = e); samp <- e$samp; save(samp, file = f)

  expect_warning(o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 100),
                 "chains contributed unequal numbers of draws")
  expect_equal(nrow(o), 3)
  expect_equal(as.vector(table(o[, "chain2"])), c(2L, 1L))
})
