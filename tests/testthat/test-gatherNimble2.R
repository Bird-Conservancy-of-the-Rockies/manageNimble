# gatherNimble2() no longer shells out to lsof at all (finding F22 - blocks are
# now written via temp-file + atomic rename). These tests don't need lsof on PATH.
#
# gatherNimble2() now derives its kept rows directly from gatherNimble()'s own
# retained.iters selection (see NEWS.md and both functions' own comments),
# rather than gathering/thinning independently, so nearly every test below
# calls the real gatherNimble() first to get a genuine retained.iters to feed
# in - exactly how runNimble() itself threads the two together. Its return
# value is now list(out, iter.key), not a bare matrix.
#
# nt2/base.thin passed to gatherNimble2() below always match the nt/nt2 the
# fixture dumps were written at (make_dump_dir's own nt/nt2 arguments), exactly
# as runNimble() threads its own nt/nt2 through - see test-worker-script.R.

gather_both <- function(d, burnin, ni.block, nt2, base.thin = 1, max.samples.saved = NULL) {
  g <- gatherNimble(d, burnin = burnin, ni.block = ni.block, base.thin = base.thin,
                    max.samples.saved = max.samples.saved)
  o <- gatherNimble2(d, burnin = burnin, ni.block = ni.block, nt2 = nt2,
                     base.thin = base.thin, retained.iters = g$retained.iters)
  list(primary = g, secondary = o)
}

test_that("second monitor set is concatenated across blocks and chains", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  r <- gather_both(d, burnin = 0.5, ni.block = 100, nt2 = 100)
  o <- r$secondary

  expect_true(is.matrix(o$out))
  expect_equal(dim(o$out), c(4L, 2L))                   # blocks 3,4 x 2 chains
  expect_equal(sort(o$out[, "iter2"]), c(300, 300, 400, 400))
  expect_equal(sort(o$out[, "chain2"]), c(1, 1, 2, 2))
  expect_equal(sort(o$iter.key$iter), sort(o$out[, "iter2"]))   # structural key agrees with the fixture's own data
})

test_that("chains are stacked, not kept separate (documented behaviour)", {
  d <- make_dump_dir(list(1:4, 1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 50)

  r <- gather_both(d, burnin = 0.5, ni.block = 100, nt2 = 50)
  o <- r$secondary

  expect_length(dim(o$out), 2)                          # a plain matrix, not [draw, chain, param]
  expect_equal(nrow(o$out), 3 * 2 * 2)                  # 3 chains x 2 blocks x 2 draws
  expect_null(attr(o$out, "nChains"))
  expect_equal(nrow(o$iter.key), nrow(o$out))
  expect_setequal(unique(o$iter.key$chn), 1:3)
})

test_that("save.path writes an rds that round-trips", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)
  p <- file.path(tempdir(), "monitors2.rds")
  on.exit(unlink(p), add = TRUE)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1, max.samples.saved = NULL)
  o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 100, base.thin = 1,
                     retained.iters = g$retained.iters, save.path = p)

  expect_true(file.exists(p))
  expect_equal(readRDS(p), o)
  expect_setequal(names(o), c("out", "iter.key"))
})

test_that("dumps with no samp2 return NULL with a warning", {
  # The expected result when runNimble() ran without parameters2.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = FALSE)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1, max.samples.saved = NULL)
  expect_warning(o <- gatherNimble2(d, burnin = 0.5, ni.block = 100,
                                    base.thin = 1, retained.iters = g$retained.iters),
                 "no second monitor set found")
  expect_null(o)
})

test_that("no blocks left after burn-in returns NULL with a warning", {
  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  # gatherNimble() itself would stop() (not just warn) for this burnin, so
  # there's no real retained.iters to obtain here - gatherNimble2()'s own
  # "no blocks remain" check fires from countNimbleBlocks() alone, before
  # retained.iters is ever consulted, so any placeholder value is fine.
  expect_warning(o <- gatherNimble2(d, burnin = 5000, ni.block = 100, nt2 = 100,
                                    base.thin = 1, retained.iters = integer(0)),
                 "no blocks remain after burn-in")
  expect_null(o)
})

test_that("REMOVED: max.rows is no longer accepted - the cap is now inherited from gatherNimble()'s own max.samples.saved via retained.iters", {
  # See the "PROPERTY" test below for the replacement mechanism and its own
  # regression coverage.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 25)
  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1, max.samples.saved = NULL)

  expect_error(
    gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 25, base.thin = 1,
                 retained.iters = g$retained.iters, max.rows = 4),
    "unused argument"
  )
})

# ---------------------------------------------------------------------------
# Regression tests for fixed bugs (see NEWS.md for the full write-up)
# ---------------------------------------------------------------------------

test_that("FIXED (F8): residual burn-in is now applied to the second set too", {
  # gatherNimble2 used to drop whole burn-in blocks but, unlike gatherNimble,
  # never dropped the remaining burn-in draws inside the first retained block.
  # Now that its kept rows are the intersection with gatherNimble()'s own
  # retained.iters (which already excludes burn-in, including the residual),
  # this is handled correctly by construction rather than via a separate
  # residual-trim computation. burnin = 0.5 * 100 * 5 = 250, realized 200, so
  # base.thin = 1 retains 251:500 for the primary set; nt2 = 25 intersected
  # with that leaves 275, 300, ..., 500 (225 and 250 are both <= 250 and so
  # excluded).
  d <- make_dump_dir(list(1:5, 1:5), ni.block = 100, with_samp2 = TRUE, nt2 = 25)

  r <- gather_both(d, burnin = 0.5, ni.block = 100, nt2 = 25)
  kept <- sort(unique(r$secondary$out[, "iter2"]))

  expect_true(all(kept > 250))
  expect_equal(kept[1], 275)

  # gatherNimble, on the same directory at its own (finer) thinning rate,
  # independently agrees that nothing at or before iteration 250 survives.
  expect_equal(min(r$primary$out[, "iter"]), 251)
})

test_that("FIXED: fractional burnin.needed no longer silently empties every chain", {
  # Reproduces the exact reported scenario: nb = 0.4, ni.block = 2010,
  # nt2 = ni.block, nblks not a multiple of 5 - see gatherNimble()'s own
  # regression test of the same name for the original root cause. Handled here
  # via gatherNimble()'s (already-fixed) retained.iters, intersected exactly.
  d <- make_dump_dir(list(1:7, 1:7), ni.block = 2010, with_samp2 = TRUE, nt2 = 2010)

  r <- gather_both(d, burnin = 0.4, ni.block = 2010, nt2 = 2010)
  o <- r$secondary

  expect_false(is.null(o))
  expect_gt(nrow(o$out), 0)
  expect_equal(nrow(o$out), 10)                      # 2 chains x 5 surviving blocks
  expect_equal(as.vector(table(o$out[, "chain2"])), c(5L, 5L))
  expect_true(all(o$out[, "iter2"] > 5628))          # every retained draw is past the cutoff
})

test_that("residual burn-in trimming drops exactly the draws before the true cutoff (floor > 0)", {
  # residual/nt2 = 120/100 = 1.2 at base.thin = 1: the first retained block
  # samples samp2 at iterations 100 and 200; 100 is before the true cutoff
  # (120) and must be dropped, 200 is past it and must survive.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 250, with_samp2 = TRUE, nt2 = 100)

  r <- gather_both(d, burnin = 120, ni.block = 250, nt2 = 100)
  o <- r$secondary

  expect_false(100 %in% o$out[, "iter2"])
  expect_true(200 %in% o$out[, "iter2"])
  expect_equal(nrow(o$out), 2 * (4 * 2 - 1))          # 2 chains x (4 blocks x 2 draws - 1 dropped)
})

test_that("FIXED (F9): chains contributing unequal numbers of draws now warns", {
  # If one block dump lacks samp2 - an older dump, or a block written before
  # parameters2 was added - that chain contributes fewer draws and the stacked
  # matrix is unbalanced across chains. Previously silent.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  f <- file.path(d, "mod_chn2_4.RData")
  e <- new.env(); load(f, envir = e); samp <- e$samp; save(samp, file = f)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1, max.samples.saved = NULL)
  expect_warning(
    o <- gatherNimble2(d, burnin = 0.5, ni.block = 100, nt2 = 100, base.thin = 1,
                       retained.iters = g$retained.iters),
    "chains contributed unequal numbers of draws"
  )
  expect_equal(nrow(o$out), 3)
  expect_equal(as.vector(table(o$iter.key$chn)), c(2L, 1L))
})

test_that("nc == 1 (a single chain) works the same way", {
  d <- make_dump_dir(list(1:4), ni.block = 100, with_samp2 = TRUE, nt2 = 100)

  r <- gather_both(d, burnin = 0.5, ni.block = 100, nt2 = 100)
  o <- r$secondary

  expect_equal(nrow(o$out), 2)                # blocks 3,4
  expect_equal(sort(o$out[, "iter2"]), c(300, 400))
  expect_equal(unique(o$iter.key$chn), 1)
  joined <- merge(o$iter.key, r$primary$iter.key, by = c("chn", "iter"))
  expect_equal(nrow(joined), nrow(o$iter.key))
})

test_that("nt2 not a positive integer multiple of base.thin is a clear error", {
  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100, with_samp2 = TRUE, nt2 = 30)
  expect_error(
    gatherNimble2(d, burnin = 0, ni.block = 100, nt2 = 30, base.thin = 7,
                 retained.iters = 1:100),
    "positive integer multiple"
  )
  expect_error(
    gatherNimble2(d, burnin = 0, ni.block = 100, nt2 = 0, base.thin = 1,
                 retained.iters = 1:100),
    "positive integer multiple"
  )
})

test_that("PROPERTY: with max.samples.saved-driven additional thinning on the primary set, the second set's row ratio stays exact and every row traces to a primary row", {
  # This is the scenario the original bug hit: many retry blocks (20/chain,
  # 2000 primary iterations/chain) comfortably exceed max.samples.saved (50),
  # so gatherNimble() compresses the primary set via a systematic subsample
  # (additional.thin.rate > 1). Before this fix, gatherNimble2() had no
  # equivalent cap at all and could retain MORE rows than the (now-compressed)
  # primary set - the opposite of what nt2 > base.thin is supposed to
  # guarantee. Confirms both required properties hold here, not just in the
  # uncapped case the other tests above exercise.
  ni.block <- 100; base.thin <- 1; nt2 <- 5; burnin <- 0
  d <- make_dump_dir(list(1:20, 1:20), ni.block = ni.block, nt = base.thin,
                     with_samp2 = TRUE, nt2 = nt2)

  g <- gatherNimble(d, burnin = burnin, ni.block = ni.block, base.thin = base.thin,
                    max.samples.saved = 50)
  expect_gt(g$additional.thin.rate, 1)   # confirms the cap actually engaged

  o <- gatherNimble2(d, burnin = burnin, ni.block = ni.block, nt2 = nt2,
                     base.thin = base.thin, retained.iters = g$retained.iters)

  # Property 1 (fixed ratio, exact count): a primary-retained iteration
  # survives into the second set iff it also happens to be nt2-aligned - not
  # an approximate ratio, an exact, independently-computable prediction.
  expected.iters <- g$retained.iters[g$retained.iters %% nt2 == 0]
  expect_equal(nrow(o$out), 2 * length(expected.iters))                # 2 chains
  expect_equal(as.vector(table(o$iter.key$chn)), rep(length(expected.iters), 2))

  # Property 2 (verifiable per-row correspondence): every retained second-set
  # row's (chn, iter) is a real join key against the primary set's own
  # iter.key, not an assumed positional pairing.
  joined <- merge(o$iter.key, g$iter.key, by = c("chn", "iter"))
  expect_equal(nrow(joined), nrow(o$iter.key))                         # no orphaned second-set rows

  # Cross-check against the fixture's own encoded true-iteration data,
  # independent of the new iter.key mechanism entirely.
  expect_true(all(o$out[, "iter2"] %in% g$out[, "iter"]))
})
