# Every dump cell records its true global iteration index, so the assertions
# below pin down exactly which iterations survive burn-in and thinning.
#
# gatherNimble() no longer shells out to lsof at all (finding F22 - blocks are
# now written via temp-file + atomic rename, so a file visible under its final
# name is always complete). These tests don't need lsof on PATH.

test_that("happy path: correct iterations retained, chains aligned", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = NULL)

  expect_equal(g$nblks, 4)
  expect_equal(g$additional.thin.rate, 1)
  expect_equal(attr(g$out, "nChains"), 2)
  expect_equal(nrow(g$out), 400)

  it <- g$out[, "iter"]; ch <- g$out[, "chain"]
  expect_equal(range(it[ch == 1]), c(201, 400))
  expect_equal(range(it[ch == 2]), c(201, 400))
  expect_equal(sort(unique(it[ch == 1])), 201:400)
  expect_equal(sort(unique(it[ch == 1])), sort(unique(it[ch == 2])))
})

test_that("residual burn-in inside a block is applied on top of whole-block dropping", {
  # 5 blocks of 100, nb = 0.5 -> burnin 250. Whole blocks only reach 200, so
  # 50 further draws must be dropped from the front.
  d <- make_dump_dir(list(1:5, 1:5), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = NULL)

  it <- g$out[, "iter"]; ch <- g$out[, "chain"]
  expect_equal(range(it[ch == 1]), c(251, 500))
  expect_equal(sum(ch == 1), 250)
  expect_equal(sum(ch == 2), 250)
})

test_that("residual burn-in is expressed in draws, not iterations, when thinning", {
  d <- make_dump_dir(list(1:5, 1:5), ni.block = 100, nt = 5)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 5,
                    max.samples.saved = NULL)

  it <- g$out[, "iter"]; ch <- g$out[, "chain"]
  # burnin 250, burnin.realized 200, so 50/5 = 10 draws dropped from the front
  expect_equal(min(it[ch == 1]), 255)
  expect_equal(sum(ch == 1), 50)
})

test_that("max.samples.saved caps each chain and reports the extra thinning", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = 50)

  expect_equal(as.vector(table(g$out[, "chain"])), c(50L, 50L))
  expect_equal(g$additional.thin.rate, 4)
  # systematic subsample, both chains on identical indices
  it1 <- g$out[g$out[, "chain"] == 1, "iter"]
  it2 <- g$out[g$out[, "chain"] == 2, "iter"]
  expect_equal(it1, it2)
  expect_equal(head(it1, 3), c(201, 205, 209))
})

test_that("max.samples.saved above the available length is a no-op", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = 100000)

  expect_equal(g$additional.thin.rate, 1)
  expect_equal(nrow(g$out), 400)
})

test_that("uneven chains are truncated before gathering", {
  d <- make_dump_dir(list(1:6, 1:4), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = NULL)

  it <- g$out[, "iter"]; ch <- g$out[, "chain"]
  expect_equal(range(it[ch == 1]), c(201, 400))
  expect_equal(range(it[ch == 2]), c(201, 400))
  expect_equal(sum(ch == 1), sum(ch == 2))
})

test_that("gatherNimble works with no lsof on PATH at all", {
  # Regression guard for F22: confirms the lsof dependency is genuinely gone,
  # not just stubbed out by a test fixture.
  withr::local_path(character(0), action = "replace")
  expect_equal(Sys.which("lsof"), c(lsof = ""))

  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100)
  expect_no_error(gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                               max.samples.saved = NULL))
})

# ---------------------------------------------------------------------------
# Regression tests for fixed bugs (see NEWS.md for the full write-up)
# ---------------------------------------------------------------------------

test_that("FIXED (F3): a gapped chain is refused upstream, before gatherNimble ever runs", {
  # Chain 1 is missing block 4. Previously countNimbleBlocks() relabelled its
  # block 5 as block 4, silently pairing it against chain 2's genuine block 4
  # and computing Rhat across chains covering different parts of the run.
  # countNimbleBlocks() now refuses to renumber (see test-countNimbleBlocks.R),
  # so gatherNimble() never gets the chance to gather misaligned data.
  d <- make_dump_dir(list(c(1, 2, 3, 5), 1:4), ni.block = 100)

  expect_error(gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                            max.samples.saved = NULL),
               "missing block dump")
})

test_that("FIXED (F19): no blocks left after burn-in gives a clear error", {
  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100)

  expect_error(gatherNimble(d, burnin = 5000, ni.block = 100, base.thin = 1,
                            max.samples.saved = NULL),
               "no blocks remain after burn-in")
})

test_that("FIXED: fractional burnin.needed no longer crashes gatherNimble() (same root cause as gatherNimble2's fix)", {
  # burnin.needed <- (burnin - burnin.realized) / base.thin used to feed
  # straight into row.start <- 1 - burnin.needed, which seeds a fractional
  # `:` sequence (block.rows) accumulated across every block in the chain.
  # ind.sav is always integer, so block.rows %in% ind.sav never matched,
  # `keep` was empty for every block, pieces ended up entirely empty,
  # do.call(rbind, list()) returned NULL, and as.mcmc(NULL) errored - not a
  # rare edge case, but the deterministic outcome whenever the residual
  # burn-in isn't an exact multiple of base.thin. Reproduces the exact
  # reported values: ni.block = 2010, nb = 0.4, base.thin = 10 (this
  # project's actual settings), at nblks = 1 (fires on the very first
  # automated convergence check) and nblks = 7 (a later non-multiple-of-5
  # case). See countNimbleBlocks()'s own accounting for why: burnin.block
  # only advances every 2010 iterations, so burnin.needed = residual/10 is
  # fractional at every block count except exact multiples of 5.
  ni.block <- 2010; nb <- 0.4; nt <- 10

  # nblks = 1: burnin = 804, burnin.realized = 0, residual = 804,
  # burnin.needed = 80.4 -> floor 80. rows.per.block = 201, so
  # chain.length.now = 201 - 80 = 121 draws/chain, all strictly after the
  # true cutoff (804), starting at the first multiple of 10 above it (810).
  d1 <- make_dump_dir(list(1, 1), ni.block = ni.block, nt = nt)
  g1 <- gatherNimble(d1, burnin = nb, ni.block = ni.block, base.thin = nt,
                     max.samples.saved = NULL)
  expect_false(is.null(g1$out))
  expect_equal(nrow(g1$out), 2 * 121)
  it1 <- g1$out[, "iter"]
  expect_true(all(it1 > 804))
  expect_equal(min(it1), 810)

  # nblks = 7: burnin.block = floor(5628/2010) = 2, so blocks 1-2 are dropped
  # whole and blocks 3-7 (5 retained blocks) survive; burnin.needed =
  # 1608/10 = 160.8 -> floor 160. chain.length.now = 5*201 - 160 = 845.
  d7 <- make_dump_dir(list(1:7, 1:7), ni.block = ni.block, nt = nt)
  g7 <- gatherNimble(d7, burnin = nb, ni.block = ni.block, base.thin = nt,
                     max.samples.saved = NULL)
  expect_false(is.null(g7$out))
  expect_equal(nrow(g7$out), 2 * 845)
  it7 <- g7$out[, "iter"]
  expect_true(all(it7 > 5628))
  expect_equal(min(it7), 5630)
  expect_equal(g7$nblks, 7)

  # nblks = 5 stays the control case: residual is exactly 0 there (an exact
  # multiple of base.thin), so burnin.needed was already an integer before
  # this fix and this case is unaffected by it.
  d5 <- make_dump_dir(list(1:5, 1:5), ni.block = ni.block, nt = nt)
  g5 <- gatherNimble(d5, burnin = nb, ni.block = ni.block, base.thin = nt,
                     max.samples.saved = NULL)
  expect_equal(nrow(g5$out), 2 * 3 * 201)
})
