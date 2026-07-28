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
