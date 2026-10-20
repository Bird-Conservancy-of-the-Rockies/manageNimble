# Every dump cell records its true global iteration index, so the assertions
# below pin down exactly which iterations survive burn-in and thinning.
#
# Blocks are written via temp file + atomic rename, so a file visible under its
# final name is always complete; gatherNimble() needs no lsof or other external
# tool on the PATH.

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

test_that("retained.iters / iter.key agree with the fixture's own encoded true iteration values", {
  # retained.iters/iter.key are gatherNimble()'s correspondence key
  # (gatherNimble2() selects its rows from retained.iters - see
  # test-gatherNimble2.R). Checked here against the fixture's own
  # independently encoded "iter" data column.
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = NULL)

  expect_equal(sort(g$retained.iters), 201:400)
  expect_equal(nrow(g$iter.key), nrow(g$out))
  expect_setequal(unique(g$iter.key$chn), 1:2)
  # iter.key is row-aligned with as.matrix(out) (chain-major, ascending
  # iteration within chain) - every row's fixture-encoded true iteration must
  # match iter.key's own claim for that same row.
  expect_equal(g$iter.key$iter, unname(g$out[, "iter"]))
  expect_equal(g$iter.key$chn, unname(g$out[, "chain"]))
})

test_that("retained.iters / iter.key stay correct when max.samples.saved forces a systematic subsample", {
  d <- make_dump_dir(list(1:4, 1:4), ni.block = 100)

  g <- gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                    max.samples.saved = 50)

  expect_equal(g$additional.thin.rate, 4)
  expect_equal(sort(unique(g$retained.iters)), sort(unique(g$out[, "iter"])))
  expect_equal(g$iter.key$iter, unname(g$out[, "iter"]))
  expect_equal(g$iter.key$chn, unname(g$out[, "chain"]))
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
  # Runs with an empty PATH, so no external tool (lsof included) is available.
  withr::local_path(character(0), action = "replace")
  expect_equal(Sys.which("lsof"), c(lsof = ""))

  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100)
  expect_no_error(gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                               max.samples.saved = NULL))
})

# ---------------------------------------------------------------------------
# Error paths and fractional residual burn-in
# ---------------------------------------------------------------------------

test_that("a gapped chain is refused upstream, before gatherNimble gathers anything", {
  # Chain 1 is missing block 4. Relabelling its block 5 as block 4 would
  # silently pair it against chain 2's genuine block 4 and compute Rhat across
  # chains covering different parts of the run; countNimbleBlocks() refuses
  # instead (see test-countNimbleBlocks.R), so gatherNimble() never gathers
  # misaligned data.
  d <- make_dump_dir(list(c(1, 2, 3, 5), 1:4), ni.block = 100)

  expect_error(gatherNimble(d, burnin = 0.5, ni.block = 100, base.thin = 1,
                            max.samples.saved = NULL),
               "missing block dump")
})

test_that("no blocks left after burn-in gives a clear error", {
  d <- make_dump_dir(list(1:2, 1:2), ni.block = 100)

  expect_error(gatherNimble(d, burnin = 5000, ni.block = 100, base.thin = 1,
                            max.samples.saved = NULL),
               "no blocks remain after burn-in")
})

test_that("a fractional residual burn-in (in draws) is floored, so gathering succeeds", {
  # (burnin - burnin.realized) / base.thin is generally fractional, and it
  # seeds row.start, the `:` sequence (block.rows) matched against ind.sav.
  # gatherNimble() floors it; left fractional, every block.rows value would be
  # fractional, never match ind.sav's integers, and nothing would be gathered.
  # Uses ni.block = 2010, nb = 0.4, base.thin = 10, at nblks = 1 (the first
  # automated convergence check) and nblks = 7: burnin.block only advances
  # every 2010 iterations, so residual/10 is fractional at every block count
  # except exact multiples of 5.
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

  # nblks = 5 is the control: the residual is exactly 0 there (an exact
  # multiple of base.thin), so burnin.needed is an integer without flooring.
  d5 <- make_dump_dir(list(1:5, 1:5), ni.block = ni.block, nt = nt)
  g5 <- gatherNimble(d5, burnin = nb, ni.block = ni.block, base.thin = nt,
                     max.samples.saved = NULL)
  expect_equal(nrow(g5$out), 2 * 3 * 201)
})

test_that("max.nblks gathers exactly as if the extra blocks had never been written", {
  # Blocks 5-6 stand in for polling overshoot in seeded manual mode.
  for (mss in list(NULL, 50)) {
    g <- gatherNimble(make_dump_dir(list(1:6, 1:6), ni.block = 100), burnin = 0.5,
                      ni.block = 100, base.thin = 1, max.samples.saved = mss, max.nblks = 4)
    ref <- gatherNimble(make_dump_dir(list(1:4, 1:4), ni.block = 100), burnin = 0.5,
                        ni.block = 100, base.thin = 1, max.samples.saved = mss)
    expect_equal(g$nblks, 4)
    expect_equal(g$out[, c("iter", "chain", "alpha")], ref$out[, c("iter", "chain", "alpha")])
    expect_identical(g$retained.iters, ref$retained.iters)
    expect_identical(g$iter.key, ref$iter.key)
    expect_equal(g$additional.thin.rate, ref$additional.thin.rate)
    expect_true(max(g$out[, "iter"]) <= 400)
  }
})
