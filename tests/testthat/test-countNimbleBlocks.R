test_that("clean run: chains equal length, proportional burn-in drops whole blocks", {
  d <- make_empty_dump_dir(list(1:4, 1:4))
  r <- countNimbleBlocks(d, burnin = 0.5, ni.block = 1000)

  expect_equal(r$nblks, 4)
  expect_equal(r$burnin, 2000)             # 0.5 * 1000 * 4
  expect_equal(r$burnin.realized, 2000)    # floor(2000/1000) * 1000
  expect_equal(sort(unique(r$m[, "chn"])), c(1L, 2L))
  expect_equal(unname(r$m[, "blk"]), c(3L, 4L, 3L, 4L))
  expect_equal(rownames(r$m),
               c("mod_chn1_3.RData", "mod_chn1_4.RData",
                 "mod_chn2_3.RData", "mod_chn2_4.RData"))
})

test_that("chains are truncated to the length of the shortest", {
  d <- make_empty_dump_dir(list(1:5, 1:3))
  r <- countNimbleBlocks(d, burnin = 0.5, ni.block = 1000)

  expect_equal(r$nblks, 3)
  expect_equal(max(r$m[, "blk"]), 3)
  expect_equal(as.vector(table(r$m[, "chn"])), c(2L, 2L))
  # burn-in recomputed against the truncated length: 0.5 * 1000 * 3 = 1500
  expect_equal(r$burnin, 1500)
  expect_equal(r$burnin.realized, 1000)
})

test_that("proportional burn-in grows as the run grows, retiring earlier blocks", {
  # Same chains, two lengths. Documents the growing-burn-in semantics.
  r4 <- countNimbleBlocks(make_empty_dump_dir(list(1:4, 1:4)), 0.5, 1000)
  r8 <- countNimbleBlocks(make_empty_dump_dir(list(1:8, 1:8)), 0.5, 1000)

  expect_equal(min(r4$m[, "blk"]), 3)
  expect_equal(min(r8$m[, "blk"]), 5)   # block 3, retained at nblks=4, is now burn-in
  expect_equal(r8$burnin, 4000)
})

test_that("burnin.realized never exceeds burnin (gatherNimble relies on this)", {
  for (nblks in 1:10) {
    for (nb in c(0.1, 0.25, 0.5, 0.9)) {
      r <- countNimbleBlocks(make_empty_dump_dir(list(seq_len(nblks), seq_len(nblks))),
                             burnin = nb, ni.block = 1000)
      expect_lte(r$burnin.realized, r$burnin)
    }
  }
})

test_that("absolute burn-in is taken literally and can retain nothing", {
  d <- make_empty_dump_dir(list(1:2, 1:2))
  r <- countNimbleBlocks(d, burnin = 5000, ni.block = 1000)

  expect_equal(r$burnin, 5000)
  expect_equal(nrow(r$m), 0)   # runNimble() polls on this condition
})

test_that("two-digit chain numbers parse correctly", {
  d <- make_empty_dump_dir(rep(list(1:2), 12))
  r <- countNimbleBlocks(d, burnin = 0.5, ni.block = 1000)

  expect_equal(sort(unique(r$m[, "chn"])), 1:12)
})

test_that("non-block files in dump.path are ignored", {
  d <- make_empty_dump_dir(list(1:3, 1:3))
  file.create(file.path(d, "NimbleObjects.RData"))
  file.create(file.path(d, "ModRunScript.R"))
  file.create(file.path(d, "runNimbleDirective.txt"))
  file.create(file.path(d, "m.csv"))
  file.create(file.path(d, "Check_log.csv"))

  r <- countNimbleBlocks(d, burnin = 0.5, ni.block = 1000)
  expect_equal(r$nblks, 3)
  expect_true(all(grepl("^mod_chn", rownames(r$m))))
})

test_that("a chain with no dumps at all is returned without complaint", {
  # countNimbleBlocks only describes what is on disk; it is not its job to know
  # how many chains SHOULD exist. The "did every chain report in" check lives in
  # runNimble() (finding F12), which knows nc and compares against it.
  d <- make_empty_dump_dir(list(1:3, integer(0), 1:3))

  r <- countNimbleBlocks(d, burnin = 0.5, ni.block = 1000)
  expect_equal(sort(unique(r$m[, "chn"])), c(1L, 3L))
})

# ---------------------------------------------------------------------------
# Regression tests for fixed bugs (see NEWS.md for the full write-up)
# ---------------------------------------------------------------------------

test_that("FIXED (F3): gapped block numbers are refused rather than renumbered", {
  # Chain 1 is missing block 4. Renumbering its block 5 down to "block 4" would
  # silently pair it against chain 2's genuine block 4, misaligning the two
  # chains' iteration ranges without any warning. countNimbleBlocks now stops
  # instead, naming the chain and the missing block number.
  d <- make_empty_dump_dir(list(c(1, 2, 3, 5), 1:4))

  expect_error(countNimbleBlocks(d, burnin = 0.5, ni.block = 1000),
               "missing block dump\\(s\\) 4")
})

test_that("FIXED (F15): a single surviving row no longer drops to a vector", {
  # `m <- m[which(m[,2] <= nblks),]` used to lack drop = FALSE, so one chain
  # with one block collapsed the matrix and the next subscript failed.
  d <- make_empty_dump_dir(list(1))
  expect_no_error(countNimbleBlocks(d, burnin = 0.5, ni.block = 1000))

  d2 <- make_empty_dump_dir(list(1, 1))
  expect_no_error(countNimbleBlocks(d2, burnin = 0.5, ni.block = 1000))
})

test_that("FIXED (F23): an unparseable block number is reported, not a cryptic NA crash", {
  # Any file matching "mod_chn" whose chain/block number does not parse as an
  # integer used to put NA into m, and any(NA != ...) in the (now-removed)
  # renumbering loop was neither TRUE nor FALSE.
  d <- make_empty_dump_dir(list(1:3, 1:3))
  file.create(file.path(d, "mod_chn2_.RData"))       # e.g. a truncated write

  expect_error(countNimbleBlocks(d, burnin = 0.5, ni.block = 1000),
               "could not parse a chain or block number.*mod_chn2_\\.RData")
})
