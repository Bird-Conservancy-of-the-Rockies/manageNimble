cols <- function(x) dimnames(x)[[2]]

test_that("par.drop removes exactly the named families", {
  mo <- make_mcmcOutput(list(`a[1]` = "good", `a[2]` = "good",
                             `b[1]` = "good", d = "good"), n = 100)
  expect_equal(cols(mcmcOutputSubset(mo, par.drop = "a")), c("b[1]", "d"))
  expect_equal(cols(mcmcOutputSubset(mo, par.drop = "d")), c("a[1]", "a[2]", "b[1]"))
})

test_that("par.keep alone keeps only the named families", {
  mo <- make_mcmcOutput(list(`a[1]` = "good", `a[2]` = "good",
                             `b[1]` = "good", d = "good"), n = 100)
  expect_equal(cols(mcmcOutputSubset(mo, par.keep = "a")), c("a[1]", "a[2]"))
  expect_equal(cols(mcmcOutputSubset(mo, par.keep = "d")), "d")
})

test_that("par.keep + par.drop keeps everything not dropped, PLUS the keeps", {
  # This settles what par.dontign is for. Supplying par.keep does NOT restrict
  # the result to par.keep: the union of (par.keep) and (everything not in
  # par.drop) is returned. par.dontign therefore acts as a rescue list for
  # parameters that par.drop would otherwise remove.
  mo <- make_mcmcOutput(list(`a[1]` = "good", `a[2]` = "good",
                             `b[1]` = "good", `c[1]` = "good", d = "good"), n = 100)

  expect_setequal(cols(mcmcOutputSubset(mo, par.keep = "b", par.drop = "a")),
                  c("b[1]", "c[1]", "d"))

  # an empty par.keep gives the same answer as par.drop alone
  expect_setequal(cols(mcmcOutputSubset(mo, par.keep = character(0), par.drop = "a")),
                  cols(mcmcOutputSubset(mo, par.drop = "a")))

  # a family in both is rescued, and warns
  expect_warning(r <- mcmcOutputSubset(mo, par.keep = "a", par.drop = "a"),
                 "will be summarized even if they appear")
  expect_setequal(cols(r), c("a[1]", "a[2]", "b[1]", "c[1]", "d"))
})

test_that("both arguments empty is an error", {
  mo <- make_mcmcOutput(list(a = "good"), n = 100)
  expect_error(mcmcOutputSubset(mo), "need to be specified")
})

test_that("chain structure is preserved", {
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100, nchain = 3)
  r <- mcmcOutputSubset(mo, par.drop = "a")

  expect_equal(attr(r, "nChains"), 3)
  expect_equal(nrow(r), nrow(mo))
  # chain 2 of the subset is chain 2 of the original, same rows
  expect_equal(unname(r[101:200, "b[1]"]), unname(mo[101:200, "b[1]"]))
})

# ---------------------------------------------------------------------------
# Design characteristics that remain by choice (not bugs) - see the matching
# section in test-checkNimble.R for the rescue mechanism (par.dontign, F2)
# ---------------------------------------------------------------------------

test_that("par.drop patterns are unanchored by design, and over-match without par.dontign", {
  # "p" is meant to select p[1], p[2]. The grep patterns are "p\\[" and "p$"
  # with no leading ^, so they also match beta.p[1] and mu.p, which are
  # therefore also removed. Documented, intentional (man/checkNimble.Rd:
  # par.ignore matches names "containing" the given strings) - par.dontign is
  # the designed rescue for exactly this, and it works correctly (see
  # checkNimble()'s tests).
  mo <- make_mcmcOutput(list(`p[1]` = "good", `p[2]` = "good",
                             `beta.p[1]` = "good", `mu.p` = "good",
                             `psi[1]` = "good", alpha = "good"), n = 100)

  expect_equal(cols(mcmcOutputSubset(mo, par.drop = "p")), c("psi[1]", "alpha"))

  # what an anchored match would give instead, for contrast:
  # c("beta.p[1]", "mu.p", "psi[1]", "alpha")
})

test_that("par.keep is over-matched the same way, by the same design", {
  mo <- make_mcmcOutput(list(`p[1]` = "good", `beta.p[1]` = "good",
                             `mu.p` = "good", alpha = "good"), n = 100)
  expect_setequal(cols(mcmcOutputSubset(mo, par.keep = "p")),
                  c("p[1]", "beta.p[1]", "mu.p"))
})

test_that("FIXED (F17): a single-chain object can now be subset", {
  # m <- chn.iter * matrix(c(0, rep(1:(nc-1), each = 2), nc), nrow = 2) + c(1,0)
  # For nc = 1, 1:(nc-1) was 1:0 = c(1, 0), producing a 2x3 index matrix with a
  # descending range. nc == 1 is now handled directly.
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100, nchain = 1)
  expect_equal(cols(mcmcOutputSubset(mo, par.drop = "a")), "b[1]")

  expect_no_error(mcmcOutputSubset(
    make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100, nchain = 2),
    par.drop = "a"))
})

test_that("FIXED (F18): par.keep matching nothing gives a diagnosable message", {
  # A typo'd or absent parameter name used to yield zero columns and an
  # obscure "subscript out of bounds", rather than naming the problem.
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100)
  expect_error(mcmcOutputSubset(mo, par.keep = "typo"), "par.keep matched no columns")

  # an empty par.keep (as opposed to a non-empty one matching nothing) is
  # still valid - it's the same as supplying par.drop alone.
  expect_no_error(mcmcOutputSubset(mo, par.keep = character(0), par.drop = "a"))
})

test_that("column order is permuted by par.keep (documented behaviour, not a bug)", {
  # Kept families are moved to the front, so summary() row order no longer
  # matches the original parameter order. Not one of the findings selected for
  # a fix; documented here so a future change to this is deliberate.
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good",
                             `c[1]` = "good", d = "good"), n = 100)
  expect_equal(cols(mcmcOutputSubset(mo, par.keep = "d", par.drop = "a")),
               c("d", "b[1]", "c[1]"))
})
