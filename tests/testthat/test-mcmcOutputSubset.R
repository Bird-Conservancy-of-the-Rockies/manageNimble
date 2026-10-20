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
# Intentional design characteristics (see test-checkNimble.R for the
# par.dontign rescue) and edge cases
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
})

test_that("par.keep is over-matched the same way, by the same design", {
  mo <- make_mcmcOutput(list(`p[1]` = "good", `beta.p[1]` = "good",
                             `mu.p` = "good", alpha = "good"), n = 100)
  expect_setequal(cols(mcmcOutputSubset(mo, par.keep = "p")),
                  c("p[1]", "beta.p[1]", "mu.p"))
})

test_that("a single-chain object can be subset", {
  # The chain-boundary index matrix is built from 1:(nc-1), which for nc = 1 is
  # 1:0 = c(1, 0), a descending range; mcmcOutputSubset() handles nc == 1
  # directly.
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100, nchain = 1)
  expect_equal(cols(mcmcOutputSubset(mo, par.drop = "a")), "b[1]")

  expect_no_error(mcmcOutputSubset(
    make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100, nchain = 2),
    par.drop = "a"))
})

test_that("par.keep matching nothing gives a diagnosable message", {
  # A typo'd or absent parameter name is named in the error, rather than
  # surfacing as zero columns and an obscure "subscript out of bounds".
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good"), n = 100)
  expect_error(mcmcOutputSubset(mo, par.keep = "typo"), "par.keep matched no columns")

  # an empty par.keep (as opposed to a non-empty one matching nothing) is
  # valid - it's the same as supplying par.drop alone.
  expect_no_error(mcmcOutputSubset(mo, par.keep = character(0), par.drop = "a"))
})

test_that("column order is permuted by par.keep (documented behaviour)", {
  # Kept families are moved to the front, so summary() rows don't follow the
  # input parameter order. Pinned here so a change to this is deliberate.
  mo <- make_mcmcOutput(list(`a[1]` = "good", `b[1]` = "good",
                             `c[1]` = "good", d = "good"), n = 100)
  expect_equal(cols(mcmcOutputSubset(mo, par.keep = "d", par.drop = "a")),
               c("d", "b[1]", "c[1]"))
})
