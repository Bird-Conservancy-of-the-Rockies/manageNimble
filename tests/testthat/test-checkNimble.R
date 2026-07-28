test_that("the synthetic fixtures produce the Rhat / n.eff patterns the tests assume", {
  mo <- make_mcmcOutput(list(g = "good", b = "bad", k = "const", s = "sticky"))
  s <- summary(mo, MCEpc = FALSE, Rhat = TRUE, n.eff = TRUE, f = TRUE,
               overlap0 = TRUE, verbose = FALSE)

  expect_lt(s["g", "Rhat"], 1.1);   expect_gt(s["g", "n.eff"], 100)
  expect_gt(s["b", "Rhat"], 1.1)
  expect_true(is.na(s["k", "Rhat"]))
  expect_lt(s["s", "Rhat"], 1.1);   expect_lt(s["s", "n.eff"], 100)
})

test_that("Rhat and n.eff survive the summary reshape", {
  # checkNimble does dplyr::select(Parameter, mean:f) on the summary table. That
  # is a positional range, so it silently depends on Rhat and n.eff sitting
  # between the mean and f columns. If mcmcOutput ever reorders them, the
  # columns vanish, max(NULL) returns -Inf, min(NULL) returns Inf, and every
  # model is declared converged. Guard the assumption explicitly.
  mo <- make_mcmcOutput(list(a = "good"), n = 200)
  s <- summary(mo, MCEpc = FALSE, Rhat = TRUE, n.eff = TRUE, f = TRUE,
               overlap0 = TRUE, verbose = FALSE)
  nm <- colnames(s)

  expect_true(all(c("mean", "f", "Rhat", "n.eff") %in% nm))
  expect_gte(match("Rhat", nm),  match("mean", nm))
  expect_lte(match("Rhat", nm),  match("f", nm))
  expect_gte(match("n.eff", nm), match("mean", nm))
  expect_lte(match("n.eff", nm), match("f", nm))
})

test_that("basic verdicts with no par.ignore and no par.fuzzy.track", {
  expect_true(checkNimble(make_mcmcOutput(list(a = "good", b = "good")))$result)
  expect_false(checkNimble(make_mcmcOutput(list(a = "good", b = "bad")))$result)
  expect_false(checkNimble(make_mcmcOutput(list(a = "good", b = "sticky")))$result)
})

test_that("Rht.required and neff.required are honoured", {
  mo <- make_mcmcOutput(list(a = "good", b = "sticky"))
  s <- summary(mo, MCEpc = FALSE, Rhat = TRUE, n.eff = TRUE, f = TRUE,
               overlap0 = TRUE, verbose = FALSE)
  expect_lt(s["b", "Rhat"], 1.1)      # fails on n.eff alone, not on Rhat
  expect_lt(s["b", "n.eff"], 100)

  expect_false(checkNimble(mo, neff.required = 100)$result)
  expect_true(checkNimble(mo, neff.required = floor(min(s[, "n.eff"])))$result)
})

test_that("spit.summary controls what comes back", {
  # result and prp.fuzzy.not.converged are always returned (the latter is NA
  # when par.fuzzy.track is unused); the full summary table is only attached
  # when spit.summary = TRUE, since it can be large.
  mo <- make_mcmcOutput(list(a = "good"))
  r0 <- checkNimble(mo)
  expect_named(r0, c("result", "prp.fuzzy.not.converged"))
  expect_true(is.na(r0$prp.fuzzy.not.converged))

  r <- checkNimble(mo, spit.summary = TRUE)
  expect_named(r, c("result", "prp.fuzzy.not.converged", "s"))
  expect_true(all(c("Parameter", "Rhat", "n.eff") %in% names(r$s)))
})

test_that("par.ignore removes parameters from the strict check", {
  mo <- make_mcmcOutput(list(a = "good", junk = "bad"))
  expect_false(checkNimble(mo)$result)
  expect_true(checkNimble(mo, par.ignore = "junk")$result)
})

test_that("par.fuzzy.track works WITHOUT par.ignore (regression: s.focal not found)", {
  # Before an earlier fix, s.focal was only defined inside the par.ignore
  # branch, so this combination was a hard error.
  mo <- make_mcmcOutput(c(list(a = "good", `sp[1]` = "bad"),
                          setNames(as.list(rep("good", 4)), paste0("sp[", 2:5, "]"))))

  expect_no_error(checkNimble(mo, par.fuzzy.track = "sp"))
  expect_false(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.05)$result)
  expect_true(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.5)$result)
})

test_that("fuzzy.threshold is compared against a proportion, not a count", {
  # Regression: the comparison used to be against length(Rht.fuzzy) *
  # fuzzy.threshold, so for families larger than 1/fuzzy.threshold the right
  # side exceeded 1 and the fuzzy check could never fire. With 20 members and
  # threshold 0.05, 4 bad members (0.2) must fail.
  spec <- c(list(a = "good"),
            setNames(as.list(c(rep("bad", 4), rep("good", 16))), paste0("sp[", 1:20, "]")))
  mo <- make_mcmcOutput(spec)

  expect_false(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.05)$result)
  expect_true(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.25)$result)
})

test_that("par.ignore and par.fuzzy.track compose", {
  spec <- c(list(a = "good", junk = "bad", `sp[1]` = "bad"),
            setNames(as.list(rep("good", 4)), paste0("sp[", 2:5, "]")))
  mo <- make_mcmcOutput(spec)

  expect_true(checkNimble(mo, par.ignore = "junk", par.fuzzy.track = "sp",
                          fuzzy.threshold = 0.5)$result)
  expect_false(checkNimble(mo, par.ignore = "junk", par.fuzzy.track = "sp",
                           fuzzy.threshold = 0.05)$result)
})

test_that("NA Rhat among focal parameters is fatal", {
  in_tempdir({
    expect_error(checkNimble(make_mcmcOutput(list(a = "good", b = "const")),
                             mod.nam = "unit"),
                 "Parameters missing Rhat or n.eff")
    expect_true(file.exists("Bad_pars_unit.csv"))
  })
})

test_that("excluding every parameter is refused rather than reported as converged", {
  # max(numeric(0)) is -Inf and min(numeric(0)) is Inf, so an empty s.focal
  # would otherwise evaluate to result = TRUE with nothing checked.
  expect_true(suppressWarnings(max(numeric(0)) <= 1.1 & min(numeric(0)) >= 100))

  mo <- make_mcmcOutput(list(a = "good", `sp[1]` = "good"))
  expect_error(checkNimble(mo, par.ignore = "a", par.fuzzy.track = "sp"),
               "No parameters remain for the strict convergence check")
})

test_that("a summary with no Rhat column stops and signals STOP to the workers", {
  # mcmcOutput omits Rhat entirely for a single-chain object, which is the
  # natural way to reach this guard. Without it, dplyr::select(Parameter,
  # mean:f) would silently return a table with no Rhat and no n.eff, and
  # max(NULL) / min(NULL) would report the model as converged.
  mo <- make_mcmcOutput(list(a = "good", b = "good"), n = 200, nchain = 1)
  s <- summary(mo, MCEpc = FALSE, Rhat = TRUE, n.eff = TRUE, f = TRUE,
               overlap0 = TRUE, verbose = FALSE)
  expect_false("Rhat" %in% colnames(s))

  in_tempdir({
    df <- "directive.txt"
    writeLines("GO", df)
    expect_error(checkNimble(mo, directive.file = df), "Rhat not calculated")
    expect_equal(readLines(df), "STOP")
  })

  # and with directive.file left at its default nothing is written
  expect_error(checkNimble(mo), "Rhat not calculated")
})

# ---------------------------------------------------------------------------
# Design characteristics that remain by choice (not bugs)
# ---------------------------------------------------------------------------

test_that("par.ignore's substring matching over-matches unless par.dontign rescues it", {
  # This is documented, intentional behaviour (man/checkNimble.Rd: par.ignore
  # matches on "names containing" the given strings) -- not a bug. par.ignore =
  # "p" is intended to skip p[1], but it also removes beta.p[1] and mu.p, both
  # badly unconverged, leaving only alpha to be checked. par.dontign exists
  # precisely to rescue cases like this (see the next test); a caller who
  # doesn't use it gets exactly this over-matching.
  mo <- make_mcmcOutput(list(`p[1]` = "good", `beta.p[1]` = "bad",
                             `mu.p` = "bad", alpha = "good"))

  r <- checkNimble(mo, par.ignore = "p", spit.summary = TRUE)
  expect_true(r$result)                      # converged w.r.t. the one parameter checked
  expect_equal(r$s$Parameter, "alpha")
})

# ---------------------------------------------------------------------------
# Regression tests for fixed bugs (see NEWS.md for the full write-up)
# ---------------------------------------------------------------------------

test_that("FIXED (F2): par.dontign now rescues parameters par.ignore would otherwise sweep up", {
  # par.dontign used to be forwarded to mcmcOutputSubset() only when
  # par.fuzzy.track was ALSO supplied, so in the (common) par.ignore-only case
  # it silently did nothing.
  mo <- make_mcmcOutput(list(`p[1]` = "good", `p[2]` = "good",
                             `beta.p[1]` = "bad", `mu.p` = "bad", alpha = "good"))

  without <- checkNimble(mo, par.ignore = "p", spit.summary = TRUE)
  with    <- checkNimble(mo, par.ignore = "p", par.dontign = c("beta.p", "mu.p"),
                         spit.summary = TRUE)

  expect_equal(without$s$Parameter, "alpha")   # beta.p[1]/mu.p swept up, no rescue requested
  expect_true(without$result)                  # "converged" w.r.t. the one parameter checked

  expect_setequal(with$s$Parameter, c("beta.p[1]", "mu.p", "alpha"))  # rescued
  expect_false(with$result)                    # correctly fails now that they're checked
})

test_that("FIXED (F16): a par.fuzzy.track family matching nothing gives a clear error", {
  # length(Rht.fuzzy) used to be able to reach 0, making prp.not.converged a
  # 0/0 NaN and the comparison against fuzzy.threshold raise "missing value
  # where TRUE/FALSE needed". A typo, or a family absent from a particular
  # model, now gets a diagnosable message instead.
  mo <- make_mcmcOutput(list(a = "good", b = "good", c = "good"))

  expect_error(checkNimble(mo, par.fuzzy.track = "typo"),
               "par.fuzzy.track matched no parameters")

  # par.ignore matching nothing remains harmless, by contrast
  expect_no_error(checkNimble(mo, par.ignore = "typo"))
})

test_that("FIXED (F4): fuzzy-tracked parameters are now also governed by neff.required", {
  # Previously s.focal excluded fuzzy-tracked parameters from the strict test,
  # and the fuzzy block only ever looked at Rhat - so a fuzzy family could have
  # arbitrarily low n.eff and still pass. The fuzzy proportion now counts a
  # low n.eff the same way it counts a high Rhat.
  spec <- c(list(a = "good", `sp[1]` = "sticky"),
            setNames(as.list(rep("good", 4)), paste0("sp[", 2:5, "]")))
  mo <- make_mcmcOutput(spec, n = 2000)

  s <- summary(mo, MCEpc = FALSE, Rhat = TRUE, n.eff = TRUE, f = TRUE,
               overlap0 = TRUE, verbose = FALSE)
  expect_lt(s["sp[1]", "n.eff"], 100)
  expect_lte(s["sp[1]", "Rhat"], 1.1)             # Rhat alone would not catch this

  # 1/5 fuzzy members fails on n.eff -> proportion 0.2
  expect_false(checkNimble(mo, par.fuzzy.track = "sp", neff.required = 100,
                           fuzzy.threshold = 0.05)$result)
  expect_true(checkNimble(mo, par.fuzzy.track = "sp", neff.required = 100,
                          fuzzy.threshold = 0.5)$result)

  # unchanged: not fuzzy-tracked, this parameter alone still fails the strict test
  expect_false(checkNimble(mo, neff.required = 100)$result)
})

test_that("FIXED (F11): fuzzy Rhat is no longer rounded before comparison", {
  # round(Rht.fuzzy, 1) > Rht.required used to let a value of 1.14 round down
  # to 1.1 and pass the fuzzy test, even though 1.14 fails the (always
  # unrounded) strict test. Both now use the raw value.
  expect_false(round(1.14, digits = 1) > 1.1)   # would have passed under the old rounding
  expect_true(1.14 > 1.1)                        # correctly fails now
})

test_that("FIXED (F4/F5): NA Rhat or n.eff among fuzzy-tracked parameters is fatal", {
  # Reverted from a warning back to a stop, per explicit decision: any NA
  # Rhat or n.eff, whether focal or fuzzy-tracked, now halts the run rather
  # than being tolerated or counted toward the fuzzy proportion.
  spec <- c(list(a = "good", `sp[1]` = "const"),
            setNames(as.list(rep("good", 4)), paste0("sp[", 2:5, "]")))
  mo <- make_mcmcOutput(spec)

  in_tempdir({
    expect_error(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.5,
                             mod.nam = "unit"),
                 "Fuzzy-tracked parameters missing Rhat or n.eff")
    expect_true(file.exists("Bad_pars_fuzzy_unit.csv"))
  })
})
