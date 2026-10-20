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

test_that("par.fuzzy.track works without par.ignore", {
  # s.focal is built whether or not par.ignore is supplied.
  mo <- make_mcmcOutput(c(list(a = "good", `sp[1]` = "bad"),
                          setNames(as.list(rep("good", 4)), paste0("sp[", 2:5, "]"))))

  expect_no_error(checkNimble(mo, par.fuzzy.track = "sp"))
  expect_false(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.05)$result)
  expect_true(checkNimble(mo, par.fuzzy.track = "sp", fuzzy.threshold = 0.5)$result)
})

test_that("fuzzy.threshold is compared against a proportion, not a count", {
  # The unconverged fraction of the family is compared with fuzzy.threshold:
  # with 20 members and threshold 0.05, 4 bad members (0.2) must fail.
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
# Intentional design characteristics
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
# par.dontign rescue, fuzzy tracking and NA handling
# ---------------------------------------------------------------------------

test_that("par.dontign rescues parameters par.ignore would otherwise sweep up", {
  # par.dontign applies whenever par.ignore does, including when
  # par.fuzzy.track isn't supplied.
  mo <- make_mcmcOutput(list(`p[1]` = "good", `p[2]` = "good",
                             `beta.p[1]` = "bad", `mu.p` = "bad", alpha = "good"))

  without <- checkNimble(mo, par.ignore = "p", spit.summary = TRUE)
  with    <- checkNimble(mo, par.ignore = "p", par.dontign = c("beta.p", "mu.p"),
                         spit.summary = TRUE)

  expect_equal(without$s$Parameter, "alpha")   # beta.p[1]/mu.p swept up, no rescue requested
  expect_true(without$result)                  # "converged" w.r.t. the one parameter checked

  expect_setequal(with$s$Parameter, c("beta.p[1]", "mu.p", "alpha"))  # rescued
  expect_false(with$result)                    # fails, since they're checked
})

test_that("par.ignore'd parameters stay out of $s even when a par.fuzzy.track entry substring-matches them", {
  # "sd.pnt_occupancy" is in par.ignore, but "pnt_occupancy" (a different,
  # legitimately fuzzy-tracked parameter) is in par.fuzzy.track and
  # unanchored-matches "sd.pnt_occupancy[1]" as a substring, rescuing it into
  # par.keep. checkNimble() re-applies the exact base-name drop, so it is kept
  # out of both the convergence check and the returned $s.
  mo <- make_mcmcOutput(list(`pnt_occupancy[1]` = "good", `sd.pnt_occupancy[1]` = "good",
                             alpha = "good"))

  r <- checkNimble(mo, par.ignore = "sd.pnt_occupancy",
                   par.fuzzy.track = "pnt_occupancy",
                   fuzzy.threshold = 0.5, spit.summary = TRUE)

  expect_false("sd.pnt_occupancy[1]" %in% r$s$Parameter)
  expect_true("pnt_occupancy[1]" %in% r$s$Parameter)
})

test_that("a par.fuzzy.track family matching nothing gives a clear error", {
  # With no matching members the fuzzy proportion would be 0/0 = NaN, and
  # comparing it with fuzzy.threshold would fail with "missing value where
  # TRUE/FALSE needed". A typo, or a family absent from a particular model,
  # gets a diagnosable message instead.
  mo <- make_mcmcOutput(list(a = "good", b = "good", c = "good"))

  expect_error(checkNimble(mo, par.fuzzy.track = "typo"),
               "par.fuzzy.track matched no parameters")

  # par.ignore matching nothing is harmless, by contrast
  expect_no_error(checkNimble(mo, par.ignore = "typo"))
})

test_that("fuzzy-tracked parameters are governed by neff.required as well as Rht.required", {
  # Fuzzy-tracked parameters are excluded from the strict test, so the fuzzy
  # proportion counts a low n.eff the same way it counts a high Rhat - a fuzzy
  # family can't pass with arbitrarily low n.eff.
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

  # not fuzzy-tracked, the same parameter fails the strict test
  expect_false(checkNimble(mo, neff.required = 100)$result)
})

test_that("NA Rhat or n.eff among fuzzy-tracked parameters is fatal", {
  # Any NA Rhat or n.eff, whether focal or fuzzy-tracked, halts the run rather
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
