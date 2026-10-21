# max.samples.saved's >100-in-automated-mode requirement is validated at the
# top of runNimble(), before dump.path/model.path/nimble/GNU parallel are ever
# touched - like nt2's own validation (test-runNimble-nt2-validation.R), this
# can be tested with a genuine, executed call. Every other required argument
# is left unsupplied on purpose: R's lazy argument evaluation means it's never
# forced before this check's stop() fires, and the last two tests below
# confirm that directly (a call that passes this check fails on the next
# missing argument - inits - instead, so this isn't just an unconditional
# error).

skip_if_not_installed("nimble")   # runNimble() attaches it before reaching inits

test_that("max.samples.saved must exceed 100 in automated mode", {
  # mcmcOutput's summary() only calculates Rhat above 100 draws per chain, and
  # max.samples.saved is a hard per-chain cap applied after burn-in and
  # thinning - checkNimble() would otherwise fail opaquely at the first
  # convergence check, well after workers have launched and sampled.
  expect_error(
    runNimble(parameters = "a", check.freq = 1, max.samples.saved = 100, sav.model = TRUE),
    "max.samples.saved \\(100\\) must exceed 100 in automated mode"
  )
  expect_error(
    runNimble(parameters = "a", check.freq = 1, max.samples.saved = 50, sav.model = TRUE),
    "must exceed 100 in automated mode"
  )
})

test_that("the check does not apply in manual mode", {
  # check.freq = NULL (the default) is manual mode, where max.samples.saved is
  # instead validated only to be a positive number (see its own stop() just
  # above this one in runNimble()) - 50 passes that and proceeds further.
  expect_error(
    runNimble(parameters = "a", max.samples.saved = 50, sav.model = TRUE),
    "inits.*missing"
  )
})

test_that("max.samples.saved = NULL in automated mode proceeds past this check", {
  # NULL means "no cap" to gatherNimble() in automated mode, so this check
  # doesn't apply to it.
  expect_error(
    runNimble(parameters = "a", check.freq = 1, max.samples.saved = NULL, sav.model = TRUE),
    "inits.*missing"
  )
})

test_that("a valid max.samples.saved proceeds past this check", {
  expect_error(
    runNimble(parameters = "a", check.freq = 1, max.samples.saved = 101, sav.model = TRUE),
    "inits.*missing"
  )
})
