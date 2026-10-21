# `seed` and `inits` are validated - and a generator `inits` function actually
# run - inside runNimble() itself, before dump.path is touched or any worker
# launches. Like test-runNimble-nt2-validation.R, these are genuine runNimble()
# calls rather than static checks: every argument the function never reaches
# is left unsupplied (R's lazy evaluation means it's never forced), and each
# test confirms dump.path was never created.

skip_if_not_installed("nimble")   # runNimble() attaches it before preparing inits

call_runNimble <- function(..., dump.path) {
  runNimble(parameters = "mu", max.samples.saved = 10, dump.path = dump.path, ...)
}

test_that("seed must be NULL or a single whole number", {
  for (bad in list("1", c(1, 2), NA, NA_integer_, 1.5, Inf, 1e10, list(1))) {
    d <- file.path(tempdir(), "never_created_seed")
    expect_error(call_runNimble(seed = bad, dump.path = d),
                 "seed must be NULL or a single whole number")
    expect_false(dir.exists(d))
  }
})

test_that("malformed inits stop runNimble() before anything is written", {
  d <- file.path(tempdir(), "never_created_inits")
  expect_error(call_runNimble(inits = list(list(mu = 0)), nc = 2, dump.path = d),
               "unnamed list of length 1.*nc = 2")
  expect_error(call_runNimble(inits = list(mu = 0, 1), dump.path = d),
               "some named and some unnamed")
  expect_error(call_runNimble(inits = list(chain1 = list(mu = 0), chain2 = list(mu = 1)),
                              dump.path = d),
               "list-valued element\\(s\\): chain1, chain2")
  expect_false(dir.exists(d))
})

test_that("an inits function is run in the parent, so its failures stop runNimble() before anything is written", {
  d <- file.path(tempdir(), "never_created_fn")
  expect_error(call_runNimble(inits = function() 1, dump.path = d),
               "return value for chain 1 must be a named list")
  expect_error(call_runNimble(inits = function() list(mu = undefined.in.session), dump.path = d),
               "inits function failed for chain 1:.*undefined.in.session")
  expect_false(dir.exists(d))
})

test_that("valid seed and inits get past validation, with the inits function run once per chain", {
  # Confirms the checks above discriminate rather than always erroring.
  # dump.path's parent directory doesn't exist, so the first write after
  # validation (the directive file) fails - deterministically, and well before
  # any worker could launch. In a temp working directory because runNimble()'s
  # on.exit() handler writes <mod.nam>_logs/ there.
  calls <- integer(0)
  f <- function(chn) { calls <<- c(calls, chn); list(mu = chn) }
  in_tempdir({
    expect_error(suppressWarnings(
                   call_runNimble(inits = f, nc = 3, seed = 1,
                                  dump.path = file.path("no_such_parent", "dump"))),
                 "cannot open")
  })
  expect_equal(calls, 1:3)
})
