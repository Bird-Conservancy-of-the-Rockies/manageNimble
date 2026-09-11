# nt2's positive-integer-multiple-of-nt requirement is validated at the very
# top of runNimble(), before dump.path/model.path/nimble/GNU parallel are ever
# touched - so, unlike nearly everything else in runNimble(), this can be
# tested with a genuine, executed call rather than only statically. Every
# other required argument is left unsupplied on purpose: R's lazy argument
# evaluation means they're never forced before this check's stop() fires, and
# the last test below confirms that directly (a valid nt2 proceeds past this
# check and fails on the next missing argument instead, so this isn't just an
# unconditional error).

test_that("nt2 not a positive integer multiple of nt is a clear, early error", {
  expect_error(
    runNimble(parameters = "a", parameters2 = "b", nt = 5, nt2 = 3, sav.model = TRUE),
    "nt2 \\(3\\) must be a positive integer multiple of nt \\(5\\)"
  )
  expect_error(
    runNimble(parameters = "a", parameters2 = "b", nt = 5, nt2 = 0, sav.model = TRUE),
    "positive integer multiple"
  )
  expect_error(
    runNimble(parameters = "a", parameters2 = "b", nt = 5, nt2 = 2.5, sav.model = TRUE),
    "positive integer multiple"
  )
})

test_that("nt2 validation is skipped entirely when parameters2 is empty", {
  # nt2 is forced to 1 internally when parameters2 is empty (see runNimble()'s
  # own top-of-function logic), so an explicit bad nt2 is simply never
  # consulted - proceeds past this check to the next missing argument instead.
  expect_error(
    runNimble(parameters = "a", parameters2 = character(), nt = 5, nt2 = 3, sav.model = TRUE),
    "max.samples.saved"
  )
})

test_that("a valid nt2 proceeds past this check (fails later, on an unrelated missing argument)", {
  # Confirms the check discriminates rather than always erroring regardless
  # of input.
  expect_error(
    runNimble(parameters = "a", parameters2 = "b", nt = 5, nt2 = 10, sav.model = TRUE),
    "max.samples.saved"
  )
})
