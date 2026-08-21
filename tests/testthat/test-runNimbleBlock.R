# These tests compile real NIMBLE models and take a few minutes. They need a
# working C++ toolchain but nothing else - no GNU parallel, no lsof - so they
# run on Windows. Set MANAGENIMBLE_TEST_NIMBLE=false to skip.

skip_if_not_installed("nimble")
skip_if(identical(tolower(Sys.getenv("MANAGENIMBLE_TEST_NIMBLE", "true")), "false"),
        "MANAGENIMBLE_TEST_NIMBLE=false")
skip_on_cran()

toy <- function() {
  code <- nimble::nimbleCode({
    mu ~ dnorm(0, sd = 10)
    sigma ~ dunif(0, 10)
    for (i in 1:N) {
      y[i] ~ dnorm(mu, sd = sigma)
      z[i] <- y[i] - mu
    }
  })
  set.seed(42)
  list(code = code, cons = list(N = 8), dat = list(y = rnorm(8, 3, 1)),
       ini = list(mu = 0, sigma = 1), pars = c("mu", "sigma"))
}

new_tmp <- function() {
  d <- file.path(tempdir(), paste0("nb", paste(sample(letters, 8, TRUE), collapse = "")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

test_that("without monitors2 the dump contains only samp, at the primary thin rate", {
  m <- toy()
  f <- file.path(new_tmp(), "b1.RData")

  suppressMessages(
    runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars),
                   n.iter = 200, n.thin = 2, tmp.path = new_tmp(),
                   dump.file.path = f)
  )

  e <- new.env(); load(f, envir = e)
  expect_setequal(ls(e), "samp")
  expect_equal(dim(e$samp), c(100L, 2L))
  expect_setequal(colnames(e$samp), c("mu", "sigma"))
})

test_that("monitors2 adds samp2 at its own thin rate, and continuation blocks do not duplicate", {
  m <- toy()
  tp <- new_tmp()
  f1 <- file.path(new_tmp(), "b1.RData")
  f2 <- file.path(new_tmp(), "b2.RData")

  cm <- suppressMessages(
    runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars),
                   n.iter = 200, n.thin = 2, tmp.path = tp, dump.file.path = f1,
                   monitors2 = "z", n.thin2 = 50)
  )

  e1 <- new.env(); load(f1, envir = e1)
  expect_setequal(ls(e1), c("samp", "samp2"))
  expect_equal(dim(e1$samp), c(100L, 2L))       # 200 / 2
  expect_equal(dim(e1$samp2), c(4L, 8L))        # 200 / 50, z[1:8]
  expect_setequal(colnames(e1$samp2), paste0("z[", 1:8, "]"))

  # Continuation block, called exactly as the generated worker script does:
  # comp.mcmc supplied, n.thin and n.thin2 omitted, monitors2 still passed.
  suppressMessages(
    runNimbleBlock(comp.mcmc = cm, n.iter = 200, tmp.path = tp,
                   dump.file.path = f2, monitors2 = "z")
  )

  e2 <- new.env(); load(f2, envir = e2)
  # resetMV = TRUE must clear mvSamples2 as well as mvSamples, otherwise every
  # block would re-emit all previous draws and gatherNimble2 would duplicate.
  expect_equal(dim(e2$samp), c(100L, 2L))
  expect_equal(dim(e2$samp2), c(4L, 8L))
  expect_false(isTRUE(all.equal(unname(e1$samp2), unname(e2$samp2))))
  expect_false(isTRUE(all.equal(unname(e1$samp), unname(e2$samp))))
})

test_that("exactly one of mod.lst and comp.mcmc must be supplied", {
  m <- toy()
  expect_error(runNimbleBlock(n.iter = 10, tmp.path = new_tmp(),
                              dump.file.path = tempfile()))
  expect_error(runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars),
                              comp.mcmc = 1, n.iter = 10, tmp.path = new_tmp(),
                              dump.file.path = tempfile()))
})

test_that("runNimbleBlock() only sources SamplerSourcePath as a top-level argument, never from inside mod.lst", {
  # runNimbleBlock() itself was never buggy here - finding F1 was entirely
  # about runNimble()'s worker-script generation calling it with
  # SamplerSourcePath buried inside mod.lst (see test-worker-script.R's
  # "FIXED (F1)" test, which confirms the generated script no longer does
  # that). This test documents the contract runNimbleBlock() has always had
  # and that fix now correctly relies on: mod.lst has no slot for it.
  m <- toy()
  sf <- file.path(new_tmp(), "sampler.R")
  writeLines("assign('MN_SAMPLER_SOURCED', TRUE, envir = globalenv())", sf)

  # buried in mod.lst - ignored by design (mod.lst[[6]] is simply never read)
  assign("MN_SAMPLER_SOURCED", FALSE, envir = globalenv())
  on.exit(suppressWarnings(rm("MN_SAMPLER_SOURCED", envir = globalenv())), add = TRUE)
  suppressMessages(
    runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars,
                                  SamplerSourcePath = sf),
                   n.iter = 50, n.thin = 1, tmp.path = new_tmp(),
                   dump.file.path = file.path(new_tmp(), "b.RData"))
  )
  expect_false(get("MN_SAMPLER_SOURCED", envir = globalenv()))   # <- never sourced

  # the shape that would have worked
  assign("MN_SAMPLER_SOURCED", FALSE, envir = globalenv())
  suppressMessages(
    runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars),
                   SamplerSourcePath = sf,
                   n.iter = 50, n.thin = 1, tmp.path = new_tmp(),
                   dump.file.path = file.path(new_tmp(), "b.RData"))
  )
  expect_true(get("MN_SAMPLER_SOURCED", envir = globalenv()))
})

test_that("FIXED: a SamplerSourcePath file's nm.conf$removeSamplers()/addSampler() calls actually reach the nm.conf used to build the MCMC", {
  # source(SamplerSourcePath) used to default to local = FALSE, evaluating the
  # sourced file in .GlobalEnv - but nm.conf is local to this function's own
  # call frame. In the normal case (a fresh worker process with no pre-existing
  # global nm.conf) that meant an immediate "object 'nm.conf' not found" error
  # on every real (non-trivial) SamplerSourcePath file; the test above never
  # caught this because its fixture never references nm.conf at all. This test
  # uses a genuine sampler-swap file instead, and confirms both that it no
  # longer errors and that the swap lands on the real nm.conf: the sourced file
  # records nm.conf$getSamplers() (nimble's own sampler-config introspection)
  # into a global side-channel immediately after the swap, which only shows
  # "slice" for mu if the local nm.conf here was the one actually mutated - a
  # stale/unrelated global object would leave it at the model's default
  # ("conjugate_dnorm_dnorm_identity" for this toy model).
  m <- toy()
  sf <- file.path(new_tmp(), "sampler.R")
  writeLines(c(
    "nm.conf$removeSamplers('mu')",
    "nm.conf$addSampler(target = 'mu', type = 'slice')",
    "assign('MN_SAMPLER_NAMES', sapply(nm.conf$getSamplers(), function(s) s$name), envir = globalenv())"
  ), sf)

  assign("MN_SAMPLER_NAMES", NULL, envir = globalenv())
  on.exit(suppressWarnings(rm("MN_SAMPLER_NAMES", envir = globalenv())), add = TRUE)

  expect_no_error(suppressMessages(
    runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars),
                   SamplerSourcePath = sf,
                   n.iter = 50, n.thin = 1, tmp.path = new_tmp(),
                   dump.file.path = file.path(new_tmp(), "b.RData"))
  ))

  samplers <- get("MN_SAMPLER_NAMES", envir = globalenv())
  names(samplers) <- NULL
  expect_setequal(samplers, c("RW", "slice"))   # sigma unchanged, mu swapped from conjugate to slice
})
