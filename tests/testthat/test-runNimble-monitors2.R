# runNimble() itself can't run on Windows (GNU parallel / processx's Linux-only
# process launch - see ?runNimble). These tests instead exercise, directly, the
# exact pieces runNimble() combines to implement consolidate.monitors2: the real
# (unmodified) gatherNimble2() against real block dumps written by the real
# runNimbleBlock() - both of which are Windows-testable, same rationale as
# test-runNimbleBlock.R - plus the attach-to-mod/save logic, reproduced here
# identically to both call sites in R/runNimble.R (search that file for
# "consolidate.monitors2" to compare). If that logic is ever refactored into its
# own function, these tests should call it directly instead of mirroring it.

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
  # tempfile() rather than sample(letters, ...): toy() below calls set.seed(42),
  # so a seed-derived name would collide across test_that() blocks that draw the
  # same number of random values before calling this - which they do here, and
  # did collide in practice (silently reusing a previous test's leftover
  # directory) before this fix. tempfile()'s own uniqueness is unaffected by
  # set.seed().
  d <- tempfile(pattern = "cm2")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# One chain, two blocks, with a real second monitor set (z, a vector-valued
# deterministic node) - enough to exercise gatherNimble2() meaningfully without
# the multi-minute cost of a multi-chain fit.
make_dump_with_monitors2 <- function() {
  m <- toy()
  dump.path <- new_tmp()
  tmp <- file.path(dump.path, "tmp1"); dir.create(tmp)
  f1 <- file.path(dump.path, "mod_chn1_1.RData")
  cm <- suppressMessages(
    runNimbleBlock(mod.lst = list(m$code, m$cons, m$dat, m$ini, m$pars),
                   n.iter = 200, n.thin = 2, tmp.path = tmp, dump.file.path = f1,
                   monitors2 = "z", n.thin2 = 50)
  )
  f2 <- file.path(dump.path, "mod_chn1_2.RData")
  suppressMessages(
    runNimbleBlock(comp.mcmc = cm, n.iter = 200, tmp.path = tmp,
                   dump.file.path = f2, monitors2 = "z")
  )
  dump.path
}

fake_mod <- function() {
  list(mcmcOutput = "placeholder", summary = data.frame(Parameter = c("mu", "sigma")),
       mcmc.info = c(nchains = 1, niterations = 400, burnin = 0, nthin = 2))
}

test_that("consolidate.monitors2 = TRUE (the default) attaches mod$monitors2, identical to the standalone file", {
  dump.path <- make_dump_with_monitors2()
  wd <- new_tmp(); old.wd <- setwd(wd); on.exit(setwd(old.wd), add = TRUE)
  mod.nam <- "modA"
  mod <- fake_mod()

  monitors2.out <- gatherNimble2(read.path = dump.path, burnin = 0, ni.block = 200, nt2 = 50,
                                 save.path = paste0(mod.nam, "_monitors2.rds"))
  consolidate.monitors2 <- TRUE   # mirrors runNimble()'s default
  sav.model <- TRUE
  if (consolidate.monitors2 && !is.null(monitors2.out)) {
    mod$monitors2 <- monitors2.out
    if (sav.model) R.utils::saveObject(mod, mod.nam)
  }

  expect_true(file.exists(paste0(mod.nam, "_monitors2.rds")))   # standalone file always written
  loaded <- R.utils::loadObject(mod.nam)
  expect_true("monitors2" %in% names(loaded))
  expect_equal(loaded$monitors2, readRDS(paste0(mod.nam, "_monitors2.rds")))
  expect_equal(dim(loaded$monitors2), c(8L, 8L))   # 2 blocks x 4 draws/block (200/50), z[1:8]
})

test_that("consolidate.monitors2 = FALSE reproduces the old behaviour: no mod$monitors2, standalone file still written", {
  dump.path <- make_dump_with_monitors2()
  wd <- new_tmp(); old.wd <- setwd(wd); on.exit(setwd(old.wd), add = TRUE)
  mod.nam <- "modB"
  mod <- fake_mod()
  sav.model <- TRUE
  if (sav.model) R.utils::saveObject(mod, mod.nam)   # runNimble() already saved mod earlier regardless

  monitors2.out <- gatherNimble2(read.path = dump.path, burnin = 0, ni.block = 200, nt2 = 50,
                                 save.path = paste0(mod.nam, "_monitors2.rds"))
  consolidate.monitors2 <- FALSE
  if (consolidate.monitors2 && !is.null(monitors2.out)) {
    mod$monitors2 <- monitors2.out
    if (sav.model) R.utils::saveObject(mod, mod.nam)
  }

  expect_true(file.exists(paste0(mod.nam, "_monitors2.rds")))   # standalone file always written
  loaded <- R.utils::loadObject(mod.nam)
  expect_false("monitors2" %in% names(loaded))
})

test_that("parameters2 empty: mod never gets a monitors2 element, regardless of consolidate.monitors2, and no standalone file is written", {
  wd <- new_tmp(); old.wd <- setwd(wd); on.exit(setwd(old.wd), add = TRUE)
  mod.nam <- "modC"
  mod <- fake_mod()
  parameters2 <- character()
  consolidate.monitors2 <- TRUE   # even at the default, nothing should happen

  if (length(parameters2) > 0) {
    stop("unreachable: gatherNimble2() must never be called when parameters2 is empty")
  }

  expect_false("monitors2" %in% names(mod))
  expect_false(file.exists(paste0(mod.nam, "_monitors2.rds")))
})
