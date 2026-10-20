# End-to-end check of `seed` and per-chain `inits` through the REAL generated
# worker script (ModRunScript.R), run in separate Rscript processes exactly as
# GNU parallel would run them - just launched here with processx instead, since
# GNU parallel isn't available on Windows. The parent-side preparation uses the
# same helpers runNimble() calls (deriveChainStreams(), prepareChainInits()),
# and the objects saved for the workers are read off runNimble()'s own
# save(list = ...) call, so the test can't drift from the real thing.
#
# Each run: 2 chains, automated mode with check.freq = 2, so every worker
# samples an initial block and one continuation block, then reports STOP. The
# test model file also saves the `inits` each worker actually used.
#
# Compiles 6 models in parallel; takes a few minutes. Set
# MANAGENIMBLE_TEST_NIMBLE=false to skip.

skip_if_not_installed("nimble")
skip_if(identical(tolower(Sys.getenv("MANAGENIMBLE_TEST_NIMBLE", "true")), "false"),
        "MANAGENIMBLE_TEST_NIMBLE=false")
skip_on_cran()

# The `list = c(...)` of runNimble()'s save() call - the objects workers load.
saved_object_names <- function() {
  found <- NULL
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("save")) && "list" %in% names(x)) found <<- eval(x[["list"]])
      for (i in seq_along(x)) if (!is.null(x[[i]])) try(walk(x[[i]]), silent = TRUE)
    }
  }
  walk(body(manageNimble::runNimble))
  stopifnot(!is.null(found))
  found
}

# When this package was loaded from source (devtools::test()), a subprocess's
# require(manageNimble) would find whatever older version is installed; load
# the source tree into the worker first so it runs the code under test. Under
# R CMD check the installed build is the code under test, so nothing to do.
dev_source_path <- function() {
  p <- getNamespaceInfo("manageNimble", "path")
  if (file.exists(file.path(p, "R", "runNimble.R"))) normalizePath(p, winslash = "/") else ""
}

# Builds dump.path exactly as runNimble() does before launching workers.
prepare_run <- function(seed, inits, nc = 2) {
  d <- normalizePath(new_dir("repro"), winslash = "/")
  model.file <- file.path(d, "model.R")
  writeLines(c(
    "model <<- nimble::nimbleCode({",
    "  mu ~ dnorm(0, sd = 10)",
    "  sigma ~ dunif(0, 10)",
    "  for (i in 1:N) y[i] ~ dnorm(mu, sd = sigma)",
    "})",
    "saveRDS(inits, file.path(dump.path, paste0('inits_seen_', chn, '.rds')))"
  ), model.file)

  e <- new.env()
  e$model.path <- model.file
  e$constants <- list(N = 8)
  e$data <- list(y = c(2.1, 3.4, 2.9, 3.8, 1.7, 3.1, 2.6, 3.3))
  e$parameters <- c("mu", "sigma")
  e$ni <- 200; e$nt <- 2
  e$dump.path <- d
  e$SamplerSourcePath <- NA
  e$check.freq <- 2
  e$automate.convergence.checks <- TRUE
  e$directive.file <- file.path(d, "runNimbleDirective.txt")
  e$parameters2 <- character(); e$nt2 <- 1
  e$rng.streams <- if (!is.null(seed)) deriveChainStreams(seed, nc) else NULL
  ci <- prepareChainInits(inits, nc, e$rng.streams)
  e$inits <- ci$inits
  e$inits.per.chain <- ci$per.chain
  save(list = saved_object_names(), envir = e, file = file.path(d, "NimbleObjects.RData"))

  writeLines(worker_script_lines(), file.path(d, "ModRunScript.R"))
  writeLines("GO", e$directive.file)
  for (cn in seq_len(nc)) writeLines("GO", file.path(d, paste0("block", cn, "Status.txt")))

  dev <- dev_source_path()
  writeLines(c(
    if (nzchar(dev)) paste0("pkgload::load_all(", deparse(dev), ", quiet = TRUE, export_all = FALSE)"),
    paste0("source(", deparse(file.path(d, "ModRunScript.R")), ")")
  ), file.path(d, "worker.R"))
  d
}

launch_workers <- function(d, nc = 2) {
  rscript <- file.path(R.home("bin"), "Rscript")
  lapply(seq_len(nc), function(cn)
    processx::process$new(rscript, c(file.path(d, "worker.R"), as.character(cn),
                                     file.path(d, "NimbleObjects.RData")),
                          stdout = file.path(d, paste0("worker", cn, ".out")),
                          stderr = "2>&1"))
}

# Waits for every chain in every run to report STOP (both blocks written),
# then releases the workers with a STOP directive.
finish_runs <- function(dirs, procs, nc = 2, timeout = 1200) {
  t0 <- Sys.time()
  repeat {
    errs <- unlist(lapply(dirs, list.files, pattern = "_error\\.log$", full.names = TRUE))
    if (length(errs) > 0) stop("Worker error: ", paste(readLines(errs[1]), collapse = " | "))
    dead <- vapply(procs, function(p) !p$is_alive() && !identical(p$get_exit_status(), 0L), logical(1))
    if (any(dead)) stop("A worker exited unexpectedly; see worker*.out under ", paste(dirs, collapse = ", "))
    status <- unlist(lapply(dirs, function(d) vapply(seq_len(nc), function(cn) {
      x <- tryCatch(readLines(file.path(d, paste0("block", cn, "Status.txt"))), error = function(e) "")
      if (length(x) == 1) x else ""
    }, character(1))))
    if (all(status == "STOP")) break
    if (as.numeric(Sys.time() - t0, units = "secs") > timeout) stop("Timed out waiting for workers.")
    Sys.sleep(2)
  }
  for (d in dirs) writeLines("STOP", file.path(d, "runNimbleDirective.txt"))
  for (p in procs) p$wait(60000)
}

block <- function(d, cn, b) {
  e <- new.env(); load(file.path(d, paste0("mod_chn", cn, "_", b, ".RData")), envir = e); e$samp
}
seen <- function(d, cn) readRDS(file.path(d, paste0("inits_seen_", cn, ".rds")))

test_that("seed + inits: same seed reproduces every draw and initial value; different seed/chain differs", {
  jitter <- function() list(mu = rnorm(1, 0, 2), sigma = runif(1, 0.5, 2))
  per.chain <- list(list(mu = -5, sigma = 1), list(mu = 5, sigma = 2))

  a1 <- prepare_run(123, jitter)
  a2 <- prepare_run(123, jitter)
  b  <- prepare_run(456, per.chain)
  procs <- c(launch_workers(a1), launch_workers(a2), launch_workers(b))
  on.exit(for (p in procs) if (p$is_alive()) p$kill(), add = TRUE)
  finish_runs(c(a1, a2, b), procs)

  for (cn in 1:2) for (blk in 1:2) {
    # same seed, separate processes: identical draws, initial and continuation blocks
    expect_identical(block(a1, cn, blk), block(a2, cn, blk))
  }
  expect_false(identical(block(a1, 1, 1), block(a1, 2, 1)))   # chains differ
  expect_false(identical(block(a1, 1, 1), block(a1, 1, 2)))   # stream continues, not reseeded
  expect_false(identical(block(a1, 1, 1), block(b, 1, 1)))    # seed matters
  expect_equal(nrow(block(a1, 1, 1)), 100)                    # ni / nt

  # generated inits: what each worker used is what the parent generated for
  # that chain, reproducibly, and the chains got different values
  gen <- prepareChainInits(jitter, 2, deriveChainStreams(123, 2))$inits
  for (cn in 1:2) {
    expect_identical(seen(a1, cn), gen[[cn]])
    expect_identical(seen(a2, cn), gen[[cn]])
  }
  expect_false(identical(seen(a1, 1), seen(a1, 2)))

  # per-chain list: each worker got its own element
  for (cn in 1:2) expect_identical(seen(b, cn), per.chain[[cn]])
})
