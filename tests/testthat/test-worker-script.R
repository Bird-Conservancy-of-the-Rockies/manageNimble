# The worker script is assembled as a character vector inside runNimble() and
# only ever executed inside a spawned Rscript process, where a syntax error or
# an undefined symbol is easy to miss. These tests check it statically.

test_that("the generated ModRunScript.R is syntactically valid R", {
  txt <- worker_script_lines()
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(txt, f)

  expect_no_error(parse(f))
})

test_that("the multi-line runNimbleBlock() calls survive the writeLines split", {
  # These calls are split across several elements of the character vector; a
  # misplaced comma fails only at runtime inside a worker.
  txt <- worker_script_lines()
  p <- parse(text = paste(txt, collapse = "\n"))

  calls <- character(0)
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("runNimbleBlock"))) {
        calls <<- c(calls, paste(deparse(x), collapse = " "))
      }
      for (i in seq_along(x)) if (!is.null(x[[i]])) try(walk(x[[i]]), silent = TRUE)
    }
  }
  for (e in p) walk(e)

  expect_equal(length(calls), 3)                       # 1 initial + 2 continuation
  expect_true(any(grepl("n\\.thin2", calls)))          # initial block configures thin2
  expect_true(all(grepl("monitors2", calls)))          # all three pass monitors2
})

test_that("every top-level symbol the worker reads is either saved to disk, assigned, or a local helper", {
  saved <- c("model.path", "constants", "data", "inits", "parameters", "ni", "nt",
             "dump.path", "SamplerSourcePath", "check.freq",
             "automate.convergence.checks", "directive.file", "parameters2", "nt2")
  assigned <- c("chn", "path.NimbleWorkspace", "i", "dump.file.path", "mod.comp",
                "status.file", "status.chain", "i.stop", "directive", "error.log")
  # read_state()/write_state() are defined in-script (finding F26); e/path/text/
  # tmp/x are their own parameters and locals, picked up by the crude
  # all.vars() scan below because it doesn't understand function scoping.
  helpers <- c("read_state", "write_state", "e", "path", "text", "tmp", "x")
  from_model_file <- "model"
  pkgs <- c("nimble", "manageNimble")

  p <- parse(text = paste(worker_script_lines(), collapse = "\n"))
  syms <- unique(unlist(lapply(p, all.vars)))
  undefined <- setdiff(syms, c(saved, assigned, helpers, from_model_file, pkgs))

  # FIXED (F14): this used to always include "cn" - the "undefined condition
  # reached" branch referenced the PARENT's loop variable, which does not
  # exist in the worker, so that branch died with "object 'cn' not found"
  # instead of ever reporting the error. It now raises stop() with the
  # worker's own `chn`, caught by the enclosing tryCatch like any other error.
  expect_equal(undefined, character(0))
})

test_that("the tryCatch/error-log convention actually catches and logs a worker error", {
  # Not just a text-pattern check: splices a forced stop() into the real
  # generated wrapper (keeping the actual error handler unchanged) and runs it
  # in a real Rscript subprocess, since the handler ends with quit() which
  # would otherwise kill the test process.
  txt <- worker_script_lines()
  i.try   <- which(txt == "tryCatch({")
  i.catch <- which(txt == "}, error = function(e) {")
  i.error.log <- which(grepl("^error\\.log <-", txt))
  stopifnot(length(i.try) == 1, length(i.catch) == 1, length(i.error.log) == 1,
            i.catch > i.try, i.error.log < i.try)

  d <- new_dir("errlogtest")
  script <- c(
    "chn <- '9'",
    paste0("dump.path <- ", deparse(d)),
    txt[i.error.log],
    "tryCatch({",
    "  stop('deliberate test error for F25/F27 log verification')",
    txt[i.catch:length(txt)]
  )
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(script, f)

  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(rscript, args = shQuote(f), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")

  log.file <- file.path(d, "chn9_error.log")
  expect_true(file.exists(log.file))
  contents <- readLines(log.file)
  expect_true(any(grepl("^Time: ", contents)))
  expect_true(any(grepl("^Chain: 9$", contents)))
  expect_true(any(grepl("deliberate test error", contents)))
  expect_equal(status, 1L)
})

test_that("the objects runNimble saves are exactly the ones the worker loads", {
  # The save() call and the worker script have to be kept in step by hand.
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl('"parameters2", "nt2"', b, fixed = TRUE))

  txt <- paste(worker_script_lines(), collapse = "\n")
  expect_true(grepl("monitors2 = parameters2", txt, fixed = TRUE))
  expect_true(grepl("n.thin2 = nt2", txt, fixed = TRUE))
})

test_that("the whole worker body is wrapped in tryCatch, logging to chn<c>_error.log", {
  txt <- paste(worker_script_lines(), collapse = "\n")

  expect_true(grepl("error.log <- paste0(dump.path, '/chn', chn, '_error.log')", txt, fixed = TRUE))
  expect_true(grepl("tryCatch({", txt, fixed = TRUE))
  expect_true(grepl("}, error = function(e) {", txt, fixed = TRUE))
  expect_true(grepl("quit(save = 'no', status = 1)", txt, fixed = TRUE))
})

test_that("defaults leave the second monitor set switched off", {
  fm <- formals(manageNimble::runNimble)
  expect_equal(eval(fm$parameters2), character(0))
  expect_null(fm$nt2)

  fb <- formals(manageNimble::runNimbleBlock)
  expect_equal(eval(fb$monitors2), character(0))
  expect_equal(fb$n.thin2, 1)
})

# ---------------------------------------------------------------------------
# Regression tests for fixed bugs (see NEWS.md for the full write-up)
# ---------------------------------------------------------------------------

test_that("FIXED (F1): SamplerSourcePath is passed as a real argument to runNimbleBlock", {
  # Previously runNimble built
  #   runNimbleBlock(mod.lst = list(model, constants, data, inits, parameters,
  #                                 SamplerSourcePath = SamplerSourcePath), ...)
  # so it became the sixth element of mod.lst. runNimbleBlock's own
  # SamplerSourcePath argument kept its default NA, and neither
  # require(nimbleHMC) nor source(SamplerSourcePath) ever ran.
  txt <- paste(worker_script_lines(), collapse = "\n")

  expect_false(grepl("list(model, constants, data, inits, parameters, SamplerSourcePath = SamplerSourcePath)",
                     txt, fixed = TRUE))
  expect_true(grepl("list(model, constants, data, inits, parameters)", txt, fixed = TRUE))

  p <- parse(text = txt)
  found <- NULL
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("runNimbleBlock")) && "mod.lst" %in% names(x)) found <<- x
      for (i in seq_along(x)) if (!is.null(x[[i]])) try(walk(x[[i]]), silent = TRUE)
    }
  }
  for (e in p) walk(e)

  expect_true("SamplerSourcePath" %in% names(found))            # now a real argument of the call
  expect_false("SamplerSourcePath" %in% names(found$mod.lst))    # no longer buried in mod.lst
})

test_that("FIXED (F26): worker coordination reads retry, writes are atomic", {
  # readLines() on a file being rewritten by another process can return
  # character(0) or more than one line. Previously every conditional in the
  # worker treated that as a real state and crashed; the parent already
  # wrapped its own reads in try()+retry, but the worker had no equivalent.
  # read_state()'s own definition legitimately contains a bare readLines()
  # call (that's its implementation) - what must be gone is every OTHER call
  # site that used to read the directive/status files directly.
  txt <- paste(worker_script_lines(), collapse = "\n")

  expect_true(grepl("read_state <- function", txt, fixed = TRUE))
  expect_true(grepl("write_state <- function", txt, fixed = TRUE))
  expect_false(grepl("readLines(status.file)", txt, fixed = TRUE))
  expect_false(grepl("readLines(directive.file)", txt, fixed = TRUE))
  expect_false(grepl("writeLines('STOP', status.file)", txt, fixed = TRUE))
  expect_false(grepl("writeLines('GO', status.file)", txt, fixed = TRUE))
  expect_true(grepl("read_state(status.file)", txt, fixed = TRUE))
  expect_true(grepl("read_state(directive.file)", txt, fixed = TRUE))
  expect_true(grepl("write_state('STOP', status.file)", txt, fixed = TRUE))
  expect_true(grepl("write_state('GO', status.file)", txt, fixed = TRUE))

  # what a torn read used to do to a bare conditional, for context
  empty <- character(0)
  expect_error(if (empty == "GO") TRUE, "argument is of length zero")
  expect_error(while (empty != "STOP") break, "argument is of length zero")
  expect_error(if (c("GO", "GO") == "GO") TRUE, "condition has length > 1")

  # the parent's directive/status writes are atomic too now
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("write_state <- function(text, path)", b, fixed = TRUE))
  expect_false(grepl('writeLines("GO", directive.file)', b, fixed = TRUE))
})

test_that("FIXED (F10): rtrn.model assigns the model object under mod.nam, not the name itself", {
  b <- fn_text(manageNimble::runNimble)

  expect_true(grepl("assign(mod.nam, mod, envir = .GlobalEnv)", b, fixed = TRUE))
  expect_false(grepl('assign("mod", mod.nam, envir = .GlobalEnv)', b, fixed = TRUE))
})

test_that("FIXED (F21): dplyr/stringr/coda/mcmcOutput/R.utils are declared Imports, not just Suggests", {
  # Previously these were listed under Suggests: with no require() calls for
  # most of them, so runNimble() failed inside the search-path fallback the
  # moment a caller hadn't happened to attach dplyr/stringr themselves - and
  # by the time it failed, process$new() had already launched nc workers that
  # would then sample forever against nobody watching.
  #
  # The fix here is NAMESPACE imports (DESCRIPTION Imports: + importFrom() in
  # NAMESPACE), not require() calls in the function body: package code
  # resolves bare names through its own declared imports regardless of what
  # the caller has attached, so runNimble() needs no require(dplyr) of its
  # own for this to be reliable - and R CMD check flags require() calls for
  # anything already available via NAMESPACE as bad practice.
  desc <- read.dcf(system.file("DESCRIPTION", package = "manageNimble"))
  imports <- trimws(strsplit(desc[1, "Imports"], ",")[[1]])
  imports <- sub("\\s*\\(.*\\)$", "", imports)
  expect_true(all(c("dplyr", "stringr", "coda", "mcmcOutput", "R.utils") %in% imports))

  imp <- getNamespaceImports("manageNimble")
  expect_true("filter" %in% imp$dplyr)
  expect_true("str_detect" %in% imp$stringr)
  expect_true("as.mcmc" %in% imp$coda)
  expect_true("mcmcOutput" %in% imp$mcmcOutput)

  # bare require() calls for these are no longer in runNimble()'s own body -
  # package code resolves them through the imports above regardless of the
  # caller's session, so a require() here would only be redundant and would
  # trip R CMD check's "declared dependencies" warning.
  b <- fn_text(manageNimble::runNimble)
  expect_false(grepl("require(dplyr)", b, fixed = TRUE))
  expect_false(grepl("require(stringr)", b, fixed = TRUE))
  expect_false(grepl("require(coda)", b, fixed = TRUE))
})

test_that("FIXED (F29): the initial wait counts distinct chains, not raw file count", {
  # Previously `sum(str_detect(list.files(dump.path), "mod_chn")) < nc` was
  # satisfied by nc files from a single fast chain.
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("count_distinct_chains()", b, fixed = TRUE))
  expect_false(grepl('sum(str_detect(list.files(dump.path), "mod_chn")) < nc', b, fixed = TRUE))

  files <- c("mod_chn1_1.RData", "mod_chn1_2.RData")     # chain 2 has written nothing
  expect_gte(sum(stringr::str_detect(files, "mod_chn")), 2)
  expect_lt(length(unique(stringr::str_split(files, "_", simplify = TRUE)[, 2])), 2)
})

test_that("FIXED (F12): expecting the wrong number of chains is caught in both directions", {
  # Previously only `> nc` was checked; a chain that silently wrote nothing
  # (`< nc`) passed uncaught.
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("length(unique(check.blocks$m[, 1])) != nc", b, fixed = TRUE))
  expect_false(grepl("length(unique(check.blocks$m[, 1])) > nc", b, fixed = TRUE))
})

test_that("FIXED (F25/F27/F30): worker failures are captured and workers can't be orphaned", {
  b <- fn_text(manageNimble::runNimble)

  expect_true(grepl("check_workers_alive <- function", b, fixed = TRUE))
  expect_true(grepl("--results", b, fixed = TRUE))
  expect_true(grepl("--joblog", b, fixed = TRUE))
  expect_true(grepl("on.exit({", b, fixed = TRUE))
  expect_true(grepl("kill_tree()", b, fixed = TRUE))
  expect_true(grepl("copy_logs <- function", b, fixed = TRUE))
})

test_that("FIXED (F28): proc is a local variable, not assigned into .GlobalEnv", {
  b <- fn_text(manageNimble::runNimble)
  expect_false(grepl("proc <<-", b, fixed = TRUE))
  expect_true(grepl("proc <- NULL", b, fixed = TRUE))
  expect_true(grepl("proc <- process$new", b, fixed = TRUE))
})
