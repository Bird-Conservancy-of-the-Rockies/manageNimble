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
             "automate.convergence.checks", "directive.file", "parameters2", "nt2",
             "inits.per.chain", "rng.streams")
  assigned <- c("chn", "path.NimbleWorkspace", "i", "dump.file.path", "mod.comp",
                "status.file", "status.chain", "i.stop", "directive", "error.log")
  # read_state()/write_state() are defined in-script; e/path/text/
  # tmp/x are their own parameters and locals, picked up by the crude
  # all.vars() scan below because it doesn't understand function scoping.
  helpers <- c("read_state", "write_state", "e", "path", "text", "tmp", "x")
  from_model_file <- "model"
  pkgs <- c("nimble", "manageNimble")

  p <- parse(text = paste(worker_script_lines(), collapse = "\n"))
  syms <- unique(unlist(lapply(p, all.vars)))
  undefined <- setdiff(syms, c(saved, assigned, helpers, from_model_file, pkgs))

  # An undefined name would only fail when its line runs inside a worker. The
  # "undefined condition reached" branch, for example, raises stop() with the
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
    "  stop('deliberate test error for error-log verification')",
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
  expect_true(grepl('"parameters2", "nt2", "inits.per.chain", "rng.streams"', b, fixed = TRUE))

  txt <- paste(worker_script_lines(), collapse = "\n")
  expect_true(grepl("monitors2 = parameters2", txt, fixed = TRUE))
  expect_true(grepl("n.thin2 = nt2", txt, fixed = TRUE))
  expect_true(grepl("rng.streams[[as.integer(chn)]]", txt, fixed = TRUE))
  expect_true(grepl("if(inits.per.chain) inits <- inits[[as.integer(chn)]]", txt, fixed = TRUE))
})

test_that("the worker installs its RNG stream and selects its inits before anything can draw a random number", {
  # Order matters: the stream must be in place before the model file is
  # sourced, the model built/compiled, or the first block sampled; per-chain
  # inits must be selected before the first runNimbleBlock() call reads them.
  # Both inside the tryCatch, so a failure is logged like any other.
  txt <- worker_script_lines()
  i.try    <- which(txt == "tryCatch({")
  i.seed   <- which(txt == "if(!is.null(rng.streams)) {")
  i.sel    <- which(txt == "if(inits.per.chain) inits <- inits[[as.integer(chn)]]")
  i.source <- which(txt == "source(model.path)")
  i.block  <- which(grepl("runNimbleBlock(", txt, fixed = TRUE))[1]
  expect_length(i.seed, 1)
  expect_length(i.sel, 1)
  expect_true(i.try < i.seed && i.seed < i.sel && i.sel < i.source && i.source < i.block)

  # The worker's RNGkind() call must name the same generator, normal and
  # sample kinds the parent derived the streams under.
  i.kind <- which(grepl("^RNGkind\\(", txt))
  expect_length(i.kind, 1)
  call <- parse(text = txt[i.kind])[[1]]
  expect_equal(unname(vapply(as.list(call)[-1], as.character, character(1))),
               unname(RNG_KIND))
})

test_that("defaults leave seeding off", {
  expect_null(formals(manageNimble::runNimble)$seed)
})

test_that("seeded manual mode caps the gathered block count", {
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("nblks.cap <- if (!is.null(seed) && !automate.convergence.checks)", b, fixed = TRUE))
  expect_true(grepl("if (!is.null(nblks.cap) && nblks >= nblks.cap) nblks <- nblks.cap", b, fixed = TRUE))
  expect_true(grepl("max.nblks = nblks.cap", b, fixed = TRUE))
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
# Worker coordination and runNimble() structure
# ---------------------------------------------------------------------------

test_that("SamplerSourcePath is passed as a real argument to runNimbleBlock, not inside mod.lst", {
  # runNimbleBlock() only sources SamplerSourcePath (and loads nimbleHMC) when
  # it arrives as its own argument; mod.lst has no slot for it.
  txt <- paste(worker_script_lines(), collapse = "\n")
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

  expect_true("SamplerSourcePath" %in% names(found))
  expect_false("SamplerSourcePath" %in% names(found$mod.lst))
})

test_that("worker coordination reads retry, writes are atomic", {
  # readLines() on a file being rewritten by another process can return
  # character(0) or more than one line, and a bare conditional on either
  # errors (demonstrated below). Every worker read of the directive/status
  # files goes through read_state(), which retries until it gets exactly one
  # line; every write goes through write_state() (temp file + rename).
  # read_state()'s own definition contains the one bare readLines() call.
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

  # a bare conditional on a torn read
  empty <- character(0)
  expect_error(if (empty == "GO") TRUE, "argument is of length zero")
  expect_error(while (empty != "STOP") break, "argument is of length zero")
  expect_error(if (c("GO", "GO") == "GO") TRUE, "condition has length > 1")

  # the parent's directive/status writes are atomic too
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("write_state <- function(text, path)", b, fixed = TRUE))
  expect_false(grepl('writeLines("GO", directive.file)', b, fixed = TRUE))
})

test_that("rtrn.model assigns the model object under mod.nam", {
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("assign(mod.nam, mod, envir = .GlobalEnv)", b, fixed = TRUE))
})

test_that("dplyr/stringr/coda/mcmcOutput/R.utils are declared Imports", {
  # Package code resolves bare names through its own declared imports
  # (DESCRIPTION Imports: + importFrom() in NAMESPACE) regardless of what the
  # caller has attached, so runNimble() works in a session that has not
  # attached dplyr/stringr. require() calls for these in runNimble()'s body
  # would be redundant and trip R CMD check's "declared dependencies" warning.
  desc <- read.dcf(system.file("DESCRIPTION", package = "manageNimble"))
  imports <- trimws(strsplit(desc[1, "Imports"], ",")[[1]])
  imports <- sub("\\s*\\(.*\\)$", "", imports)
  expect_true(all(c("dplyr", "stringr", "coda", "mcmcOutput", "R.utils") %in% imports))

  imp <- getNamespaceImports("manageNimble")
  expect_true("filter" %in% imp$dplyr)
  expect_true("str_detect" %in% imp$stringr)
  expect_true("as.mcmc" %in% imp$coda)
  expect_true("mcmcOutput" %in% imp$mcmcOutput)

  b <- fn_text(manageNimble::runNimble)
  expect_false(grepl("require(dplyr)", b, fixed = TRUE))
  expect_false(grepl("require(stringr)", b, fixed = TRUE))
  expect_false(grepl("require(coda)", b, fixed = TRUE))
})

test_that("the initial wait counts distinct chains, not raw file count", {
  # nc files can come from a single fast chain, so runNimble() waits until nc
  # distinct chain numbers appear among the block dumps.
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("count_distinct_chains()", b, fixed = TRUE))

  files <- c("mod_chn1_1.RData", "mod_chn1_2.RData")     # chain 2 has written nothing
  expect_gte(sum(stringr::str_detect(files, "mod_chn")), 2)
  expect_lt(length(unique(stringr::str_split(files, "_", simplify = TRUE)[, 2])), 2)
})

test_that("count_distinct_chains' own pattern excludes an in-progress .tmp dump", {
  # Extracts the actual pattern string from the parsed call (evaluating the
  # literal AST node rather than deparse()'d text, which would double-escape
  # its backslash). A readdir caught mid-rename can list a block's .tmp
  # alongside its already-renamed final file - an unanchored "mod_chn" match
  # would count both as the same chain, which happens to still pass the
  # distinct-chain-count check, but would feed a duplicate block number into
  # countNimbleBlocks() downstream.
  found <- NULL
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("str_detect")) && identical(x[[2]], as.name("fl")))
        found <<- eval(x[[3]])
      for (i in seq_along(x)) if (!is.null(x[[i]])) try(walk(x[[i]]), silent = TRUE)
    }
  }
  walk(body(manageNimble::runNimble))
  expect_false(is.null(found))

  expect_true(stringr::str_detect("mod_chn1_1.RData", found))
  expect_false(stringr::str_detect("mod_chn1_1.RData.tmp", found))
})

test_that("a chain count different from nc is caught in both directions", {
  # Catches both an unexpected extra chain and a chain that wrote nothing.
  b <- fn_text(manageNimble::runNimble)
  expect_true(grepl("length(unique(check.blocks$m[, 1])) != nc", b, fixed = TRUE))
})

test_that("worker failures are captured and workers can't be orphaned", {
  b <- fn_text(manageNimble::runNimble)

  expect_true(grepl("check_workers_alive <- function", b, fixed = TRUE))
  expect_true(grepl("--results", b, fixed = TRUE))
  expect_true(grepl("--joblog", b, fixed = TRUE))
  expect_true(grepl("on.exit({", b, fixed = TRUE))
  expect_true(grepl("kill_tree()", b, fixed = TRUE))
  expect_true(grepl("copy_logs <- function", b, fixed = TRUE))
})

test_that("proc is a local variable, not assigned into .GlobalEnv", {
  b <- fn_text(manageNimble::runNimble)
  expect_false(grepl("proc <<-", b, fixed = TRUE))
  expect_true(grepl("proc <- NULL", b, fixed = TRUE))
  expect_true(grepl("proc <- process$new", b, fixed = TRUE))
})
