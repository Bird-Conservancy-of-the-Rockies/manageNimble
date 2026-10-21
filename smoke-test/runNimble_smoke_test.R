# End-to-end smoke test of manageNimble 0.4.0 on a Linux server: real
# runNimble() calls (GNU parallel, separate worker processes, the parent's
# polling loop) on a tiny toy model, checking per-chain inits, seeding, and
# that the second monitor set and its iter.key join work. Takes roughly
# 10-20 minutes (9 runs, each compiling the model once per chain).
#
# Run from the shell, not from inside an interactive R session, so the
# library you set is the one the worker processes see too:
#
#   R_LIBS=$HOME/R/mn-0.4.0 Rscript smoke-test/runNimble_smoke_test.R [output.dir]
#
# output.dir defaults to manageNimble_smoke_test/ next to the repo (outside
# it); each run writes to its own timestamped subdirectory there. Exits with
# status 0 if every check passes, 1 otherwise; results are printed and saved
# as smoke_test_results.csv.

expected.version <- "0.4.0"
nc <- 2
run.timeout <- 1800   # seconds; a run that hangs is killed (with its workers) and reported as FAIL

# --- Locations ----------------------------------------------------------------
script.path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script.path) != 1) stop("Run this with Rscript from the shell - see the header.")
repo.dir <- dirname(dirname(normalizePath(script.path)))
args <- commandArgs(trailingOnly = TRUE)
out.dir <- if (length(args) >= 1) args[[1]] else file.path(dirname(repo.dir), "manageNimble_smoke_test")
out.dir <- normalizePath(out.dir, mustWork = FALSE)
if (out.dir == repo.dir || startsWith(out.dir, paste0(repo.dir, "/")))
  stop("output.dir must be outside the repo (", repo.dir, ").")
run.dir <- file.path(out.dir, format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)
cat("Results go to:", run.dir, "\n\n")

# --- Recording ----------------------------------------------------------------
results <- data.frame(check = character(), status = character(), detail = character())
record <- function(check, status, detail = "") {
  if (is.logical(status)) status <- if (isTRUE(status)) "PASS" else "FAIL"
  results[nrow(results) + 1, ] <<- list(check, status, as.character(detail))
  cat(sprintf("[%s] %s%s\n", status, check, if (nzchar(detail)) paste0(" - ", detail) else ""))
}
finish <- function() {
  write.csv(results, file.path(run.dir, "smoke_test_results.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(run.dir, "sessionInfo.txt"))
  failed <- sum(results$status == "FAIL")
  cat("\n", if (failed == 0) "SMOKE TEST PASSED" else paste("SMOKE TEST FAILED:", failed, "check(s)"),
      "\nResults:", file.path(run.dir, "smoke_test_results.csv"), "\n")
  quit(save = "no", status = if (failed == 0) 0 else 1)
}

# --- Preflight: is this the right package, for the parent AND the workers? ----
cat("== Preflight ==\n")
record("Linux/unix platform", .Platform$OS.type == "unix", R.version$platform)
gnu.parallel <- nzchar(Sys.which("parallel")) &&
  any(grepl("GNU parallel", suppressWarnings(system2("parallel", "--version", stdout = TRUE, stderr = TRUE))))
record("GNU parallel on PATH", gnu.parallel, Sys.which("parallel"))
record("nimble installed", requireNamespace("nimble", quietly = TRUE))

parent.version <- tryCatch(as.character(packageVersion("manageNimble")), error = function(e) "not installed")
parent.path <- tryCatch(find.package("manageNimble"), error = function(e) "")
record("this session loads manageNimble 0.4.0", parent.version == expected.version,
       paste(parent.version, "from", parent.path))

# Workers are launched by GNU parallel as plain `Rscript` from the PATH, in new
# processes that only see the library through the environment (R_LIBS), not
# this session's .libPaths(). Ask such a process directly.
rscript <- Sys.which("Rscript")
probe <- suppressWarnings(system2(rscript, c("-e", shQuote(paste0(
  "cat(tryCatch(as.character(packageVersion(\"manageNimble\")), error = function(e) \"not installed\"), ",
  "tryCatch(find.package(\"manageNimble\"), error = function(e) \"\"), R.version.string, sep = \"\\n\")"))),
  stdout = TRUE, stderr = FALSE))
worker.version <- if (length(probe) >= 1) probe[1] else "no output"
record("a fresh Rscript (as a worker) loads manageNimble 0.4.0", worker.version == expected.version,
       paste(worker.version, "from", if (length(probe) >= 2) probe[2] else "?", "using", rscript))
record("workers load the same installed copy as this session",
       length(probe) >= 2 && identical(normalizePath(probe[2], mustWork = FALSE), normalizePath(parent.path, mustWork = FALSE)))
record("workers run the same R version as this session",
       length(probe) >= 3 && identical(probe[3], R.version.string), if (length(probe) >= 3) probe[3] else "")

if (any(results$status == "FAIL")) {
  cat("\nPreflight failed - not running anything. See the header for how to set R_LIBS.\n")
  finish()
}

# --- Toy model ----------------------------------------------------------------
setwd(run.dir)   # runNimble() saves mod.nam, *_monitors2.rds and *_logs/ in the working directory
model.file <- file.path(run.dir, "toy_model.R")
writeLines(c(
  "model <<- nimble::nimbleCode({",
  "  mu ~ dnorm(0, sd = 10)",
  "  sigma ~ dunif(0, 10)",
  "  for (i in 1:N) {",
  "    y[i] ~ dnorm(mu, sd = sigma)",
  "    z[i] <- y[i] - mu",
  "  }",
  "})"
), model.file)
y <- c(2.1, 3.4, 2.9, 3.8, 1.7, 3.1, 2.6, 3.3, 2.2, 4.0,
       2.8, 3.6, 1.9, 3.0, 2.5, 3.9, 2.4, 3.2, 2.7, 3.5)

# Automated mode; max.samples.saved forces additional thinning, so the second
# monitor set's correspondence with the primary set is exercised under a cap.
# ni, nt and nb give (ni - ni * nb) / nt = 240 pre-cap draws per chain (350
# under D's absolute nb = 100 below), comfortably above max.samples.saved = 150
# so the cap actually engages; max.samples.saved itself clears mcmcOutput's own
# threshold for calculating Rhat (more than 100 draws per chain) once applied.
automated <- list(model.path = model.file, data = list(y = y), constants = list(N = 20),
                  parameters = c("mu", "sigma"), parameters2 = "z", nt2 = 20,
                  nc = nc, ni = 800, nt = 2, nb = 0.4, check.freq = 1, max.tries = 3,
                  max.samples.saved = 150, Rht.required = 1.1, neff.required = 50,
                  inits = list(mu = 0, sigma = 1), sav.model = TRUE, rtrn.model = FALSE)
# Sets each argument whole, so a list-valued one (such as a per-chain inits
# list) replaces the base value outright instead of being merged into it.
with_args <- function(base, ...) {
  v <- list(...)
  for (n in names(v)) base[n] <- list(v[[n]])
  base
}
# Manual mode: no convergence checking; stops once enough draws exist.
manual <-with_args(automated, check.freq = NULL, max.samples.saved = 300)

jitter <- function(chn) list(mu = rnorm(1, 0, 2), sigma = runif(1, 0.5, 2))
per.chain <- list(list(mu = -3, sigma = 1), list(mu = 3, sigma = 2))

run_case <- function(name, args) {
  args$mod.nam <- paste0("mod_", name)
  args$dump.path <- file.path(run.dir, paste0("dump_", name))
  cat("\n== Run", name, "==\n")
  t0 <- Sys.time()
  res <- tryCatch(
    R.utils::withTimeout(do.call(manageNimble::runNimble, args),
                         timeout = run.timeout, onTimeout = "error"),
    error = function(e) e)
  secs <- round(as.numeric(Sys.time() - t0, units = "secs"))
  if (inherits(res, "error")) {
    record(paste0(name, ": runNimble() completed"), FALSE,
           paste0(conditionMessage(res), " (after ", secs, " s; see ", args$mod.nam, "_logs/)"))
    return(NULL)
  }
  record(paste0(name, ": runNimble() completed"), TRUE, paste0(secs, " s"))
  R.utils::loadObject(args$mod.nam)
}

# Which top-level elements of two saved objects differ.
differing <- function(a, b) {
  nms <- union(names(a), names(b))
  nms[!vapply(nms, function(n) identical(a[[n]], b[[n]]), logical(1))]
}
check_identical <- function(label, a, b) {
  if (is.null(a) || is.null(b)) return(record(label, FALSE, "a run failed"))
  d <- differing(a, b)
  record(label, length(d) == 0, if (length(d)) paste("differs in:", paste(d, collapse = ", ")) else "")
}
check_monitors2 <- function(label, mod) {
  if (is.null(mod)) return(invisible())
  ok <- !is.null(mod$monitors2) && !is.null(mod$iter.key) &&
    nrow(mod$iter.key) == nrow(mod$mcmcOutput) &&
    nrow(merge(mod$monitors2$iter.key, mod$iter.key, by = c("chn", "iter"))) == nrow(mod$monitors2$iter.key)
  record(paste0(label, ": every monitors2 row joins to a primary row on iter.key"), ok,
         if (!is.null(mod$monitors2)) paste(nrow(mod$monitors2$out), "monitors2 rows,",
                                            nrow(mod$mcmcOutput), "primary rows") else "no monitors2")
}

# --- A: the existing usage - single inits list, no seed -----------------------
a <- run_case("A_baseline", automated)
if (!is.null(a)) {
  record("A: saved object has no new elements without seed or per-chain inits",
         !any(c("rng", "inits") %in% names(a)), paste(names(a), collapse = ", "))
  check_monitors2("A", a)
  record("A: max.samples.saved forced additional thinning (correspondence tested under a cap)",
         "INFO", paste("mcmc.info nthin =", a$mcmc.info["nthin"]))
}

# --- B/C: seed + inits function -----------------------------------------------
b.args <- with_args(automated, inits = jitter, seed = 2024)
b1 <- run_case("B_seed_fn_1", b.args)
b2 <- run_case("B_seed_fn_2", b.args)
check_identical("B: same seed -> identical saved objects", b1, b2)
if (!is.null(b1)) {
  record("B: mod$rng records the seed and one stream per chain",
         identical(b1$rng$seed, 2024L) && identical(dim(b1$rng$streams), c(as.integer(nc), 7L)))
  record("B: mod$inits records the generated inits, different per chain",
         length(b1$inits) == nc && !identical(b1$inits[[1]], b1$inits[[2]]))
  check_monitors2("B", b1)
}
c1 <- run_case("C_seed_fn_other", with_args(b.args, seed = 2025))
if (!is.null(b1) && !is.null(c1))
  record("C: different seed -> different draws", !identical(b1$mcmcOutput, c1$mcmcOutput))

# --- D: seed + per-chain list + absolute burn-in ------------------------------
d.args <- with_args(automated, inits = per.chain, nb = 100, seed = 77)
d1 <- run_case("D_seed_perchain_1", d.args)
d2 <- run_case("D_seed_perchain_2", d.args)
check_identical("D: same seed -> identical saved objects", d1, d2)
if (!is.null(d1)) {
  record("D: mod$inits is exactly the per-chain list supplied", identical(d1$inits, per.chain))
  record("D: absolute burn-in recorded", unname(d1$mcmc.info["burnin"]) == 100,
         paste("mcmc.info burnin =", d1$mcmc.info["burnin"]))
}

# --- E/F: manual mode -----------------------------------------------------------
# Workers sample freely until the parent's check sees enough draws; on a toy
# model they overshoot by many blocks. A seeded run should gather exactly the
# minimum number needed regardless.
cap <- manageNimble:::minGatherBlocks(manual$ni, manual$nb, manual$nt, manual$max.samples.saved)
e1 <- run_case("E_manual_seed_1", with_args(manual, seed = 99))
e2 <- run_case("E_manual_seed_2", with_args(manual, seed = 99))
check_identical("E: same seed -> identical saved objects (manual mode)", e1, e2)
if (!is.null(e1)) {
  record("E: seeded manual mode gathered exactly the minimum number of blocks",
         unname(e1$mcmc.info["niterations"]) == manual$ni * cap,
         paste0("niterations = ", e1$mcmc.info["niterations"], ", expected ", manual$ni * cap,
                " (", cap, " blocks)"))
  check_monitors2("E", e1)
}
f <- run_case("F_manual_noseed", manual)
if (!is.null(f)) {
  record("F: unseeded manual mode adds no new elements", !any(c("rng", "inits") %in% names(f)))
  record("F: unseeded manual mode gathers every block available, overshoot included", "INFO",
         paste0("niterations = ", f$mcmc.info["niterations"], " vs. the seeded minimum ", manual$ni * cap))
}

finish()
