# The parent-side helpers behind runNimble()'s `inits` forms and `seed`
# (R/reproducibility.R). Pure R - no NIMBLE, no workers.

# Every test that touches the session RNG restores it itself, without relying
# on the helpers under test.
local_session_rng <- function(env = parent.frame()) {
  kind <- RNGkind()
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  withr::defer({
    suppressWarnings(RNGkind(kind[1], kind[2], kind[3]))
    if (is.null(seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", seed, envir = globalenv())
  }, envir = env)
}

# ---------------------------------------------------------------------------
# classifyInits()
# ---------------------------------------------------------------------------

test_that("a fully named list, list(), NULL and non-lists are a single inits list", {
  expect_equal(classifyInits(list(mu = 0, sigma = 1), nc = 3), "single")
  expect_equal(classifyInits(list(), nc = 3), "single")
  expect_equal(classifyInits(NULL, nc = 3), "single")
  expect_equal(classifyInits(c(mu = 0), nc = 3), "single")        # passed through, NIMBLE decides
  expect_equal(classifyInits(data.frame(mu = 0), nc = 3), "single")
})

test_that("an unnamed list of nc named lists is per-chain; an empty list is a valid chain", {
  expect_equal(classifyInits(list(list(mu = 0), list(mu = 1)), nc = 2), "per.chain")
  expect_equal(classifyInits(list(list(mu = 0), list()), nc = 2), "per.chain")
  expect_equal(classifyInits(list(list(mu = 0)), nc = 1), "per.chain")
})

test_that("a function is a generator", {
  expect_equal(classifyInits(function() list(mu = 0), nc = 2), "function")
})

test_that("ambiguous or malformed inits are errors", {
  expect_error(classifyInits(list(list(mu = 0), list(mu = 1)), nc = 3),
               "unnamed list of length 2.*nc = 3")
  expect_error(classifyInits(list(mu = 0, 1), nc = 2), "some named and some unnamed")
  expect_error(classifyInits(list(1, 2), nc = 2), "inits\\[\\[1\\]\\] must be a named list")
  expect_error(classifyInits(list(list(0), list(1)), nc = 2), "every element needs a name")
})

test_that("GUARD: list-valued elements in an inits list are an error, not silently ignored", {
  # nimbleModel() quietly ignores inits names it doesn't recognise, list-valued
  # ones included, so list(chain1 = ..., chain2 = ...) would otherwise run with
  # every initial value discarded.
  expect_error(classifyInits(list(chain1 = list(mu = 0), chain2 = list(mu = 1)), nc = 2),
               "list-valued element\\(s\\): chain1, chain2")
  expect_error(classifyInits(list(list(mu = 0, z = list(1)), list(mu = 1)), nc = 2),
               "inits\\[\\[1\\]\\] has list-valued element\\(s\\): z")
})

# ---------------------------------------------------------------------------
# prepareChainInits()
# ---------------------------------------------------------------------------

test_that("a single inits list passes through untouched, with nothing to record", {
  x <- list(mu = 0, sigma = matrix(1:4, 2))
  r <- prepareChainInits(x, nc = 3)
  expect_identical(r$inits, x)
  expect_false(r$per.chain)
  expect_null(r$used)
})

test_that("a per-chain list passes through as-is and is recorded", {
  x <- list(list(mu = 0), list(mu = 5))
  r <- prepareChainInits(x, nc = 2)
  expect_identical(r$inits, x)
  expect_true(r$per.chain)
  expect_identical(r$used, x)
})

test_that("a zero-argument generator is called once per chain", {
  calls <- 0
  f <- function() { calls <<- calls + 1; list(mu = calls) }
  r <- prepareChainInits(f, nc = 3)
  expect_equal(calls, 3)
  expect_true(r$per.chain)
  expect_equal(vapply(r$inits, function(i) i$mu, numeric(1)), c(1, 2, 3))
  expect_identical(r$used, r$inits)
})

test_that("a generator taking an argument receives the chain number", {
  r <- prepareChainInits(function(chn) list(mu = chn * 10), nc = 3)
  expect_equal(vapply(r$inits, function(i) i$mu, numeric(1)), c(10, 20, 30))
})

test_that("a generator's return value is validated, and its errors name the chain", {
  expect_error(prepareChainInits(function() 1, nc = 2),
               "return value for chain 1 must be a named list")
  expect_error(prepareChainInits(function(chn) if (chn == 2) list(chains = list(1)) else list(mu = 0), nc = 2),
               "return value for chain 2 has list-valued element")
  expect_error(prepareChainInits(function(chn) if (chn == 2) stop("boom") else list(mu = 0), nc = 2),
               "inits function failed for chain 2: boom")
})

test_that("seeded generated inits are reproducible, seed-dependent and differ by chain", {
  local_session_rng()
  f <- function() list(mu = rnorm(1), sigma = runif(1))
  a1 <- prepareChainInits(f, nc = 2, deriveChainStreams(123, 2))$inits
  a2 <- prepareChainInits(f, nc = 2, deriveChainStreams(123, 2))$inits
  b  <- prepareChainInits(f, nc = 2, deriveChainStreams(456, 2))$inits
  expect_identical(a1, a2)
  expect_false(identical(a1, b))
  expect_false(identical(a1[[1]], a1[[2]]))
})

test_that("generated inits use a substream independent of the chain's sampling stream", {
  local_session_rng()
  streams <- deriveChainStreams(123, 2)
  gen <- prepareChainInits(function() list(u = runif(1)), nc = 2, streams)$inits
  RNGkind(RNG_KIND[1], RNG_KIND[2], RNG_KIND[3])
  assign(".Random.seed", streams[[1]], envir = globalenv())
  expect_false(gen[[1]]$u == runif(1))   # the chain's first sampling draw is a different number
})

test_that("generating inits leaves the session's RNG exactly as it was, seeded or not", {
  local_session_rng()
  set.seed(42)
  before <- get(".Random.seed", envir = globalenv())
  prepareChainInits(function() list(mu = rnorm(1)), nc = 3)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  prepareChainInits(function() list(mu = rnorm(1)), nc = 3, deriveChainStreams(1, 3))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_equal(RNGkind()[1], "Mersenne-Twister")
})

# ---------------------------------------------------------------------------
# deriveChainStreams()
# ---------------------------------------------------------------------------

test_that("streams are L'Ecuyer-CMRG states, one per chain, deterministic and distinct", {
  local_session_rng()
  s <- deriveChainStreams(123, 4)
  expect_length(s, 4)
  expect_true(all(vapply(s, function(x) is.integer(x) && length(x) == 7, logical(1))))
  expect_identical(s, deriveChainStreams(123, 4))
  expect_equal(length(unique(s)), 4)
  expect_false(identical(s, deriveChainStreams(124, 4)))
  # the first nc - 1 streams don't depend on nc
  expect_identical(deriveChainStreams(123, 2), s[1:2])

  # element 1 encodes the generator kinds a worker adopts on assigning it
  RNGkind(RNG_KIND[1], RNG_KIND[2], RNG_KIND[3]); set.seed(1)
  expect_equal(s[[1]][1], get(".Random.seed", envir = globalenv())[1])
})

test_that("an installed stream reproduces the same draws", {
  local_session_rng()
  s <- deriveChainStreams(7, 1)[[1]]
  RNGkind(RNG_KIND[1], RNG_KIND[2], RNG_KIND[3])
  assign(".Random.seed", s, envir = globalenv()); x <- c(runif(3), rnorm(3))
  assign(".Random.seed", s, envir = globalenv()); y <- c(runif(3), rnorm(3))
  expect_identical(x, y)
})

test_that("deriving streams leaves the session's RNG kind and state exactly as they were", {
  local_session_rng()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  deriveChainStreams(123, 3)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_equal(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))
})

test_that("no .Random.seed before means none after", {
  local_session_rng()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  deriveChainStreams(123, 2)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("a session's non-default normal/sample kinds don't leak into the streams", {
  local_session_rng()
  reference <- deriveChainStreams(123, 2)
  suppressWarnings(RNGkind("Mersenne-Twister", "Box-Muller", "Rounding"))
  expect_identical(deriveChainStreams(123, 2), reference)
  expect_equal(RNGkind(), c("Mersenne-Twister", "Box-Muller", "Rounding"))
})

# ---------------------------------------------------------------------------
# minGatherBlocks()
# ---------------------------------------------------------------------------

test_that("minGatherBlocks() is the smallest block count passing manual mode's stopping rule", {
  # The rule as written in runNimble()'s manual branch (do.gather.check).
  passes <- function(nblks, ni, nb, nt, mss) {
    nb.now <- ifelse(nb<1, nb*ni*nblks, nb)
    ni.now <- ni*nblks
    (ni.now - nb.now) >= mss * nt
  }
  grid <- expand.grid(ni = c(200, 1000, 2010), nb = c(0, 0.25, 0.4, 0.5, 3000),
                      nt = c(1, 5, 10), mss = c(1, 100, 333, 1000))
  for (k in seq_len(nrow(grid))) {
    g <- grid[k, ]
    n <- minGatherBlocks(g$ni, g$nb, g$nt, g$mss)
    expect_true(passes(n, g$ni, g$nb, g$nt, g$mss))
    if (n > 1) expect_false(passes(n - 1, g$ni, g$nb, g$nt, g$mss))
  }
})
