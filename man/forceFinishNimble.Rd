\name{forceFinishNimble}
\alias{forceFinishNimble}
\title{
  Manually finish a slow-to-converge \code{runNimble} run from its partial block dumps
}
\description{
  Rescues a \code{\link{runNimble}} run that has not (yet) met its convergence targets, or that
  you don't want to wait on any longer, by gathering whatever sample blocks it has already
  checkpointed to \code{dump.path} and applying a manually-specified burn-in. This is safe to do
  while the run is still live: each block is written once via an atomic temp-file-then-rename (see
  \code{\link{runNimbleBlock}}) and never touched again, so there is nothing to race against when
  only reading. In preview mode (\code{commit = FALSE}, the default) nothing is written anywhere -
  it just reports Rhat/n.effective for the burn-in you specified and plots the result, so you can
  try a few burn-in values, including against a still-growing set of blocks, before deciding.  In
  commit mode (\code{commit = TRUE}) it assembles and (by default) saves the same list structure
  \code{runNimble} itself produces, so the result can be used anywhere a normal \code{runNimble}
  output would be.

  ***Note: like \code{\link{runNimble}}, the live-process safety check below relies on \code{pgrep}
  and so is only meaningful on Linux; it degrades to a warning (and no check) elsewhere.
}
\usage{
  forceFinishNimble(dump.path, burnin, mod.nam,
                    par.ignore = c(), par.dontign = c(),
                    par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                    Rht.required = 1.1, neff.required = 100,
                    max.samples.saved = NULL,
                    commit = FALSE, sav.model = TRUE,
                    delete.blocks = FALSE, plot.result = TRUE,
                    check.live = TRUE, quiet = FALSE)
}
\arguments{
  \item{dump.path}{Directory containing the sample blocks and \code{NimbleObjects.RData} written
    by an in-progress or abandoned \code{runNimble} call (its own \code{dump.path} argument).
    \code{ni} and \code{nt} are read back from \code{NimbleObjects.RData} so block files are always
    interpreted consistently with how they were written; they are not separate arguments here.}
  \item{burnin}{If < 1, proportion of the chain accumulated \emph{so far} to discard as burn-in. If
    >= 1, an absolute number of iterations to discard from the start of each chain. Same convention
    as \code{nb} in \code{\link{runNimble}}. Start small with \code{commit = FALSE} and increase
    until the diagnostics/plot look acceptable.}
  \item{mod.nam}{Character string used as the file label for the saved model object (when
    \code{sav.model = TRUE}) and to label diagnostic CSVs \code{\link{checkNimble}} writes on a
    fatal (NA or Inf/-Inf Rhat or n.effective) stop.}
  \item{par.ignore}{As in \code{\link{runNimble}}/\code{\link{checkNimble}}: parameters to exclude
    from convergence checking (unanchored substring match).}
  \item{par.dontign}{As in \code{\link{runNimble}}/\code{\link{checkNimble}}: rescue list for
    \code{par.ignore}'s substring over-matching.}
  \item{par.fuzzy.track}{As in \code{\link{runNimble}}/\code{\link{checkNimble}}: parameters
    governed by \code{fuzzy.threshold} instead of the strict Rhat/n.effective requirement.}
  \item{fuzzy.threshold}{As in \code{\link{runNimble}}/\code{\link{checkNimble}}.}
  \item{Rht.required}{Maximum Rhat allowed for \code{result} (in preview) or an unconditioned
    \code{mod} (in commit) to be reported as converged.}
  \item{neff.required}{Minimum n.effective required for the same.}
  \item{max.samples.saved}{Per-chain cap on retained draws, applied the same way as in
    \code{\link{runNimble}}. \code{NULL} (default) keeps everything the burn-in leaves behind.}
  \item{commit}{If \code{FALSE} (default), only reports diagnostics and plots the gathered chains -
    nothing is written to disk. If \code{TRUE}, assembles the final model object and (depending on
    \code{sav.model}) saves it, exactly as \code{runNimble} would have on giving up or converging.}
  \item{sav.model}{If \code{commit = TRUE}, whether to save the assembled model object to disk under
    \code{mod.nam} via \code{R.utils::saveObject}. The object is always returned by the function
    regardless of this setting, so \code{sav.model = FALSE} is only useful when you intend to save
    or use the returned value yourself.}
  \item{delete.blocks}{If \code{commit = TRUE}, whether to delete \code{dump.path} afterward.
    Default \code{FALSE}, so you can re-derive with a different burn-in later if needed.}
  \item{plot.result}{Whether to call \code{plot()} on the gathered \code{mcmcOutput} object (trace/
    density plots) before returning. Default \code{TRUE}.}
  \item{check.live}{Whether to check for still-running processes referencing \code{dump.path}
    (via \code{pgrep}) before allowing a commit. Informational only when \code{commit = FALSE};
    a hard \code{stop()} when \code{commit = TRUE} and a live process is found - committing while
    the run is still writing new blocks risks racing its own eventual save. Set to \code{FALSE} to
    skip this check entirely (not recommended; only meaningful if you have independently confirmed
    the run has stopped).}
  \item{quiet}{If \code{TRUE}, suppresses the progress/diagnostic \code{message()} calls (the
    returned value and any \code{stop()}/\code{warning()} are unaffected).}
}
\details{
  The live-process check only recognizes processes whose command line contains \code{dump.path}
  itself - this covers \code{runNimble}'s own GNU \code{parallel} supervisor and the
  \code{Rscript ModRunScript.R} chain workers it launches, since both are invoked with
  \code{dump.path} in their arguments. \code{dump.path} is matched as a literal string (its
  characters are escaped before being handed to \code{pgrep}), so punctuation in the path - "(",
  "+", "[", and similar - cannot cause it to be misread as a pattern. It has no knowledge of any
  outer process a caller's own pipeline may have wrapped around the original \code{runNimble} call
  (e.g. a per-species driver script that does other work before/after fitting); if your pipeline
  needs that process to have exited too before it will proceed, check for it separately.

  Only \code{ni} and \code{nt} are recovered automatically from \code{dump.path/NimbleObjects.RData}.
  Everything else (\code{Rht.required}, \code{neff.required}, \code{max.samples.saved}, the
  \code{par.*} arguments) is not saved there and must be supplied to match what the original
  \code{runNimble} call used, if comparable convergence behavior is wanted.

  When \code{quiet = FALSE}, the printed "Max Rhat"/"Min neff" diagnostic is computed the same way,
  and over the same set of parameters (\code{par.ignore}'d parameters excluded,
  \code{par.fuzzy.track}'d parameters included), as \code{runNimble}'s own equivalent status
  message; when \code{par.fuzzy.track} is non-empty it likewise additionally reports the proportion
  of fuzzy-tracked parameters not converged.
}
\value{
  If \code{commit = FALSE}, invisibly, a list with \code{result} (logical), \code{summary} (the
  parameter summary table), \code{mcmcOutput}, \code{nblks}, and \code{additional.thin.rate}.

  If \code{commit = TRUE}, a list with the same structure \code{runNimble} itself produces:
  \code{mcmcOutput}, \code{summary}, \code{mcmc.info}, and \code{warning} (present only if
  \code{result} was \code{FALSE}, i.e. the forced burn-in still did not meet
  \code{Rht.required}/\code{neff.required}).
}
\author{
  Quresh S. Latif and Bryan L. Nuse, Bird Conservancy of the Rockies.

  This function was originally written by Claude (Anthropic) in a Cowork session, generalizing
  an ad hoc project-specific script; see \code{NEWS.md}.
}

\keyword{nimble}
\keyword{mcmcOutput}
\keyword{runNimble}
