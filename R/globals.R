# Column names referenced via dplyr's non-standard evaluation (e.g. inside
# filter()/pull() calls on a data frame column) look like undefined globals to
# R CMD check's codetools pass, because it can't see that they resolve to data
# frame columns at the point of use. Not a bug - see Writing R Extensions,
# "Non-standard evaluation".
#
# `samp` is different: it is created by load()-ing a block dump file inside
# gatherNimble()/gatherNimble2(), which codetools can't trace either.
utils::globalVariables(c("Parameter", "Rhat", "f", "n.eff", "samp"))
