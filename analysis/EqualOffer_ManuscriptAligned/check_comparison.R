#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop("Supply one completed result directory")
run <- normalizePath(args[1],winslash="/",mustWork=TRUE)
x <- read.csv(file.path(run,"version_comparison.csv"),check.names=FALSE)
for(family in unique(x$family)) {
  z <- x[x$family==family,]
  stopifnot(max(abs(p.adjust(z$new_p_raw,"holm")-z$new_model_p_legacy_Holm))<1e-12)
}
y <- read.csv(file.path(run,"new_contrasts.csv"),check.names=FALSE)
for(ex in unique(y$experiment)) for(context in unique(y$context)) {
  z <- y[y$experiment==ex & y$context==context,]
  stopifnot(max(abs(p.adjust(z$p.value,"fdr")-z$p_BH))<1e-12)
}
writeLines("PASS: Python diagnostic Holm bridge matches R p.adjust; all BH families match R fdr alias.",
           file.path(run,"comparison_verification.txt"))
cat("Comparison correction checks passed\n")
