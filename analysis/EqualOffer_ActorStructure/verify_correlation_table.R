#!/usr/bin/env Rscript
# The legacy table labels candidate parameter count as df_diff. Preserve it and
# export an explicitly corrected, independently checked descriptive table.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Supply actor result directory")
out <- normalizePath(args[1],winslash="/",mustWork=TRUE)
target <- file.path(out,"correlation_selection_verified.csv")
if(file.exists(target)) stop("Preserve verified export")
root <- dirname(dirname(out))
.libPaths(c(file.path(root,"renv/library/R-4.3/x86_64-w64-mingw32"),.libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
x <- read.csv(file.path(out,"correlation_selection.csv"))
names(x)[names(x)=="df_diff"] <- "candidate_npar"
x$reference_npar <- vapply(x$experiment,function(ex) {
  m <- readRDS(file.path(out,paste0(ex,"_selection_stages/step3_selected.rds")))
  as.numeric(attr(logLik(m),"df"))
},numeric(1))
x$df_difference <- x$candidate_npar-x$reference_npar
x$p_recomputed <- pchisq(x$chisq,df=x$df_difference,lower.tail=FALSE)
stopifnot(all(x$df_difference>0),all(abs(x$p_recomputed-x$p_value)<1e-12))
readr::write_csv(x,target)
writeLines(c("Verified correlation LRT p-values using candidate minus reference parameter count.",
  "Legacy df_diff field contains candidate parameter count; source engine uses the actual anova p-value for selection.",
  "Original raw table and engine preserved. This descriptive-label correction changes no fit or selection decision."),
  file.path(out,"CORRELATION_TABLE_AUDIT.txt"))
