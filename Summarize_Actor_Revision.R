# Post-hoc summary of the saved actor-sensitivity fits (no refitting).
renv_lib <- normalizePath("renv/library/R-4.3/x86_64-w64-mingw32", mustWork = TRUE)
.libPaths(c(renv_lib, .libPaths()))
options(contrasts = c("contr.sum", "contr.poly"))
suppressPackageStartupMessages({ library(glmmTMB); library(lme4); library(lmerTest); library(emmeans) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the actor revision output directory.")
out <- normalizePath(args[1], mustWork = TRUE)

bind_fill <- function(l) {
  l <- Filter(Negate(is.null), l); if (!length(l)) return(data.frame())
  nm <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) { for (n in setdiff(nm, names(x))) x[[n]] <- NA; x[, nm, drop = FALSE] }))
}
vc_table <- function(m, experiment, outcome, label) {
  vc <- if (inherits(m, "glmmTMB")) VarCorr(m)$cond else VarCorr(m)
  rows <- list()
  for (i in seq_along(vc)) {
    g <- names(vc)[i]; v <- diag(as.matrix(vc[[i]]))
    for (j in seq_along(v)) rows[[length(rows) + 1]] <- data.frame(experiment, outcome, model = label, group = g,
                                                                 term = names(v)[j], variance = v[j], sd = sqrt(v[j]))
  }
  do.call(rbind, rows)
}
contrast_tbl <- function(m, experiment, outcome, label, by_offer) {
  grid <- if (by_offer) emmeans(m, ~ emotion | offer_type, type = "link", lmerTest.limit = 1e6, lmer.df = "satterthwaite")
          else emmeans(m, ~ emotion, type = "link", lmerTest.limit = 1e6, lmer.df = "satterthwaite")
  raw <- as.data.frame(summary(pairs(grid, adjust = "none"), infer = c(TRUE, TRUE)))
  fdr <- as.data.frame(summary(pairs(grid, adjust = "fdr"), infer = c(TRUE, TRUE)))
  raw$p_fdr <- fdr$p.value; raw$CI_low_fdr <- fdr$lower.CL; raw$CI_high_fdr <- fdr$upper.CL
  names(raw)[names(raw) %in% c("z.ratio", "t.ratio")] <- "ratio"
  if (!is.null(raw$df)) raw$df <- as.numeric(raw$df)
  cbind(experiment, outcome, model = label, family = if (by_offer) "emotion_within_offer" else "emotion_marginal", raw)
}
vcs <- list(); cons <- list()
for (f in list.files(out, pattern = "^model_.*\\.rds$")) {
  p <- strsplit(sub("^model_(E[12])_(choice|rt_main)_(.*)\\.rds$", "\\1|\\2|\\3", f), "\\|")[[1]]
  m <- readRDS(file.path(out, f))
  valid <- isTRUE(attr(m, "diag")$valid)
  vcs[[f]] <- cbind(vc_table(m, p[1], p[2], p[3]), valid = valid)
  if (valid) { cons[[paste(f, "m")]] <- contrast_tbl(m, p[1], p[2], p[3], FALSE); cons[[paste(f, "o")]] <- contrast_tbl(m, p[1], p[2], p[3], TRUE) }
  rm(m); invisible(gc())
}
write.csv(bind_fill(vcs), file.path(out, "variance_components.csv"), row.names = FALSE, na = "")
write.csv(bind_fill(cons), file.path(out, "emotion_contrasts_comparison.csv"), row.names = FALSE, na = "")
message("done")
