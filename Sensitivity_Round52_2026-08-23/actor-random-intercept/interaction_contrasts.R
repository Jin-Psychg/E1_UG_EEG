.libPaths(c("C:/Code/UG_ERP_Project/renv/library/R-4.3/x86_64-w64-mingw32", .libPaths()))
options(contrasts = c("contr.sum", "contr.poly"))
suppressPackageStartupMessages({library(glmmTMB); library(emmeans)})
out <- list()
for (e in c("E1","E2")) {
  m <- readRDS(sprintf("model_%s_choice_canonical_refit.rds", e))
  em <- emmeans(m, ~ emotion * offer_type, type = "link")
  ic <- as.data.frame(summary(contrast(em, interaction = "pairwise", adjust = "none"), infer = c(TRUE, TRUE)))
  ic$experiment <- e; out[[e]] <- ic
}
res <- do.call(rbind, out)
write.csv(res, "choice_interaction_contrasts_canonical.csv", row.names = FALSE)
print(res[res$emotion_pairwise %in% c("dis - dom","aff - dom","dom - enj","aff - enj"), c("experiment","emotion_pairwise","offer_type_pairwise","estimate","SE","asymp.LCL","asymp.UCL","p.value")], digits = 3)
