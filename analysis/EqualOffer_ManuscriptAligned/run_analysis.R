#!/usr/bin/env Rscript
# Run from the project root: Rscript --vanilla analysis/EqualOffer_ManuscriptAligned/run_analysis.R --out=results/EqualOffer_Aligned_NEW
# Existing output directories are rejected by the canonical audit.
flags <- commandArgs(trailingOnly = TRUE)
if (any(!grepl("^--(out|root)=", flags))) stop("Allowed arguments: --out=, --root=")
file_flag <- grep("^--file=", commandArgs(), value = TRUE)
start <- if (length(file_flag)) dirname(sub("^--file=", "", file_flag[1])) else getwd()
root <- normalizePath(start, winslash = "/", mustWork = TRUE)
while (!file.exists(file.path(root, "UG_ERP_Project.Rproj"))) {
  parent <- dirname(root)
  if (identical(parent, root)) stop("Project root not found")
  root <- parent
}
code <- file.path(root, "analysis", "EqualOffer_ManuscriptAligned")
# Only the existing audit stage runs; no original or Holm-version models are refitted.
audit <- new.env(parent = globalenv())
sys.source(file.path(root, "analysis", "EqualOffer_Exploratory", "run_analysis.R"), envir = audit)
out <- audit$out
writeLines("RUNNING: manuscript-aligned selection", file.path(out, "STATUS.txt"))
required <- c("buildmer", "lme4", "lmerTest", "glmmTMB", "dplyr", "stringr", "emmeans", "car", "R.utils", "jsonlite", "digest")
if (any(!vapply(required, requireNamespace, logical(1), quietly = TRUE))) stop("Restore existing project dependencies first")
suppressPackageStartupMessages(pacman::p_load(char = required, install = FALSE))
source(file.path(code, "manuscript_engine.R"))
options(contrasts = c("contr.sum", "contr.poly"), digits = 17)
set.seed(2024)
dir.create(file.path(out, "aligned_code"))
file.copy(list.files(code, full.names = TRUE), file.path(out, "aligned_code"))
writeLines(capture.output(sessionInfo()), file.path(out, "ALIGNED_SESSION_INFO.txt"))
save_csv <- function(x, name) readr::write_csv(x, file.path(out, name), na = "NA")
provenance <- jsonlite::read_json(file.path(code, "engine_provenance.json"))
stopifnot(identical(digest::digest(file = file.path(root, "Sta_Behaviour_E1_E2_Integrative.Rmd"), algo = "sha256"), provenance$source_sha256))
stopifnot(identical(digest::digest(file = file.path(code, "manuscript_engine.R"), algo = "sha256"), provenance$snapshot_sha256))

contrasts_all <- probabilities <- omnibus <- diagnostics <- main_effects <- list()
for (ex in c("E1", "E2")) {
  message("MANUSCRIPT-ALIGNED: ", ex)
  d <- audit$data_by_exp[[ex]]
  log_file <- file(file.path(out, paste0(ex, "_aligned_selection.log")), "wt")
  sink(log_file, split = TRUE)
  selected <- tryCatch(fit_stage1_bates(d, "reject_binary", "emotion * allocation",
    c("emotion", "allocation", "emotion:allocation"), model_family = "binomial",
    random_max_vars = c("emotion", "allocation"), timeout_sec = 1800), error = identity)
  sink(); close(log_file)
  if (inherits(selected, "error")) {
    diagnostics[[ex]] <- data.frame(experiment = ex, valid = FALSE, reason = conditionMessage(selected))
    save_csv(dplyr::bind_rows(diagnostics), "aligned_diagnostics.csv")
    next
  }
  m <- selected
  stopifnot(isTRUE(mixed_model_diagnostics(m)$valid))
  saveRDS(m, file.path(out, paste0(ex, "_fair_choice.rds")))
  saveRDS(attr(m, "bates_diag"), file.path(out, paste0(ex, "_aligned_selection.rds")))
  jsonlite::write_json(attr(m, "bates_diag"), file.path(out, paste0(ex, "_aligned_selection.json")), pretty = TRUE, auto_unbox = TRUE, na = "null")
  diagnostics[[ex]] <- data.frame(experiment = ex, valid = TRUE, reason = format_model_diagnostics(mixed_model_diagnostics(m)),
    formula = paste(deparse(formula(m)), collapse = " "), optimizer = attr(m, "final_refit_optimizer"))
  save_csv(dplyr::bind_rows(diagnostics), "aligned_diagnostics.csv")
  tab <- as.data.frame(car::Anova(m, type = 3)); tab$term <- rownames(tab); rownames(tab) <- NULL
  tab$experiment <- ex; omnibus[[ex]] <- tab
  p_int <- tab$`Pr(>Chisq)`[tab$term == "emotion:allocation"]
  gate <- if (p_int < .05) "original_main_route" else if (p_int < .10) "original_supplement_route" else "original_not_triggered"
  eg <- emmeans::emmeans(m, ~ emotion * allocation)
  grid <- as.data.frame(eg)
  pairs <- combn(levels(d$emotion), 2, simplify = FALSE)
  pair_list <- setNames(pairs, vapply(pairs, paste, character(1), collapse = " - "))
  weights <- audit$pair_weights(grid, pair_list, "allocation", c("5:5", "6:4"))
  ct <- as.data.frame(summary(emmeans::contrast(eg, method = weights), infer = c(TRUE, TRUE), adjust = "none"))
  ct$experiment <- ex; ct$model <- "primary"; ct$scale <- "log_odds"
  ct$context <- sub("^.* \\| ", "", ct$contrast)
  ct$pair <- sub(" \\|.*$", "", ct$contrast)
  ct <- ct %>% dplyr::group_by(context) %>% dplyr::mutate(p_BH = p.adjust(p.value, "BH"), family_n = dplyr::n()) %>% dplyr::ungroup()
  stopifnot(all(ct$family_n == 10))
  ct$OR_or_ratio_of_OR <- exp(ct$estimate); ct$OR_lower <- exp(ct$asymp.LCL); ct$OR_upper <- exp(ct$asymp.UCL)
  ct$original_omnibus_route <- gate; ct$status <- paste0("exploratory_", gate)
  contrasts_all[[ex]] <- ct
  pr <- as.data.frame(summary(eg, type = "response", infer = c(TRUE, FALSE)))
  pr$experiment <- ex; probabilities[[ex]] <- pr
  me <- as.data.frame(summary(pairs(emmeans::emmeans(m, ~ emotion)), infer = c(TRUE, TRUE), adjust = "none"))
  me$p_BH <- p.adjust(me$p.value, "BH"); me$experiment <- ex
  main_effects[[ex]] <- me
  allocation_main <- as.data.frame(summary(pairs(emmeans::emmeans(m, ~ allocation)), infer = c(TRUE, TRUE), adjust = "tukey"))
  save_csv(allocation_main, paste0(ex, "_allocation_main.csv"))
  save_csv(dplyr::bind_rows(contrasts_all), "new_contrasts.csv")
  save_csv(dplyr::bind_rows(probabilities), "new_probabilities.csv")
  save_csv(dplyr::bind_rows(omnibus), "new_omnibus.csv")
  save_csv(dplyr::bind_rows(main_effects), "emotion_main_contrasts.csv")
  writeLines(c("Retrospective exploratory fair-only Bernoulli-logit GLMM.",
    paste("Formula:", paste(deparse(formula(m)), collapse = " ")),
    "Exact manuscript five-step selection functions; blocked LRT alpha=.20; final dual-optimizer validation.",
    "BH correction for each experiment/context's ten expression-pair tests; no cross-experiment family.",
    "All participants retained. Original response/RT exclusions preserved.",
    "All contrasts are retained for audit; standard_reporting enforces the original omnibus gate.",
    "Unadjusted 95% Wald CIs; selection uncertainty not included; predictions set random effects to zero."),
    file.path(out, paste0("Methods_paragraph_aligned_", ex, ".md")))
  rm(m, selected); gc()
}
valid <- all(vapply(diagnostics, function(x) isTRUE(x$valid), logical(1)))
writeLines(if (valid) "FITTED: independent verification pending" else "INCOMPLETE: manuscript acceptance gate failed; see aligned_diagnostics.csv", file.path(out, "STATUS.txt"))
writeLines(c("New manuscript-aligned exploratory analysis; prior Holm outputs retained in their original run.",
  "See ANALYSIS_PLAN.md and aligned_code for prospectively recorded scope of this revision.",
  "The existing_* files and original REPRODUCIBILITY.txt are legacy audit-stage artifacts only.",
  "Use new_contrasts.csv p_BH for this aligned analysis. No claim of preregistration or independent replication.",
  "Primary fits only: prior actor/deletion checks apply to the old models, not these new fits."), file.path(out, "ALIGNED_REPRODUCIBILITY.txt"))
message("Aligned output: ", out)
if (!valid) quit(status = 2)
