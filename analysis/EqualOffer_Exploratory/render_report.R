#!/usr/bin/env Rscript
# Rscript --vanilla analysis/EqualOffer_Exploratory/render_report.R results/EqualOffer_full_20260928_v4
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 1) stop("Supply one verified output directory")
run <- normalizePath(a[1], winslash = "/", mustWork = TRUE)
root <- run
while (!file.exists(file.path(root, "UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if (parent == root) stop("Project root not found"); root <- parent
}
lib <- file.path(root, "renv", "library", paste0("R-", R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]), R.version$platform)
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr)})
if (!file.exists(file.path(run, "VERIFICATION.txt"))) stop("Run verify_results.R first")
if (file.exists(file.path(run, "REPORT.md"))) stop("Report already exists; not overwriting")
read <- function(f) readr::read_csv(file.path(run, f), show_col_types = FALSE)
ct <- read("new_contrasts.csv"); cells <- read("descriptive_cells.csv")
pc <- read("participant_cells.csv"); pr <- read("new_probabilities.csv")
loo <- read("participant_influence.csv"); diag <- read("sensitivity_diagnostics.csv")
primary <- ct %>% filter(model == "primary", scale == "log_odds")
valid_loo <- loo %>% filter(valid, scale == "log_odds")
# Recalculate the designated full-family Holm adjustment for each deletion:
# replace this experiment's contrasts, keep the other experiment's primary contrasts.
ids <- unique(valid_loo[c("experiment", "omitted_participant_id_internal")])
loo_adjusted <- list()
for (i in seq_len(nrow(ids))) {
  z <- valid_loo %>% filter(experiment == ids$experiment[i],
    omitted_participant_id_internal == ids$omitted_participant_id_internal[i])
  z$p_holm <- rep(NA_real_, nrow(z))
  for (fam in c("focal", "secondary")) {
    rows <- which(z$family == fam)
    other <- primary %>% filter(experiment != ids$experiment[i], family == fam)
    z$p_holm[rows] <- head(p.adjust(c(z$p.value[rows], other$p.value), method = "holm",
      n = ifelse(fam == "focal", 12L, 6L)), length(rows))
  }
  loo_adjusted[[i]] <- z
}
la <- bind_rows(loo_adjusted)
influence <- la %>% group_by(experiment, family, contrast) %>%
  summarise(valid_deletions = n(), minimum_estimate = min(estimate), maximum_estimate = max(estimate),
    maximum_p_holm = max(p_holm), minimum_p_holm = min(p_holm),
    sign_always_positive = all(estimate > 0), .groups = "drop")
write_csv(influence, file.path(run, "participant_influence_summary.csv"))
write_csv(la, file.path(run, "participant_influence_holm.csv"))
audit <- ct %>% mutate(source_file = "new_contrasts.csv", source_data_row = row_number(),
  claim = paste(experiment, model, contrast, scale, sep = " | "))
write_csv(audit, file.path(run, "stat_source_audit.csv"))

# Visualization-only update: use the manuscript's Python paired-facet workflow.
# This call reads exported CSVs and never fits or adjusts a statistical model.
python <- Sys.getenv("UG_PYTHON", unset = Sys.which("python"))
if (!nzchar(python)) stop("Python is required for manuscript-style figures; set UG_PYTHON. See README.md.")
plot_script <- file.path(root, "analysis", "EqualOffer_Exploratory", "plot_paired.py")
plot_status <- system2(python, c(shQuote(plot_script), "--run", shQuote(run)))
if (plot_status != 0L) stop("Manuscript-style plotting failed; numerical outputs are preserved.")

md_table <- function(x) {
  x[] <- lapply(x, function(z) gsub("|", "/", as.character(z), fixed = TRUE))
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
fmt <- function(x) formatC(x, digits = 3, format = "f")
format_contrasts <- function(x) data.frame(Experiment = x$experiment, Contrast = x$contrast,
  `Log odds [95% CI]` = paste0(fmt(x$estimate), " [", fmt(x$asymp.LCL), ", ", fmt(x$asymp.UCL), "]"),
  `OR or ratio of OR [95% CI]` = paste0(fmt(x$odds_ratio_or_ratio_of_OR), " [", fmt(x$OR_lower), ", ", fmt(x$OR_upper), "]"),
  `Holm p` = formatC(x$p_holm, format = "g", digits = 5), check.names = FALSE)
des <- cells %>% mutate(`Observed rejection` = paste0(rejections, "/", valid_trials, " (", fmt(100 * trial_weighted_rate), "%)"),
  `Participants rejecting` = paste0(participants_with_rejections, "/", participants_with_trials)) %>%
  select(experiment, allocation, emotion, `Observed rejection`, `Participants rejecting`)
focal <- primary %>% filter(family == "focal")
secondary <- primary %>% filter(family == "secondary")
inf_table <- influence %>% mutate(across(c(minimum_estimate, maximum_estimate), fmt),
  across(c(maximum_p_holm, minimum_p_holm), function(x) formatC(x, format = "g", digits = 5)))
old <- read("existing_omnibus.csv") %>% filter(term != "(Intercept)") %>%
  mutate(test = ifelse(outcome == "GLMM_rejection", "Wald chi-square", "Satterthwaite F"),
    statistic = ifelse(outcome == "GLMM_rejection", Chisq, `F value`),
    numerator_df = ifelse(outcome == "GLMM_rejection", Df, NumDF),
    denominator_df = ifelse(outcome == "GLMM_rejection", NA_real_, DenDF),
    p = ifelse(outcome == "GLMM_rejection", `Pr(>Chisq)`, `Pr(>F)`)) %>%
  select(experiment, outcome, term, test, statistic, numerator_df, denominator_df, p) %>%
  mutate(across(c(statistic, numerator_df, denominator_df), fmt), p = formatC(p, format = "g", digits = 5))
invalid <- loo %>% filter(!valid) %>% select(any_of(c("experiment", "omitted_participant_id_internal", "reason"))) %>% distinct()
sim_lines <- unlist(lapply(c("E1", "E2"), function(ex) {
  z <- read(paste0(ex, "_parametric_model_check.csv"))
  flagged <- z$cell[z$observed_rejections < z$simulated_lower | z$observed_rejections > z$simulated_upper]
  paste0("- ", ex, ": observed summaries outside the 200-simulation central 95% ranges: ",
    ifelse(length(flagged), paste(flagged, collapse = ", "), "none"),
    ". These checks condition on fitted parameters and are not formal goodness-of-fit tests.")
}))
report <- c(
  "# Equal-offer exploratory analysis: results and verification", "",
  "## Material Passport", "", "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate", "- Origin Date: 2026-09-28", "- Verification Status: VERIFIED (computational scope only)",
  "- Human verification: false", "- Analysis status: retrospective exploratory", "",
  "## Scope", "",
  "The original manuscript analyses were preserved. Canonical trial inputs were fingerprint-matched to the explicit formal runs and matched to saved choice model frames in order. New fair-only models separate 5:5 from 6:4, retaining all expressions and all participants. See README.md in the script folder for the frozen contrasts and fitting rules.", "",
  "These data are preprocessing outputs. Trial counts do not reconstruct exclusions upstream of that export. All probabilities below must be distinguished from subjective WTR, AV, status threat, and intentions, which were not measured by this analysis.", "",
  "## Existing manuscript models: separate choice and RT evidence", "", md_table(old), "",
  "Source: existing_omnibus.csv and existing_model_sources.csv. Omnibus tests do not by themselves identify which expression contrast varies with pooled fairness. Effect estimates and intervals for the specified expression comparisons are in existing_targeted_contrasts.csv. No pooled-experiment model was refitted or substituted.", "",
  "## Observed data", "", md_table(des), "",
  "Source: descriptive_cells.csv; locate rows by experiment, allocation and emotion. Participant-level distribution: participant_cells.csv. Missing cells and zero-rejection participants are retained.", "",
  "![Observed participant pairs and saved model estimates](paired_figures/equal_offers_paired.png)", "",
  "Thin lines connect each participant's 5:5 and 6:4 observed rates. Diamonds show saved model predictions and unadjusted 95% CIs; dashed lines mark Neutral-6:4 predictions. See paired_figures/caption.txt for the full legend. Palette and paired-facet design match the manuscript.", "",
  "## Focal exploratory contrasts", "", md_table(format_contrasts(focal)), "",
  "Holm correction covers all 12 focal log-odds tests across both experiments. Intervals are unadjusted 95% Wald intervals. Interaction orientation is 6:4 minus 5:5; exponentiated interaction estimates are ratios of ORs.", "",
  "## Secondary Disgust-Reward contrasts", "", md_table(format_contrasts(secondary)), "",
  "These six log-odds tests form a separately disclosed Holm family. Secondary does not mean preregistered or theoretically unimportant.", "",
  "## Model probability estimates", "", "Model predictions appear as diamonds in the paired figure above.", "",
  "Random effects are set to zero. These are conditional model predictions, not integrated population means and not the observed percentages above. Probability differences and their unadjusted intervals are in new_contrasts.csv (scale=probability_difference); their p-values are not the designated tests.", "",
  "## Diagnostics and sensitivity", "", md_table(diag), "",
  "Actor-intercept sensitivity estimates are fully reported in new_contrasts.csv. Numerical validity is necessary but does not establish the adequacy of normal random effects or precision with sparse rejection events.", "",
  "### Participant influence", "", md_table(inf_table), "",
  "Each valid deletion uses the selected participant formula unchanged. For its adjusted p-values, the other experiment's primary contrasts remain in the same full family. No deletion is adopted as the primary analysis. Changes in significance are diagnostics, not grounds for excluding participants.", "",
  if (nrow(invalid)) c("Invalid or boundary deletion fits (not interpreted):", "", md_table(invalid), "") else
    "All participant-deletion fits passed the recorded numerical and random-boundary checks.", "",
  "### Parametric model checks", "", sim_lines, "",
  "## Independent verification", "",
  "verify_results.R independently rebuilt trial counts and IDs using base R, verified every participant-condition count, recomputed contrasts and standard errors by direct model-matrix algebra, checked confidence intervals and Holm adjustments, and refitted the selected primary models from canonical inputs. See verification_checks.csv and independent_refit_comparison.csv. Exact reproduction on this machine does not guarantee bitwise equality on a different software stack.", "",
  "## Statistical fallacy scan: 11/11 assessed", "",
  "| Check | Assessment within this analysis |", "|---|---|",
  "| Simpson's paradox | Experiments and allocations are reported separately; pooled descriptions are not used as conditional effect estimates. |",
  "| Ecological inference | Trial models account for participants; aggregate rates do not establish person-level motives. |",
  "| Selection/Berkson bias | Preprocessing selection is inherited and reported; excluded upstream trials are not reconstructed. |",
  "| Collider adjustment | No new post-treatment covariate adjustment; original RT filtering remains a design limitation. |",
  "| Base-rate neglect | Event counts and absolute probabilities accompany odds ratios. |",
  "| Regression to mean | No extreme-score selection or pre/post improvement claim. |",
  "| Survivorship bias | All available participants retained, including zero rejecters; upstream selection remains explicit. |",
  "| Look-elsewhere effect | Complete focal/secondary families and Holm corrections reported. |",
  "| Forking paths | Exploratory status, candidate reduction rules, code snapshots and all numerical diagnostics disclosed. |",
  "| Correlation/causation | No claim that these contrasts identify an unmeasured psychological mechanism. |",
  "| Reverse causality | No causal chain among unmeasured appraisals is inferred. |", "",
  "The skill's generic effect-size labels and significance-based confidence categories are not used: the author's instructions require context-specific interpretation, uncertainty and a distinction between numerical verification and scientific reliability.", "",
  "## Interpretation boundaries", "",
  "Expression differences at 5:5 would require an account beyond monetary inequality alone. They do not establish favorable subjective WTR or prove status defense, intention reciprocity, or AV. Nonsignificant allocation contrasts do not establish equivalence or additivity. No theoretical hierarchy or manuscript wording was changed.", "",
  "## Provenance and limitations", "",
  "All reported new statistics map to new_contrasts.csv and stat_source_audit.csv by experiment, model, contrast and scale. Existing choice/RT evidence is in existing_omnibus.csv and existing_targeted_contrasts.csv. Software, fingerprints, model selection and formulas are preserved in SESSION_INFO.txt, input_manifest.csv, selection logs and Methods_paragraph_E*.md.", "",
  "Wald inference assumes adequate likelihood approximation. Sparse events, random-effect distribution assumptions and random-structure selection uncertainty remain limitations despite convergence and exact refits. No formal equivalence test or preregistered theory test was performed.")
writeLines(report, file.path(run, "REPORT.md"))
file.copy(file.path(root, "analysis", "EqualOffer_Exploratory", c("verify_results.R", "render_report.R", "README.md", "plot_paired.py", "manuscript_style.json")), run)
message("Report and diagnostic figures saved in ", run)
