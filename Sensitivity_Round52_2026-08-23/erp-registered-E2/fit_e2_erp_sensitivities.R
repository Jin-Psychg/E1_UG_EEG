# Experiment 2 offer-locked ERP sensitivities for Round 52 (roadmap R1.5, R1.6)
#
# Inputs (read only): e2_epoch_features.csv (extracted by extract_e2_epoch_features.py from the
# canonical single-trial epoch files) and the canonical trials.csv / final-model REPRODUCIBILITY
# formulas of the live analysis tree C:/Code/UG_ERP_Project (2026-08-22 run).
#
# Part A (R1.6): the confirmatory P3-LPP (LPP_pre, 400-600 ms) and early fronto-central (FRN_pre,
#   192-244 ms) models refitted with the canonical final formula on (i) the canonical trial set
#   (validation against the 2026-08-22 outputs) and (ii) the subset of epochs that also pass the
#   preregistered +/-100 uV peak-to-peak threshold.
# Part B (R1.5): the Experiment 2 registered windows N400 (350-450 ms; Fz Cz CPz Pz), LPP 500-800 ms
#   (Cz C1 C2 CP1 CP2 Pz) and N2 (200-300 ms; ROI [VALUE NEEDED], provisional FRN ROI) fitted with the
#   Alday baseline-as-covariate fixed structure and the canonical final random structure of the
#   matching confirmatory component (LPP_pre structure for the two late windows; FRN_pre structure
#   for N2). This is a sensitivity analysis, not the canonical buildmer-selected pipeline.
# Acceptance gate = pipeline mixed_model_diagnostics (code 0, PD Hessian, finite SE/likelihood,
# nonsingular). Outputs only under draft/sensitivity/erp-registered-E2_2026-08-23/.

renv_lib <- "C:/Code/UG_ERP_Project/renv/library/R-4.3/x86_64-w64-mingw32"
.libPaths(c(renv_lib, .libPaths()))
options(contrasts = c("contr.sum", "contr.poly"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(emmeans) })
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

out <- "E:/MyPaper_withARS/draft/sensitivity/erp-registered-E2_2026-08-23"
feat <- read.csv(file.path(out, "e2_epoch_features.csv"), stringsAsFactors = FALSE)
trials <- read.csv("C:/Code/UG_ERP_Project/data/02_Pipeline_Output_E2/Method_Regression/Stimulus_Locked/trials.csv",
                   stringsAsFactors = FALSE)
stopifnot(nrow(feat) == nrow(trials))
feat$ptp_max_uV <- feat$ptp_max_uV_excl_A2   # validated channel set (A2 is not an EEG-typed channel; validation.txt)
d <- merge(trials, feat[, c("participant_id", "index", "ptp_max_uV", "N2_reg", "N400_reg", "LPP_reg",
                            "Baseline_N2_reg", "Baseline_N400_reg", "Baseline_LPP_reg", "LPP_pre_recomputed")],
           by = c("participant_id", "index"))
stopifnot(nrow(d) == nrow(trials))

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000), calc.derivs = TRUE)

diag_lmm <- function(m) {
  singular <- lme4::isSingular(m, tol = 1e-4)
  code <- m@optinfo$conv$opt %||% NA_integer_
  msgs <- m@optinfo$conv$lme4$messages %||% character(0)
  H <- m@optinfo$derivs$Hessian
  pd <- if (!is.null(H)) { ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values; all(is.finite(ev)) && min(ev) > 0 } else FALSE
  se_ok <- all(is.finite(sqrt(diag(vcov(m)))))
  ll_ok <- is.finite(AIC(m))
  hard <- "failed to converge|unable to evaluate scaled gradient|degenerate.*Hessian|negative eigenvalue"
  list(valid = isTRUE(code == 0) && !any(grepl(hard, msgs, ignore.case = TRUE)) && !singular && se_ok && ll_ok && pd,
       code = code, messages = paste(msgs, collapse = " | "), singular = singular, pd_hessian = pd)
}

prep <- function(d, dv, baseline_col, subset_100 = FALSE) {
  x <- d[d$Offers_Other %in% c(5, 6, 8, 9) & d$reaction != 0 & d$RT >= 300 & d$RT <= 3000, ]
  x$offer_type <- factor(ifelse(x$Offers_Other %in% c(5, 6), "fair", "unfair"), levels = c("fair", "unfair"))
  x$emotion <- factor(x$emotion, levels = c("neu", "aff", "dis", "dom", "enj"))
  x <- x[!is.na(x[[dv]]) & abs(x[[dv]]) < 200, ]            # pipeline QC filter (also drops rejected epochs)
  x <- x[!is.na(x$LPP_pre), ]                                # canonical retained-epoch set (200 uV gate)
  if (subset_100) x <- x[x$ptp_max_uV <= 100, ]
  x$Baseline_Raw <- x[[baseline_col]]
  x$Baseline_c <- x$Baseline_Raw - mean(x$Baseline_Raw, na.rm = TRUE)
  x <- x[!is.na(x$Baseline_c), ]
  x$participant_id <- factor(x$participant_id)
  mm <- model.matrix(~ offer_type, x, contrasts.arg = list(offer_type = contr.sum(2)))[, -1, drop = FALSE]
  x$RE_offer_type_1 <- mm[, 1]
  x
}

F_LPP <- "%s ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + emotion:Baseline_c + offer_type:Baseline_c + (1 | participant_id) + (0 + RE_offer_type_1 | participant_id)"
F_FRN <- "%s ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + (1 | participant_id)"

specs <- list(
  list(label = "LPP_pre canonical set (validation)",  dv = "LPP_pre",  bl = "Baseline_LPP_pre",  f = F_LPP, s100 = FALSE),
  list(label = "LPP_pre +/-100 uV subset (R1.6)",     dv = "LPP_pre",  bl = "Baseline_LPP_pre",  f = F_LPP, s100 = TRUE),
  list(label = "FRN_pre canonical set (validation)",  dv = "FRN_pre",  bl = "Baseline_FRN_pre",  f = F_FRN, s100 = FALSE),
  list(label = "FRN_pre +/-100 uV subset (R1.6)",     dv = "FRN_pre",  bl = "Baseline_FRN_pre",  f = F_FRN, s100 = TRUE),
  list(label = "LPP_reg 500-800 ms registered (R1.5)", dv = "LPP_reg",  bl = "Baseline_LPP_reg",  f = F_LPP, s100 = FALSE),
  list(label = "N400_reg 350-450 ms registered (R1.5)", dv = "N400_reg", bl = "Baseline_N400_reg", f = F_LPP, s100 = FALSE),
  list(label = "N2_reg 200-300 ms registered, provisional ROI (R1.5)", dv = "N2_reg", bl = "Baseline_N2_reg", f = F_FRN, s100 = FALSE),
  list(label = "LPP_reg 500-800 ms, +/-100 uV subset", dv = "LPP_reg", bl = "Baseline_LPP_reg", f = F_LPP, s100 = TRUE)
)

omni <- list(); con <- list(); dg <- list(); emo <- list()
for (s in specs) {
  x <- prep(d, s$dv, s$bl, s$s100)
  fo <- as.formula(sprintf(s$f, s$dv))
  m <- lmerTest::lmer(fo, data = x, REML = TRUE, control = ctrl)
  di <- diag_lmm(m)
  dg[[length(dg) + 1]] <- data.frame(analysis = s$label, dv = s$dv, subset100 = s$s100, n_trials = nrow(x),
                                     n_participants = nlevels(x$participant_id), valid = di$valid, code = di$code,
                                     messages = di$messages, singular = di$singular, pd_hessian = di$pd_hessian,
                                     formula = paste(deparse(fo), collapse = " "))
  message(sprintf("%s: n=%d valid=%s", s$label, nrow(x), di$valid))
  if (!di$valid) next
  a <- as.data.frame(anova(m, type = 3, ddf = "Satterthwaite"))
  omni[[length(omni) + 1]] <- data.frame(analysis = s$label, dv = s$dv, subset100 = s$s100, term = rownames(a),
                                         F = a$`F value`, df1 = a$NumDF, df2 = a$DenDF, p = a$`Pr(>F)`, row.names = NULL)
  em <- emmeans(m, ~ offer_type, at = list(Baseline_c = 0), lmer.df = "satterthwaite", lmerTest.limit = 1e6)
  cc <- as.data.frame(summary(contrast(em, method = list("fair - unfair" = c(1, -1))), infer = c(TRUE, TRUE)))
  con[[length(con) + 1]] <- cbind(analysis = s$label, dv = s$dv, subset100 = s$s100, cc)
  ee <- emmeans(m, ~ emotion, at = list(Baseline_c = 0), lmer.df = "satterthwaite", lmerTest.limit = 1e6)
  pe <- as.data.frame(summary(pairs(ee, adjust = "fdr"), infer = c(TRUE, TRUE)))
  emo[[length(emo) + 1]] <- cbind(analysis = s$label, dv = s$dv, subset100 = s$s100, pe)
}
bind <- function(l) if (length(l)) do.call(rbind, l) else data.frame()
write.csv(bind(dg), file.path(out, "erp_sensitivity_diagnostics.csv"), row.names = FALSE, na = "")
write.csv(bind(omni), file.path(out, "erp_sensitivity_omnibus.csv"), row.names = FALSE, na = "")
write.csv(bind(con), file.path(out, "erp_sensitivity_fair_minus_unfair.csv"), row.names = FALSE, na = "")
write.csv(bind(emo), file.path(out, "erp_sensitivity_emotion_pairs_fdr.csv"), row.names = FALSE, na = "")
writeLines(c(
  "E2 offer-locked ERP sensitivities (Round 52, R1.5/R1.6)",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Input features: e2_epoch_features.csv (see extract_e2_epoch_features.py; canonical epoch files, read only).",
  "Trial filters: Offers_Other in {5,6,8,9}; reaction != 0; RT in [300,3000]; |amplitude| < 200 uV; canonical retained epochs (LPP_pre non-missing).",
  "+/-100 uV subset: epochs whose maximum EEG-channel peak-to-peak amplitude over the full -0.5..1.5 s epoch is <= 100 uV.",
  "Baseline_c: grand-mean-centred -200..0 ms ROI mean of the analysed trials (pipeline logic).",
  "Random structures: canonical 2026-08-22 final structures of LPP_pre and FRN_pre (not re-selected by buildmer).",
  "Inference: Type III Satterthwaite F; fair-unfair contrast and FDR emotion pairs at Baseline_c = 0.",
  paste0("R ", R.version.string, "; lme4 ", packageVersion("lme4"), "; lmerTest ", packageVersion("lmerTest"), "; emmeans ", packageVersion("emmeans"))
), file.path(out, "REPRODUCIBILITY.txt"))
message("done")
