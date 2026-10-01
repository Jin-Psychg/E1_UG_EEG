# Actor-sensitivity reanalysis of the confirmatory offer-locked ERP models (Round 52, R0.1)
# Read-only inputs from C:/Code/UG_ERP_Project (2026-08-22 Alday primary run); outputs under draft/sensitivity/.
renv_lib <- "C:/Code/UG_ERP_Project/renv/library/R-4.3/x86_64-w64-mingw32"
.libPaths(c(renv_lib, .libPaths()))
options(contrasts = c("contr.sum", "contr.poly"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(emmeans) })
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a
project_root <- "C:/Code/UG_ERP_Project"
out <- "E:/MyPaper_withARS/draft/sensitivity/actor-random-intercept_2026-08-23"
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000), calc.derivs = TRUE)

diag_lmm <- function(m) {
  singular <- lme4::isSingular(m, tol = 1e-4)
  code <- m@optinfo$conv$opt %||% NA_integer_
  msgs <- m@optinfo$conv$lme4$messages %||% character(0)
  H <- m@optinfo$derivs$Hessian
  pd <- if (!is.null(H)) { ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values; all(is.finite(ev)) && min(ev) > 0 } else FALSE
  hard <- "failed to converge|unable to evaluate scaled gradient|degenerate.*Hessian|negative eigenvalue"
  list(valid = isTRUE(code == 0) && !any(grepl(hard, msgs, ignore.case = TRUE)) && !singular &&
         all(is.finite(sqrt(diag(vcov(m))))) && is.finite(AIC(m)) && pd,
       code = code, messages = paste(msgs, collapse = " | "), singular = singular, pd = pd)
}
read_formula <- function(path) {
  x <- readLines(path, warn = FALSE); i <- grep("^## Final formula:", x)
  f <- x[i + 1]; j <- i + 2
  while (j <= length(x) && nchar(trimws(x[j])) > 0 && !grepl("^##", x[j])) { f <- paste(f, x[j]); j <- j + 1 }
  as.formula(f)
}
prep <- function(experiment, comp) {
  path <- file.path(project_root, "data", paste0("02_Pipeline_Output_", experiment), "Method_Regression", "Stimulus_Locked", "trials.csv")
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$Offers_Other %in% c(5, 6, 8, 9) & d$reaction != 0 & d$RT >= 300 & d$RT <= 3000, ]
  d$offer_type <- factor(ifelse(d$Offers_Other %in% c(5, 6), "fair", "unfair"), levels = c("fair", "unfair"))
  d$emotion <- factor(d$emotion, levels = c("neu", "aff", "dis", "dom", "enj"))
  d <- d[!is.na(d[[comp]]) & abs(d[[comp]]) < 200, ]
  d$Baseline_Raw <- d[[paste0("Baseline_", comp)]]
  d$Baseline_c <- d$Baseline_Raw - mean(d$Baseline_Raw, na.rm = TRUE)
  d <- d[!is.na(d$Baseline_c), ]
  d$participant_id <- factor(d$participant_id)
  d$actor_id <- factor(sub("^[LR]_(neu|aff|dis|dom|enj)", "", d$stim, perl = TRUE))
  stopifnot(nlevels(d$actor_id) == 60L)
  for (rv in c("emotion", "offer_type")) {
    mm <- model.matrix(reformulate(rv), d, contrasts.arg = setNames(list(contr.sum(nlevels(d[[rv]]))), rv))[, -1, drop = FALSE]
    for (j in seq_len(ncol(mm))) d[[sprintf("RE_%s_%d", rv, j)]] <- mm[, j]
  }
  d
}
dg <- list(); om <- list(); fu <- list()
for (experiment in c("E1", "E2")) for (comp in c("FRN_pre", "LPP_pre")) {
  rep <- file.path(project_root, "results/EEG", paste0(experiment, "_TwoStage_Bates_Alday/Stage1_TrialLevel"), comp, "REPRODUCIBILITY.txt")
  fo <- read_formula(rep); d <- prep(experiment, comp)
  for (label in c("canonical_refit", "actor_intercept")) {
    f <- if (label == "canonical_refit") fo else update(fo, . ~ . + (1 | actor_id))
    m <- lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl); di <- diag_lmm(m)
    av <- as.numeric(if ("actor_id" %in% names(VarCorr(m))) diag(as.matrix(VarCorr(m)$actor_id))[1] else NA)
    dg[[length(dg) + 1]] <- data.frame(experiment, component = comp, model = label, n_obs = nrow(d), valid = di$valid, code = di$code,
                                       messages = di$messages, singular = di$singular, pd_hessian = di$pd, logLik = as.numeric(logLik(m)),
                                       actor_intercept_var = av, residual_var = sigma(m)^2, formula = paste(deparse(f), collapse = " "))
    message(sprintf("%s %s %s valid=%s actor_var=%s", experiment, comp, label, di$valid, format(av, digits = 4)))
    if (!di$valid) next
    a <- as.data.frame(anova(m, type = 3, ddf = "Satterthwaite"))
    om[[length(om) + 1]] <- data.frame(experiment, component = comp, model = label, term = rownames(a), F = a$`F value`,
                                       df1 = a$NumDF, df2 = a$DenDF, p = a$`Pr(>F)`, row.names = NULL)
    em <- emmeans(m, ~ offer_type, at = list(Baseline_c = 0), lmer.df = "satterthwaite", lmerTest.limit = 1e6)
    cc <- as.data.frame(summary(contrast(em, method = list("fair - unfair" = c(1, -1))), infer = c(TRUE, TRUE)))
    fu[[length(fu) + 1]] <- cbind(experiment, component = comp, model = label, cc)
  }
}
write.csv(do.call(rbind, dg), file.path(out, "erp_diagnostics.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, om), file.path(out, "erp_omnibus_comparison.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, fu), file.path(out, "erp_fair_minus_unfair_comparison.csv"), row.names = FALSE, na = "")
message("done")
