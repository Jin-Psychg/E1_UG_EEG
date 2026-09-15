# E1 choice actor-sensitivity update after the random-structure revision.
# Adapted from Sensitivity_Round52_2026-08-23/actor-random-intercept/
# run_actor_sensitivity_behaviour.R, retaining its fit and acceptance policies.
# Fresh source data are checked against the revised primary model frame.
# Outputs use a new results/ActorRevision_E1_<timestamp>/ directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the formal revised run directory relative to this project.")
project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "UG_ERP_Project.Rproj"))) stop("Run from the project root.")
reference_root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
renv_lib <- file.path(project_root, "renv/library/R-4.3/x86_64-w64-mingw32")
.libPaths(c(renv_lib, .libPaths()))
options(contrasts = c("contr.sum", "contr.poly"))
suppressPackageStartupMessages({
  library(glmmTMB); library(lme4); library(lmerTest); library(emmeans); library(car)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a
output_root <- file.path(project_root, "results", paste0("ActorRevision_E1_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(output_root)) stop("Output already exists.")
dir.create(output_root, recursive = TRUE)
cat("RUN_ROOT:", output_root, "\n")

expected_md5 <- c(
  E1 = "154ae95b7a090c5a79074e197d34f9ac",   # REPRODUCIBILITY.txt, E1 GLMM_Rejection (2026-08-22)
  E2 = "7ebd8bab4d455856a0226d1d0c73307a"    # REPRODUCIBILITY.txt, E2 GLMM_Rejection (2026-08-22)
)
rt_lower_ms <- 300; rt_upper_ms <- 3000

custom_control_glmmTMB <- glmmTMB::glmmTMBControl(
  optimizer = optim, optArgs = list(method = "BFGS"), parallel = 1)
default_control_glmmTMB <- glmmTMB::glmmTMBControl(parallel = 1)
custom_control_lmer_final <- lmerControl(
  optimizer = "bobyqa", optCtrl = list(maxfun = 200000), calc.derivs = TRUE)

standardize_id_within_exp <- function(ids) {
  numeric_id <- sub("^\\D*(\\d+).*$", "\\1", ids)
  if (any(!grepl("^\\d+$", numeric_id))) stop("Non-numeric participant id.")
  sprintf("Vp%04d", as.numeric(numeric_id))
}

# Reconstruct the same main and interaction contrast columns as the revised base model.
add_random_contrast_columns <- function(data, random_vars) {
  f <- reformulate(paste(random_vars, collapse = " * "))
  mm <- model.matrix(f, data, contrasts.arg = setNames(lapply(data[random_vars], function(x) contr.sum(nlevels(x))), random_vars))
  labels <- attr(terms(f), "term.labels")
  for (k in seq_along(labels)) {
    cols <- which(attr(mm, "assign") == k)
    for (j in seq_along(cols)) data[[sprintf("RE_%s_%d", gsub(":", "_x_", labels[k], fixed = TRUE), j)]] <- mm[, cols[j]]
  }
  data
}

# Verbatim logic of the pipeline's mixed_model_diagnostics().
mixed_model_diagnostics <- function(model) {
  is_glmm <- inherits(model, "glmmTMB")
  singular <- if (is_glmm) {
    vc <- VarCorr(model)$cond
    if (is.null(vc)) FALSE else {
      vars <- unlist(lapply(vc, function(x) diag(as.matrix(x))))
      any(!is.finite(vars)) || any(vars < 1e-6)
    }
  } else lme4::isSingular(model, tol = 1e-4)
  if (is_glmm) {
    optimizer_code <- model$fit$convergence %||% NA_integer_
    msgs <- model$fit$message %||% character(0)
    hess_avail <- !is.null(model$sdr$pdHess); hess_pd <- isTRUE(model$sdr$pdHess)
    fixed_se <- suppressWarnings(sqrt(diag(vcov(model)$cond)))
  } else {
    optimizer_code <- model@optinfo$conv$opt %||% NA_integer_
    msgs <- model@optinfo$conv$lme4$messages %||% character(0)
    hessian <- model@optinfo$derivs$Hessian
    hess_avail <- !is.null(hessian)
    hess_pd <- if (hess_avail) {
      ev <- tryCatch(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) NA_real_)
      all(is.finite(ev)) && min(ev) > 0
    } else FALSE
    fixed_se <- sqrt(diag(vcov(model)))
  }
  finite_se <- length(fixed_se) > 0 && all(is.finite(fixed_se))
  finite_ll <- is.finite(AIC(model)) && is.finite(as.numeric(logLik(model)))
  hard <- "failed to converge|unable to evaluate scaled gradient|degenerate.*Hessian|negative eigenvalue"
  messages_ok <- if (is_glmm) TRUE else !any(grepl(hard, msgs, ignore.case = TRUE))
  valid <- isTRUE(optimizer_code == 0) && messages_ok && !singular && finite_se &&
    finite_ll && hess_avail && hess_pd
  list(valid = valid, optimizer_code = optimizer_code,
       convergence_messages = paste(msgs, collapse = " | "), singular = singular,
       hessian_positive_definite = hess_pd, finite_fixed_se = finite_se,
       finite_likelihood = finite_ll)
}

prepare_data <- function(experiment) {
  path <- file.path(project_root, "data", paste0("02_Pipeline_Output_", experiment),
                    "Method_Regression", "Stimulus_Locked", "trials.csv")
  if (!file.exists(path)) stop("Trial file not found: ", path)
  md5 <- as.character(tools::md5sum(path))
  if (!identical(md5, unname(expected_md5[experiment])))
    stop(experiment, " trials.csv MD5 ", md5, " != REPRODUCIBILITY value ", expected_md5[experiment])
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("participant_id", "index", "stim", "emotion", "Offers_Other", "reaction", "RT")
  if (length(setdiff(req, names(df)))) stop("Missing columns: ", paste(setdiff(req, names(df)), collapse = ", "))
  if (any(!grepl("^[LR]_(neu|aff|dis|dom|enj)(Fema|Male)\\d+$", df$stim)))
    stop(experiment, ": unexpected stim values.")
  df$participant_id_internal <- standardize_id_within_exp(df$participant_id)
  df$actor_id <- sub("^[LR]_(neu|aff|dis|dom|enj)", "", df$stim, perl = TRUE)
  n0 <- nrow(df)
  df <- df[df$Offers_Other %in% c(5, 6, 8, 9), ]
  df <- df[df$reaction != 0, ]
  df <- df[df$RT >= rt_lower_ms & df$RT <= rt_upper_ms, ]
  df$offer_type <- factor(ifelse(df$Offers_Other %in% c(5, 6), "fair", "unfair"), levels = c("fair", "unfair"))
  df$emotion <- factor(df$emotion, levels = c("neu", "aff", "dis", "dom", "enj"))
  df$reaction <- factor(ifelse(df$reaction == 1, "accept", "reject"), levels = c("accept", "reject"))
  df$reject_binary <- ifelse(df$reaction == "reject", 1, 0)
  df$logRT <- log(df$RT)
  df$participant_id_internal <- factor(df$participant_id_internal)
  df$actor_id <- factor(df$actor_id)
  if (nlevels(df$actor_id) != 60L) stop(experiment, ": not 60 actors.")
  if (anyDuplicated(df[c("participant_id_internal", "index")])) stop("Duplicate participant x index.")
  df <- add_random_contrast_columns(df, c("emotion", "offer_type"))
  message(sprintf("[%s] %d -> %d trials, %d participants, %d actors", experiment, n0, nrow(df),
                  nlevels(df$participant_id_internal), nlevels(df$actor_id)))
  df
}

fit_glmm <- function(formula, data) {
  try_one <- function(ctrl, label) {
    m <- tryCatch(glmmTMB::glmmTMB(formula, data = data, family = binomial("logit"),
                                   control = ctrl, REML = FALSE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    attr(m, "optimizer_label") <- label
    attr(m, "diag") <- mixed_model_diagnostics(m)
    m
  }
  m <- try_one(custom_control_glmmTMB, "BFGS")
  attempts <- list(m)
  if (is.null(m) || !isTRUE(attr(m, "diag")$valid)) {
    m2 <- try_one(default_control_glmmTMB, "nlminb"); attempts <- c(attempts, list(m2))
    if (!is.null(m2) && (is.null(m) || isTRUE(attr(m2, "diag")$valid))) m <- m2
  }
  list(model = m, attempts = Filter(Negate(is.null), attempts))
}

fit_lmm <- function(formula, data) {
  m <- tryCatch(lmerTest::lmer(formula, data = data, REML = TRUE,
                               control = custom_control_lmer_final), error = function(e) NULL)
  if (!is.null(m)) { attr(m, "optimizer_label") <- "bobyqa"; attr(m, "diag") <- mixed_model_diagnostics(m) }
  list(model = m, attempts = Filter(Negate(is.null), list(m)))
}

diag_row <- function(m, experiment, outcome, label, n) {
  d <- if (is.null(m)) list(valid = FALSE, optimizer_code = NA, convergence_messages = "fit error",
                             singular = NA, hessian_positive_definite = NA, finite_fixed_se = NA,
                             finite_likelihood = NA) else attr(m, "diag")
  vc_actor <- NA_real_; vc_actor_slopes <- NA_character_; vc_part_int <- NA_real_
  if (!is.null(m)) {
    vc <- if (inherits(m, "glmmTMB")) VarCorr(m)$cond else VarCorr(m)
    nm <- names(vc)
    a_idx <- which(nm == "actor_id")
    if (length(a_idx)) {
      a_vars <- unlist(lapply(vc[a_idx], function(x) diag(as.matrix(x))))
      vc_actor <- a_vars[grepl("Intercept", names(a_vars))][1]
      slope_vars <- a_vars[!grepl("Intercept", names(a_vars))]
      if (length(slope_vars)) vc_actor_slopes <- paste(sprintf("%s=%.4g", names(slope_vars), slope_vars), collapse = "; ")
    }
    p_idx <- which(nm == "participant_id_internal")[1]
    if (!is.na(p_idx)) vc_part_int <- diag(as.matrix(vc[[p_idx]]))[1]
  }
  data.frame(experiment, outcome, model = label, n_obs = n,
             optimizer = if (is.null(m)) NA else attr(m, "optimizer_label"),
             valid = d$valid, optimizer_code = d$optimizer_code, messages = d$convergence_messages,
             singular = d$singular, hessian_pd = d$hessian_positive_definite,
             finite_fixed_se = d$finite_fixed_se, finite_likelihood = d$finite_likelihood,
             logLik = if (is.null(m)) NA else as.numeric(logLik(m)),
             AIC = if (is.null(m)) NA else AIC(m),
             participant_intercept_var = vc_part_int, actor_intercept_var = vc_actor,
             actor_slope_vars = vc_actor_slopes, stringsAsFactors = FALSE)
}

omnibus_tbl <- function(m, experiment, outcome, label) {
  if (inherits(m, "glmmTMB")) {
    a <- as.data.frame(car::Anova(m, type = 3, test.statistic = "Chisq"))
    out <- data.frame(experiment, outcome, model = label, term = rownames(a),
                      statistic = a$Chisq, df1 = a$Df, df2 = NA, p = a$`Pr(>Chisq)`, test = "Wald chi-square")
  } else {
    a <- as.data.frame(anova(m, type = 3, ddf = "Satterthwaite"))
    out <- data.frame(experiment, outcome, model = label, term = rownames(a),
                      statistic = a$`F value`, df1 = a$NumDF, df2 = a$DenDF, p = a$`Pr(>F)`, test = "Satterthwaite F")
  }
  rownames(out) <- NULL; out
}

fixed_tbl <- function(m, experiment, outcome, label) {
  cf <- if (inherits(m, "glmmTMB")) summary(m)$coefficients$cond else summary(m)$coefficients
  data.frame(experiment, outcome, model = label, term = rownames(cf),
             estimate = cf[, "Estimate"], SE = cf[, "Std. Error"], row.names = NULL)
}

contrast_tbl <- function(m, experiment, outcome, label, by_offer) {
  lmer_opts <- if (inherits(m, "glmmTMB")) list() else list(lmer.df = "satterthwaite")
  grid <- if (by_offer) emmeans(m, ~ emotion | offer_type, type = "link", lmerTest.limit = 1e6, lmer.df = "satterthwaite")
          else emmeans(m, ~ emotion, type = "link", lmerTest.limit = 1e6, lmer.df = "satterthwaite")
  raw <- as.data.frame(summary(pairs(grid, adjust = "none"), infer = c(TRUE, TRUE)))
  fdr <- as.data.frame(summary(pairs(grid, adjust = "fdr"), infer = c(TRUE, TRUE)))
  raw$p_fdr <- fdr$p.value; raw$CI_low_fdr <- fdr$lower.CL; raw$CI_high_fdr <- fdr$upper.CL
  cbind(experiment, outcome, model = label, family = if (by_offer) "emotion_within_offer" else "emotion_marginal", raw)
}

all_diag <- list(); all_omni <- list(); all_fixed <- list(); all_con <- list(); all_frame <- list()

for (experiment in "E1") {
  data <- prepare_data(experiment)
  specs <- list(
    choice = list(rds = file.path(reference_root, "Behavior", paste0(experiment, "_TwoStage_Bates"),
                                  "GLMM_Rejection/final_model_GLMM_rejection.rds"), fitter = fit_glmm),
    rt_main = list(rds = file.path(reference_root, "Behavior", paste0(experiment, "_TwoStage_Bates"),
                                   "LMM_RT_main/final_model_LMM_RT_main.rds"), fitter = fit_lmm)
  )
  for (outcome in "choice") {
    spec <- specs[[outcome]]
    rds_files <- list.files(dirname(spec$rds), pattern = "\\.rds$", full.names = TRUE)
    if (!file.exists(spec$rds)) {
      if (length(rds_files) == 1L) spec$rds <- rds_files else stop("Final model rds not found for ", experiment, " ", outcome,
                                                                    "; candidates: ", paste(rds_files, collapse = ", "))
    }
    saved <- readRDS(spec$rds)
    base_formula <- formula(saved)
    message(sprintf("[%s %s] canonical formula: %s", experiment, outcome, paste(deparse(base_formula), collapse = " ")))

    # Frame audit: the saved model frame must match the freshly prepared data row for row.
    sf <- if (inherits(saved, "glmmTMB")) saved$frame else model.frame(saved)
    resp <- if (outcome == "choice") "reject_binary" else "logRT"
    all_frame[[paste(experiment, outcome)]] <- data.frame(
      experiment, outcome, saved_rows = nrow(sf), current_rows = nrow(data),
      response_mismatch = if (nrow(sf) == nrow(data)) sum(abs(sf[[resp]] - data[[resp]]) > 1e-9) else NA,
      emotion_mismatch = if (nrow(sf) == nrow(data)) sum(as.character(sf$emotion) != as.character(data$emotion)) else NA,
      participant_mismatch = if (nrow(sf) == nrow(data)) sum(as.character(sf$participant_id_internal) != as.character(data$participant_id_internal)) else NA,
      re_col_max_abs_diff = if (nrow(sf) == nrow(data)) {
        re_cols <- grep("^RE_", names(sf), value = TRUE)
        max(abs(as.matrix(sf[, re_cols, drop = FALSE]) - as.matrix(data[, re_cols, drop = FALSE])))
      } else NA)

    audit <- all_frame[[paste(experiment, outcome)]]
    stopifnot(audit$saved_rows == audit$current_rows,
              audit$response_mismatch == 0, audit$emotion_mismatch == 0,
              audit$participant_mismatch == 0, audit$re_col_max_abs_diff == 0)
    actor_int_formula <- update(base_formula, . ~ . + (1 | actor_id))
    actor_slope_formula <- update(base_formula, . ~ . + (1 | actor_id) + (0 + RE_emotion_1 | actor_id) +
                                    (0 + RE_emotion_2 | actor_id) + (0 + RE_emotion_3 | actor_id) + (0 + RE_emotion_4 | actor_id))
    models <- list(canonical_refit = base_formula, actor_intercept = actor_int_formula,
                   actor_intercept_emotion_slopes = actor_slope_formula)
    for (label in names(models)) {
      message(sprintf("  fitting %s ...", label))
      res <- spec$fitter(models[[label]], data)
      m <- res$model
      for (att in res$attempts) all_diag[[length(all_diag) + 1]] <- diag_row(att, experiment, outcome,
                                                                        paste0(label, " [attempt ", attr(att, "optimizer_label"), "]"), nrow(data))
      all_diag[[length(all_diag) + 1]] <- diag_row(m, experiment, outcome, label, nrow(data))
      if (!is.null(m)) saveRDS(m, file.path(output_root, sprintf("model_%s_%s_%s.rds", experiment, outcome, label)))
      if (!is.null(m) && isTRUE(attr(m, "diag")$valid)) {
        all_omni[[length(all_omni) + 1]] <- omnibus_tbl(m, experiment, outcome, label)
        all_fixed[[length(all_fixed) + 1]] <- fixed_tbl(m, experiment, outcome, label)
        all_con[[length(all_con) + 1]] <- contrast_tbl(m, experiment, outcome, label, FALSE)
        all_con[[length(all_con) + 1]] <- contrast_tbl(m, experiment, outcome, label, TRUE)
      } else message("    -> INVALID under the acceptance gate; estimates not extracted.")
    }
    # Saved canonical fit's own fixed effects, for comparison with the refit.
    all_fixed[[length(all_fixed) + 1]] <- fixed_tbl(saved, experiment, outcome, "canonical_saved_rds")
    rm(saved); invisible(gc())
  }
  rm(data); invisible(gc())
}

bind <- function(l) {
  if (!length(l)) return(data.frame())
  columns <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[, columns, drop = FALSE]
  }))
}
write.csv(bind(all_diag), file.path(output_root, "diagnostics.csv"), row.names = FALSE, na = "")
write.csv(bind(all_omni), file.path(output_root, "omnibus_comparison.csv"), row.names = FALSE, na = "")
write.csv(bind(all_fixed), file.path(output_root, "fixed_effects_comparison.csv"), row.names = FALSE, na = "")
write.csv(bind(all_con), file.path(output_root, "emotion_contrasts_comparison.csv"), row.names = FALSE, na = "")
write.csv(bind(all_frame), file.path(output_root, "saved_frame_audit.csv"), row.names = FALSE, na = "")

writeLines(c(
  "E1 choice actor-sensitivity update after the random-structure revision",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Inputs (READ ONLY, live analysis tree):",
  paste0("- ", project_root, "/data/02_Pipeline_Output_E1/Method_Regression/Stimulus_Locked/trials.csv (MD5 verified against REPRODUCIBILITY.txt of 2026-08-22)"),
  paste0("- ", reference_root, "/Behavior/E1_TwoStage_Bates/GLMM_Rejection/final_model_GLMM_rejection.rds"),
  "Filters: Offers_Other in {5,6,8,9}; reaction != 0; RT in [300,3000] ms (identical to pipeline).",
  "Actor derivation: stim with ^[LR]_(neu|aff|dis|dom|enj) removed (60 actors).",
  "RE_* columns: sum-contrast columns of emotion, offer_type and their interaction; checked against the revised saved frame.",
  "Models per outcome: canonical_refit (exact saved formula), actor_intercept (+ (1|actor_id)),",
  "  actor_intercept_emotion_slopes (+ (1|actor_id) + four diagonal (0+RE_emotion_k|actor_id) terms).",
  "Acceptance gate (pipeline mixed_model_diagnostics): optimizer code 0, positive-definite Hessian,",
  "  finite fixed SE, finite likelihood, no variance < 1e-6 (glmmTMB) / isSingular tol 1e-4 (lme4).",
  "GLMM optimizer order: BFGS then nlminb (pipeline stabilize_final_glmm); LMM: bobyqa, calc.derivs = TRUE.",
  "Contrasts: emmeans pairs on link scale; raw and FDR-adjusted p and CI both reported.",
  "Original input and result files were not modified; new outputs use a versioned directory.",
  paste0("Revised base-model run: ", reference_root),
  paste0("R: ", R.version.string),
  paste0("glmmTMB ", packageVersion("glmmTMB"), "; lme4 ", packageVersion("lme4"), "; lmerTest ", packageVersion("lmerTest"),
         "; emmeans ", packageVersion("emmeans"), "; car ", packageVersion("car"))
), file.path(output_root, "REPRODUCIBILITY.txt"))
message("Done: ", output_root)
