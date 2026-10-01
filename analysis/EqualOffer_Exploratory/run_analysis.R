#!/usr/bin/env Rscript
# Usage from project root:
# Rscript analysis/EqualOffer_Exploratory/run_analysis.R --stage=audit --out=results/EqualOffer_audit_v1
# Rscript analysis/EqualOffer_Exploratory/run_analysis.R --stage=all --out=results/EqualOffer_full_v1
# All inputs read-only. Output directory must not already exist.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- args[startsWith(args, paste0("--", key, "="))]
  if (length(hit) > 1L) stop("Repeated argument: ", key)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit) else default
}
if (any(!grepl("^--(stage|out|root)=", args))) stop("Allowed arguments: --stage=, --out=, --root=")
stage <- get_arg("stage", "audit")
if (!stage %in% c("audit", "all")) stop("stage must be audit or all")
find_root <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(path, "UG_ERP_Project.Rproj"))) {
    parent <- dirname(path)
    if (identical(parent, path)) stop("Project root not found; use --root=PATH")
    path <- parent
  }
  path
}
script_flag <- grep("^--file=", commandArgs(), value = TRUE)
start <- if (length(script_flag)) dirname(sub("^--file=", "", script_flag[1])) else getwd()
root <- find_root(get_arg("root", start))
code_dir <- file.path(root, "analysis", "EqualOffer_Exploratory")
source(file.path(code_dir, "config.R"))
# Use the existing project library without installing or updating anything.
rlib <- file.path(root, "renv", "library", paste0("R-", R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]), R.version$platform)
if (dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))
required <- c("pacman", "readr", "dplyr", "tidyr", "stringr", "jsonlite", "digest", "glmmTMB",
              "lme4", "lmerTest", "emmeans", "car", "ggplot2", "R.utils")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                          ". Restore the project environment first; see README.md.")
suppressPackageStartupMessages(pacman::p_load(char = required, install = FALSE))
options(contrasts = c("contr.sum", "contr.poly"), digits = 17)
set.seed(seed)
out_arg <- get_arg("out", paste0("results/EqualOffer_", stage, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
out <- if (grepl("^([A-Za-z]:|/)", out_arg)) out_arg else file.path(root, out_arg)
# Restrict outputs to a new child of the project's results directory.
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
out_parent <- normalizePath(dirname(out), winslash = "/", mustWork = TRUE)
results_root <- normalizePath(file.path(root, "results"), winslash = "/", mustWork = TRUE)
if (!(identical(tolower(out_parent), tolower(results_root)) ||
      startsWith(tolower(out_parent), paste0(tolower(results_root), "/")))) stop("Output must be inside project results/")
if (dir.exists(out) || file.exists(out)) stop("Refusing to overwrite output: ", out)
dir.create(out)
out <- normalizePath(out, winslash = "/")
write_csv <- function(x, name) readr::write_csv(x, file.path(out, name), na = "NA")
writeLines(c("RUNNING", paste("stage:", stage)), file.path(out, "STATUS.txt"))
file.copy(file.path(code_dir, c("run_analysis.R", "config.R")), out)
writeLines(c(paste("Command:", paste(commandArgs(), collapse = " ")), capture.output(sessionInfo())),
           file.path(out, "SESSION_INFO.txt"))
manifest <- list()
record_file <- function(path, role) {
  if (!file.exists(path)) stop("Required file missing: ", path)
  manifest[[length(manifest) + 1L]] <<- data.frame(role = role,
    path = normalizePath(path, winslash = "/"), bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path)), sha256 = digest::digest(file = path, algo = "sha256"))
}
for (f in c("config.R", "run_analysis.R")) record_file(file.path(code_dir, f), "analysis_code")
record_file(file.path(root, "Sta_Behaviour_E1_E2_Integrative.Rmd"), "original_analysis_code")
if (file.exists(file.path(root, "renv.lock"))) record_file(file.path(root, "renv.lock"), "environment_lock")

# 1. Canonical trial audit; keep original IDs/labels/row order and add explicit columns.
data_by_exp <- list(); checks <- list(); participant_tables <- list(); cell_tables <- list()
exclusion_tables <- list(); concentration_tables <- list()
for (ex in experiments) {
  message("Auditing ", ex)
  path <- file.path(root, "data", paste0("02_Pipeline_Output_", ex),
                    "Method_Regression", "Stimulus_Locked", "trials.csv")
  record_file(path, paste0(ex, "_canonical_trials"))
  if (unname(tools::md5sum(path)) != expected_input_md5[[ex]])
    stop(ex, ": canonical input has changed from the manuscript run. Review provenance before proceeding.")
  d <- readr::read_csv(path, show_col_types = FALSE)
  needed <- c("participant_id", "index", "stim", "Offers_Other", "Offers_You", "emotion", "reaction", "RT")
  if (length(setdiff(needed, names(d)))) stop("Missing columns: ", paste(setdiff(needed, names(d)), collapse = ", "))
  if (anyNA(d[needed])) stop(ex, ": missing values in required fields; review explicitly.")
  if (anyDuplicated(d[c("participant_id", "index")])) stop(ex, ": duplicate participant/trial keys")
  if (!all(d$emotion %in% emotion_levels) || !all(d$reaction %in% 0:2) ||
      !all(d$Offers_Other %in% 5:9) || !all(d$Offers_Other + d$Offers_You == 10) ||
      !all(is.finite(d$RT))) stop(ex, ": unexpected design or response values")
  d$source_row <- seq_len(nrow(d))
  # Reproduce, and explicitly document, the original loader's internal ID mapping.
  d$participant_id_internal <- sprintf("Vp%04d", as.numeric(stringr::str_extract(d$participant_id, "\\d+")))
  id_map <- unique(d[c("participant_id", "participant_id_internal")])
  if (anyNA(id_map) || anyDuplicated(id_map$participant_id_internal)) stop(ex, ": ID collision")
  write_csv(id_map, paste0(ex, "_id_mapping.csv"))
  d$exclusion_reason <- ifelse(!d$Offers_Other %in% c(5, 6, 8, 9), "filler_7",
    ifelse(d$reaction == 0, "timeout",
      ifelse(d$RT < rt_bounds_ms[1] | d$RT > rt_bounds_ms[2], "outside_RT_window", "included")))
  exclusion_tables[[ex]] <- d %>% count(Offers_Other, emotion, exclusion_reason, name = "trials") %>% mutate(experiment = ex)
  clean <- d[d$exclusion_reason == "included", ]
  clean$emotion <- factor(clean$emotion, levels = emotion_levels)
  clean$offer_type <- factor(ifelse(clean$Offers_Other %in% c(5, 6), "fair", "unfair"), levels = c("fair", "unfair"))
  clean$reject_binary <- as.integer(clean$reaction == 2)
  # Cross-check every included observation with the saved manuscript choice frame.
  ref_base <- file.path(root, reference_runs[[ex]])
  ref_path <- file.path(ref_base, "GLMM_Rejection", "final_model_GLMM_rejection.rds")
  record_file(ref_path, paste0(ex, "_manuscript_choice_model"))
  ref <- readRDS(ref_path); frame <- model.frame(ref)
  audit <- c(rows = nrow(frame) == nrow(clean),
    response = identical(as.numeric(frame$reject_binary), as.numeric(clean$reject_binary)),
    emotion = identical(as.character(frame$emotion), as.character(clean$emotion)),
    offer = identical(as.character(frame$offer_type), as.character(clean$offer_type)),
    participants = identical(as.character(frame$participant_id_internal), clean$participant_id_internal))
  checks[[ex]] <- data.frame(experiment = ex, check = names(audit), passed = audit)
  if (!all(audit)) stop(ex, ": saved manuscript model frame does not match canonical data")
  # Retain every participant, including all-zero rejecters; missing cells explicit.
  fair <- clean[clean$Offers_Other %in% c(5, 6), ]
  fair$allocation <- factor(fair$Offers_Other, levels = c(5, 6), labels = c("5:5", "6:4"))
  fair$participant_id_internal <- factor(fair$participant_id_internal)
  if (!all(grepl("^[LR]_(neu|aff|dis|dom|enj)(Fema|Male)[0-9]+$", fair$stim))) stop("Unrecognized actor labels")
  fair$actor_id <- factor(sub("^[LR]_(neu|aff|dis|dom|enj)", "", fair$stim))
  # Sum-coded numeric slopes avoid accidental correlated factor blocks with ||.
  mm <- model.matrix(~ emotion * allocation, fair)[, -1, drop = FALSE]
  slope_names <- sprintf("RE%02d", seq_len(ncol(mm)))
  fair[slope_names] <- mm
  write_csv(data.frame(column = slope_names, design_term = colnames(mm)), paste0(ex, "_random_design.csv"))
  data_by_exp[[ex]] <- fair
  write_csv(fair[c("source_row", "participant_id", "participant_id_internal", "index", "stim", "actor_id",
                   "Offers_Other", "Offers_You", "allocation", "emotion", "reaction", "reject_binary", "RT")],
            paste0(ex, "_analysis_trials.csv"))
  p <- fair %>% group_by(participant_id, emotion, allocation) %>%
    summarise(valid_trials = n(), rejections = sum(reject_binary), .groups = "drop") %>%
    tidyr::complete(participant_id = unique(d$participant_id), emotion, allocation,
                    fill = list(valid_trials = 0L, rejections = 0L)) %>%
    mutate(experiment = ex, rejection_rate = ifelse(valid_trials > 0, rejections / valid_trials, NA_real_),
           missing_cell = valid_trials == 0)
  participant_tables[[ex]] <- p
  cell_tables[[ex]] <- p %>% group_by(experiment, allocation, emotion) %>%
    summarise(participants_with_trials = sum(valid_trials > 0), participants_with_rejections = sum(rejections > 0),
      trial_weighted_rate = sum(rejections) / sum(valid_trials),
      participant_weighted_rate = mean(rejection_rate, na.rm = TRUE),
      largest_participant_share_of_rejections = ifelse(sum(rejections) > 0, max(rejections) / sum(rejections), NA_real_),
      valid_trials = sum(valid_trials), rejections = sum(rejections),
      .groups = "drop")
  concentration_tables[[ex]] <- p %>% group_by(experiment, participant_id, allocation) %>%
    summarise(valid_trials = sum(valid_trials), rejections = sum(rejections), .groups = "drop")
}
write_csv(bind_rows(checks), "manuscript_frame_checks.csv")
write_csv(bind_rows(exclusion_tables), "exclusions.csv")
write_csv(bind_rows(participant_tables), "participant_cells.csv")
write_csv(bind_rows(cell_tables), "descriptive_cells.csv")
write_csv(bind_rows(concentration_tables), "participant_rejection_totals.csv")

# 2. Extract existing choice / all-trial RT evidence without refitting originals.
existing_omnibus <- list(); existing_contrasts <- list(); existing_models <- list()
pair_weights <- function(grid, pairs, offer_var, offer_levels) {
  vectors <- list()
  for (nm in names(pairs)) {
    pair <- pairs[[nm]]
    emotion_weight <- as.numeric(grid$emotion == pair[1]) - as.numeric(grid$emotion == pair[2])
    for (lev in offer_levels) vectors[[paste(nm, lev, sep = " | ")]] <- emotion_weight * (grid[[offer_var]] == lev)
    vectors[[paste(nm, paste0(offer_levels[2], " minus ", offer_levels[1]), sep = " | ")]] <-
      emotion_weight * (as.numeric(grid[[offer_var]] == offer_levels[2]) - as.numeric(grid[[offer_var]] == offer_levels[1]))
  }
  vectors
}
for (ex in experiments) for (outcome in c("GLMM_rejection", "LMM_RT_main")) {
  folder <- if (outcome == "GLMM_rejection") "GLMM_Rejection" else outcome
  path <- file.path(root, reference_runs[[ex]], folder, paste0("final_model_", outcome, ".rds"))
  record_file(path, paste(ex, outcome, "reference"))
  model <- readRDS(path)
  tab <- if (outcome == "GLMM_rejection") as.data.frame(car::Anova(model, type = 3)) else
    as.data.frame(anova(lmerTest::as_lmerModLmerTest(model), type = 3, ddf = "Satterthwaite"))
  tab$term <- rownames(tab); rownames(tab) <- NULL
  tab$experiment <- ex; tab$outcome <- outcome
  existing_omnibus[[paste(ex, outcome)]] <- tab
  eg <- emmeans::emmeans(model, ~ emotion * offer_type, lmer.df = "asymptotic")
  grid <- as.data.frame(eg)
  cs <- as.data.frame(summary(emmeans::contrast(eg,
    method = pair_weights(grid, c(focal_pairs, secondary_pairs), "offer_type", c("fair", "unfair"))),
    infer = c(TRUE, TRUE), adjust = "none"))
  cs$experiment <- ex; cs$outcome <- outcome
  cs$scale <- ifelse(outcome == "GLMM_rejection", "log_odds", "log_RT")
  cs$p_holm_within_outcome_experiment <- p.adjust(cs$p.value, "holm")
  existing_contrasts[[paste(ex, outcome)]] <- cs
  existing_models[[paste(ex, outcome)]] <- data.frame(experiment = ex, outcome = outcome,
    formula = paste(deparse(formula(model)), collapse = " "), observations = nobs(model),
    source = path, md5 = unname(tools::md5sum(path)))
}
write_csv(bind_rows(existing_omnibus), "existing_omnibus.csv")
write_csv(bind_rows(existing_contrasts), "existing_targeted_contrasts.csv")
write_csv(bind_rows(existing_models), "existing_model_sources.csv")
write_csv(bind_rows(manifest), "input_manifest.csv")

# 3. New fair-only GLMM: diagnostic-based simplification, never p-value selection.
# Start with all nine diagonal slopes. Try nlminb, then BFGS if numerically invalid.
# Remove only slope SDs < .001, one at a time, with the intercept always retained.
# If no valid model results, stop inference; do not silently switch to intercept-only.
diagnose <- function(m) {
  if (inherits(m, "error")) return(list(valid = FALSE, reason = conditionMessage(m)))
  cf <- summary(m)$coefficients$cond
  numerically_valid <- isTRUE(m$sdr$pdHess) && m$fit$convergence == 0 &&
    all(is.finite(cf)) && is.finite(as.numeric(logLik(m)))
  extreme <- any(abs(cf[, "Estimate"]) > 15 | cf[, "Std. Error"] > 10)
  sd <- vapply(glmmTMB::VarCorr(m)$cond, function(v) unname(attr(v, "stddev")[1]), numeric(1))
  list(valid = numerically_valid && !extreme, hessian_positive = isTRUE(m$sdr$pdHess),
       optimizer_code = m$fit$convergence, extreme_fixed = extreme,
       random_boundary = any(sd < random_sd_boundary), minimum_random_sd = min(sd),
       reason = paste(m$fit$message, collapse = "; "))
}
make_formula <- function(slopes, actor = FALSE) {
  terms <- c("emotion * allocation", "(1 | participant_id_internal)",
             paste0("(0 + ", slopes, " | participant_id_internal)"),
             if (actor) "(1 | actor_id)")
  as.formula(paste("reject_binary ~", paste(terms, collapse = " + ")))
}
fit_one <- function(d, formula, optimizer = "nlminb") {
  control <- if (optimizer == "nlminb")
    glmmTMB::glmmTMBControl(optCtrl = list(iter.max = max_evaluations, eval.max = max_evaluations)) else
    glmmTMB::glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
                          optCtrl = list(maxit = max_evaluations))
  tryCatch(R.utils::withTimeout(glmmTMB::glmmTMB(formula, data = d, family = binomial(), control = control),
                               timeout = fit_timeout_seconds, onTimeout = "error"), error = identity)
}
select_model <- function(d, ex) {
  slopes <- grep("^RE[0-9]+$", names(d), value = TRUE)
  log <- list()
  for (iteration in seq_len(length(slopes) + 1L)) {
    fo <- make_formula(slopes)
    message(ex, " model iteration ", iteration, "; slopes=", paste(slopes, collapse = ","))
    model <- fit_one(d, fo)
    diag <- diagnose(model)
    log[[length(log) + 1L]] <- c(list(iteration = iteration, optimizer = "nlminb", formula = deparse(fo)), diag)
    if (!diag$valid) {
      alternate <- fit_one(d, fo, "BFGS")
      alt_diag <- diagnose(alternate)
      log[[length(log) + 1L]] <- c(list(iteration = iteration, optimizer = "BFGS", formula = deparse(fo)), alt_diag)
      if (alt_diag$valid || inherits(model, "error")) { model <- alternate; diag <- alt_diag }
    }
    jsonlite::write_json(log, file.path(out, paste0(ex, "_selection_log.json")), auto_unbox = TRUE, pretty = TRUE)
    if (inherits(model, "error")) return(list(model = NULL, log = log))
    vc <- glmmTMB::VarCorr(model)$cond
    sd <- vapply(vc, function(v) unname(attr(v, "stddev")[1]), numeric(1))
    term <- vapply(vc, function(v) rownames(v)[1], character(1))
    boundary <- which(term %in% slopes & is.finite(sd) & sd < random_sd_boundary)
    if (length(boundary)) {
      drop <- term[boundary[which.min(sd[boundary])]]
      log[[length(log) + 1L]] <- list(action = "remove_boundary_slope", term = drop,
        explanation = "Numerical simplification; fixed-effect tests were not inspected.")
      slopes <- setdiff(slopes, drop)
      next
    }
    if (!diag$valid || any(sd < random_sd_boundary)) return(list(model = NULL, log = log))
    return(list(model = model, slopes = slopes, log = log))
  }
  list(model = NULL, log = log)
}
extract_contrasts <- function(m, ex, label) {
  eg <- emmeans::emmeans(m, ~ emotion * allocation)
  pg <- emmeans::regrid(eg, transform = "response")
  rows <- list()
  for (family in c("focal", "secondary")) {
    pairs <- if (family == "focal") focal_pairs else secondary_pairs
    weights <- pair_weights(as.data.frame(eg), pairs, "allocation", c("5:5", "6:4"))
    for (scale in c("log_odds", "probability_difference")) {
      obj <- if (scale == "log_odds") eg else pg
      t <- as.data.frame(summary(emmeans::contrast(obj, method = weights), infer = c(TRUE, TRUE), adjust = "none"))
      t$experiment <- ex; t$model <- label; t$family <- family; t$scale <- scale
      rows[[paste(family, scale)]] <- t
    }
  }
  bind_rows(rows)
}
if (stage == "all") {
  all_contrasts <- list(); all_probabilities <- list(); diagnostic_rows <- list(); influence_rows <- list()
  for (ex in experiments) {
    d <- data_by_exp[[ex]]
    selected <- select_model(d, ex)
    if (is.null(selected$model)) stop(ex, ": no stable fair-only GLMM. Audit outputs retained; regularized-model review required.")
    model <- selected$model
    saveRDS(model, file.path(out, paste0(ex, "_fair_choice.rds")))
    diagnostic_rows[[paste(ex, "primary")]] <- data.frame(experiment = ex, model = "primary",
      valid = diagnose(model)$valid, random_boundary = diagnose(model)$random_boundary,
      minimum_random_sd = diagnose(model)$minimum_random_sd, reason = diagnose(model)$reason)
    jsonlite::write_json(selected$log, file.path(out, paste0(ex, "_selection_log.json")), auto_unbox = TRUE, pretty = TRUE)
    # Avoid a glmmTMB 1.1.14/R 4.3 print-method compatibility issue (%||%).
    vc <- glmmTMB::VarCorr(model)$cond
    variance_table <- data.frame(group = names(vc),
      term = vapply(vc, function(v) rownames(v)[1], character(1)),
      sd = vapply(vc, function(v) unname(attr(v, "stddev")[1]), numeric(1)))
    write_csv(variance_table, paste0(ex, "_random_effects.csv"))
    writeLines(c(paste("Formula:", paste(deparse(formula(model)), collapse = " ")),
      capture.output(summary(model)$coefficients$cond), capture.output(variance_table),
      capture.output(diagnose(model))), file.path(out, paste0(ex, "_model_summary.txt")))
    writeLines(c("Exploratory fair-only trial-level Bernoulli GLMM with logit link.",
      paste("Formula:", paste(deparse(formula(model)), collapse = " ")),
      "Laplace maximum likelihood; zero random-effect correlations; boundary-slope simplification.",
      "All zero-rejection participants retained. Original valid-response and RT rules retained.",
      "Fixed-effect selection was not performed. Wald intervals ignore random-structure selection uncertainty.",
      "Probabilities refer to random effects set to zero, not population-marginal rejection rates.",
      "No observation-level dispersion parameter: individual responses are Bernoulli.",
      "Input is the manuscript's preprocessing output, not all original behavioral logs."),
      file.path(out, paste0("Methods_paragraph_", ex, ".md")))
    all_contrasts[[paste(ex, "primary")]] <- extract_contrasts(model, ex, "primary")
    write_csv(bind_rows(all_contrasts), "new_contrasts_unadjusted_checkpoint.csv")
    prob <- as.data.frame(summary(emmeans::emmeans(model, ~ emotion * allocation), type = "response", infer = c(TRUE, FALSE)))
    prob$experiment <- ex; all_probabilities[[ex]] <- prob
    write_csv(bind_rows(all_probabilities), "new_probabilities.csv")
    # Parametric model checks: simulate new random effects, no parameter-uncertainty claim.
    sims <- simulate(model, nsim = 200, seed = seed + match(ex, experiments))
    # glmmTMB binomial simulations return cbind(successes, failures), even for Bernoulli.
    sims <- lapply(sims, function(y) {
      if (is.matrix(y)) {
        if (ncol(y) != 2 || any(rowSums(y) != 1)) stop("Unexpected binomial simulation size")
        y <- y[, 1]
      }
      if (length(y) != nrow(d) || any(!y %in% 0:1)) stop("Invalid Bernoulli simulation output")
      y
    })
    group <- interaction(d$allocation, d$emotion, drop = TRUE)
    observed <- tapply(d$reject_binary, group, sum)
    simulated <- vapply(sims, function(y) as.numeric(tapply(y, group, sum)), numeric(length(observed)))
    checks_sim <- data.frame(cell = names(observed), observed_rejections = as.numeric(observed),
      simulated_median = apply(simulated, 1, median),
      simulated_lower = apply(simulated, 1, quantile, probs = .025),
      simulated_upper = apply(simulated, 1, quantile, probs = .975))
    zero <- vapply(sims, function(y) sum(tapply(y, d$participant_id_internal, sum) == 0), numeric(1))
    checks_sim <- rbind(checks_sim, data.frame(cell = "participants_with_zero_fair_rejections",
      observed_rejections = sum(tapply(d$reject_binary, d$participant_id_internal, sum) == 0),
      simulated_median = median(zero), simulated_lower = unname(quantile(zero, .025)),
      simulated_upper = unname(quantile(zero, .975))))
    write_csv(checks_sim, paste0(ex, "_parametric_model_check.csv"))
    rm(sims, simulated)
    # Crossed actor intercept: preserve selected participant structure, no new selection.
    actor <- fit_one(d, make_formula(selected$slopes, actor = TRUE))
    ad <- diagnose(actor)
    if (!ad$valid) { actor <- fit_one(d, make_formula(selected$slopes, actor = TRUE), "BFGS"); ad <- diagnose(actor) }
    if (!inherits(actor, "error")) saveRDS(actor, file.path(out, paste0(ex, "_actor_sensitivity.rds")))
    ad$valid <- ad$valid && !isTRUE(ad$random_boundary)
    diagnostic_rows[[paste(ex, "actor")]] <- data.frame(experiment = ex, model = "actor_intercept", valid = ad$valid,
      random_boundary = ifelse(is.null(ad$random_boundary), NA, ad$random_boundary),
      minimum_random_sd = ifelse(is.null(ad$minimum_random_sd), NA, ad$minimum_random_sd), reason = ad$reason)
    if (ad$valid) all_contrasts[[paste(ex, "actor")]] <- extract_contrasts(actor, ex, "actor_intercept")
    # Influence: no deletion from the primary analysis; all participants assessed.
    # Keep the same formula to isolate deletion influence; invalid fits explicitly retained.
    for (id in levels(d$participant_id_internal)) {
      message(ex, " leave-one-participant-out: ", id)
      loo <- fit_one(d[d$participant_id_internal != id, ], formula(model))
      ld <- diagnose(loo)
      if (!ld$valid) { loo <- fit_one(d[d$participant_id_internal != id, ], formula(model), "BFGS"); ld <- diagnose(loo) }
      ld$valid <- ld$valid && !isTRUE(ld$random_boundary)
      if (ld$valid) {
        z <- extract_contrasts(loo, ex, "leave_one_out")
        z$omitted_participant_id_internal <- id; z$valid <- TRUE
        influence_rows[[paste(ex, id)]] <- z
      } else influence_rows[[paste(ex, id)]] <- data.frame(experiment = ex,
        omitted_participant_id_internal = id, valid = FALSE,
        reason = paste(ld$reason, "random_boundary:", ld$random_boundary))
      write_csv(bind_rows(influence_rows), "participant_influence.csv")
    }
  }
  ct <- bind_rows(all_contrasts) %>% group_by(model, family, scale) %>%
    mutate(p_holm = if (first(scale) == "log_odds") p.adjust(p.value, "holm", n = ifelse(first(family) == "focal", 12L, 6L)) else NA_real_) %>% ungroup()
  ct$odds_ratio_or_ratio_of_OR <- ifelse(ct$scale == "log_odds", exp(ct$estimate), NA_real_)
  ct$OR_lower <- ifelse(ct$scale == "log_odds", exp(ct$asymp.LCL), NA_real_)
  ct$OR_upper <- ifelse(ct$scale == "log_odds", exp(ct$asymp.UCL), NA_real_)
  write_csv(ct, "new_contrasts.csv")
  write_csv(bind_rows(all_probabilities), "new_probabilities.csv")
  write_csv(bind_rows(diagnostic_rows), "sensitivity_diagnostics.csv")
}
writeLines(c("Retrospective exploratory analysis; not preregistered.",
  paste("Stage:", stage), paste("Seed:", seed),
  "Canonical input MD5 verified against formal manuscript run; saved choice frames matched row by row.",
  "Original results are read-only; existing model evidence is extracted without refitting.",
  "Input manifest includes SHA256, MD5, file size and exact paths; code snapshots and sessionInfo included.",
  "New contrasts: focal Holm family of 12; secondary family of 6, across both experiments.",
  "All intervals are unadjusted 95% Wald intervals. Probability contrasts are companion estimates.",
  "Existing targeted contrasts are a NEW descriptive audit, with separate nine-test Holm families per experiment/outcome.",
  "No subjective WTR, AV, status threat or intention was measured by this analysis.",
  "No mechanism, equivalence, absence, or theory winner is inferred from a p-value."),
  file.path(out, "REPRODUCIBILITY.txt"))
writeLines(c("COMPLETE", paste("stage:", stage)), file.path(out, "STATUS.txt"))
message("Completed: ", out)
