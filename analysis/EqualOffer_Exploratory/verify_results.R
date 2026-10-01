#!/usr/bin/env Rscript
# Independent arithmetic/count checks plus deterministic refits of selected models.
# Rscript --vanilla analysis/EqualOffer_Exploratory/verify_results.R results/EqualOffer_full_20260928_v3
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 1) stop("Supply exactly one completed output directory")
run <- normalizePath(a[1], winslash = "/", mustWork = TRUE)
root <- run
while (!file.exists(file.path(root, "UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if (parent == root) stop("Project root not found"); root <- parent
}
lib <- file.path(root, "renv", "library", paste0("R-", R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]), R.version$platform)
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(glmmTMB); library(jsonlite)})
source(file.path(run, "config.R"))
options(contrasts = c("contr.sum", "contr.poly"), digits = 17)
if (!identical(readLines(file.path(run, "STATUS.txt"))[1], "COMPLETE")) stop("Run incomplete")
dest <- file.path(run, "verification_checks.csv")
if (file.exists(dest)) stop("Verification already exists; preserve it and choose a new run for a new validation")
read <- function(f) read.csv(file.path(run, f), check.names = FALSE, stringsAsFactors = FALSE)
tests <- list()
check <- function(name, passed, detail = "") {
  tests[[length(tests) + 1L]] <<- data.frame(check = name, passed = isTRUE(passed), detail = detail)
  write.csv(do.call(rbind, tests), dest, row.names = FALSE)
  if (!isTRUE(passed)) stop("Verification failed: ", name, " ", detail)
}
ct <- read("new_contrasts.csv"); cells <- read("descriptive_cells.csv")
pc <- read("participant_cells.csv")
manifest <- read("input_manifest.csv")
for (i in seq_len(nrow(manifest))) check(paste0("source_unchanged_", i),
  identical(unname(tools::md5sum(manifest$path[i])), manifest$md5[i]), manifest$path[i])
max_diff <- function(x, y) max(abs(as.numeric(x) - as.numeric(y)))
refit_comparison <- list()
for (ex in experiments) {
  original <- read.csv(file.path(root, "data", paste0("02_Pipeline_Output_", ex),
    "Method_Regression", "Stimulus_Locked", "trials.csv"), stringsAsFactors = FALSE)
  keep <- original$Offers_Other %in% c(5, 6) & original$reaction %in% c(1, 2) &
          original$RT >= rt_bounds_ms[1] & original$RT <= rt_bounds_ms[2]
  d <- original[keep, ]
  exported <- read(paste0(ex, "_analysis_trials.csv"))
  check(paste0(ex, "_independent_trial_keys"), identical(paste(d$participant_id, d$index),
    paste(exported$participant_id, exported$index)))
  check(paste0(ex, "_outcome_recoding"), identical(as.integer(d$reaction == 2), exported$reject_binary))
  check(paste0(ex, "_participant_cells_count"), sum(pc$experiment == ex) == length(unique(original$participant_id)) * 10)
  for (o in c(5, 6)) for (em in emotion_levels) {
    sub <- d[d$Offers_Other == o & d$emotion == em, ]
    alloc <- ifelse(o == 5, "5:5", "6:4")
    row <- cells[cells$experiment == ex & cells$allocation == alloc & cells$emotion == em, ]
    rej_by_id <- tapply(sub$reaction == 2, sub$participant_id, sum)
    rates <- tapply(as.numeric(sub$reaction == 2), sub$participant_id, mean)
    check(paste(ex, o, em, "counts", sep = "_"), nrow(row) == 1 &&
      row$valid_trials == nrow(sub) && row$rejections == sum(sub$reaction == 2) &&
      row$participants_with_trials == length(unique(sub$participant_id)) &&
      row$participants_with_rejections == sum(rej_by_id > 0))
    check(paste(ex, o, em, "rates_and_concentration", sep = "_"),
      max_diff(c(row$trial_weighted_rate, row$participant_weighted_rate, row$largest_participant_share_of_rejections),
               c(mean(sub$reaction == 2), mean(rates), max(rej_by_id) / sum(rej_by_id))) < 1e-12)
    ps <- pc[pc$experiment == ex & pc$allocation == alloc & pc$emotion == em, ]
    for (id in unique(original$participant_id)) {
      z <- sub[sub$participant_id == id, ]; pp <- ps[ps$participant_id == id, ]
      check(paste(ex, o, em, id, sep = "_"), nrow(pp) == 1 && pp$valid_trials == nrow(z) && pp$rejections == sum(z$reaction == 2))
    }
  }
  # Derive IDs independently using base R; preserve original IDs in d.
  digits <- regmatches(d$participant_id, regexpr("[0-9]+", d$participant_id))
  d$participant_id_internal <- factor(sprintf("Vp%04d", as.integer(digits)))
  d$emotion <- factor(d$emotion, levels = emotion_levels)
  d$allocation <- factor(d$Offers_Other, levels = c(5, 6), labels = c("5:5", "6:4"))
  d$reject_binary <- as.integer(d$reaction == 2)
  mm <- model.matrix(~ emotion * allocation, d)[, -1, drop = FALSE]
  d[sprintf("RE%02d", seq_len(ncol(mm)))] <- mm
  m <- readRDS(file.path(run, paste0(ex, "_fair_choice.rds")))
  check(paste0(ex, "_fit_diagnostics"), isTRUE(m$sdr$pdHess) && m$fit$convergence == 0 &&
    all(vapply(VarCorr(m)$cond, function(v) attr(v, "stddev")[1], numeric(1)) >= random_sd_boundary))
  # Repeat the selected fit from independently rebuilt canonical inputs, not its saved frame.
  message("Independent selected-model refit: ", ex)
  fitted_again <- glmmTMB(formula(m), data = d, family = binomial(),
    control = glmmTMBControl(optCtrl = list(iter.max = max_evaluations, eval.max = max_evaluations)))
  diffs <- c(coefficients = max_diff(fixef(m)$cond, fixef(fitted_again)$cond),
             covariance = max_diff(vcov(m)$cond, vcov(fitted_again)$cond),
             log_likelihood = abs(as.numeric(logLik(m)) - as.numeric(logLik(fitted_again))))
  refit_comparison[[ex]] <- data.frame(experiment = ex, metric = names(diffs), maximum_absolute_difference = diffs)
  # Same machine and engine: require exact numeric reproduction.
  check(paste0(ex, "_deterministic_refit"), all(diffs == 0), paste(names(diffs), diffs, collapse = "; "))
  saveRDS(fitted_again, file.path(run, paste0(ex, "_independent_refit.rds")))
  for (label in c("primary", "actor_intercept")) {
    path <- if (label == "primary") paste0(ex, "_fair_choice.rds") else paste0(ex, "_actor_sensitivity.rds")
    if (!any(ct$experiment == ex & ct$model == label)) next
    model <- readRDS(file.path(run, path))
    beta <- fixef(model)$cond; V <- as.matrix(vcov(model)$cond)
    grid <- expand.grid(emotion = emotion_levels, allocation = c("5:5", "6:4"))
    grid$emotion <- factor(grid$emotion, levels = emotion_levels)
    grid$allocation <- factor(grid$allocation, levels = c("5:5", "6:4"))
    X <- model.matrix(~ emotion * allocation, grid)[, names(beta), drop = FALSE]
    eta <- as.numeric(X %*% beta); pr <- plogis(eta)
    sub <- ct[ct$experiment == ex & ct$model == label, ]
    for (i in seq_len(nrow(sub))) {
      parts <- strsplit(sub$contrast[i], " | ", fixed = TRUE)[[1]]
      pair <- strsplit(parts[1], " - ", fixed = TRUE)[[1]]
      w <- as.numeric(grid$emotion == pair[1]) - as.numeric(grid$emotion == pair[2])
      w <- w * if (parts[2] == "6:4 minus 5:5")
        (as.numeric(grid$allocation == "6:4") - as.numeric(grid$allocation == "5:5")) else
        as.numeric(grid$allocation == parts[2])
      if (sub$scale[i] == "log_odds") {
        gradient <- as.numeric(w %*% X); estimate <- sum(w * eta)
      } else {
        gradient <- as.numeric((w * pr * (1 - pr)) %*% X); estimate <- sum(w * pr)
      }
      se <- sqrt(as.numeric(t(gradient) %*% V %*% gradient))
      check(paste(ex, label, i, "direct_contrast_and_SE", sep = "_"),
        max_diff(c(estimate, se), c(sub$estimate[i], sub$SE[i])) < 1e-10)
      check(paste(ex, label, i, "interval", sep = "_"),
        max_diff(estimate + c(-1, 1) * qnorm(.975) * se, c(sub$asymp.LCL[i], sub$asymp.UCL[i])) < 1e-10)
    }
  }
}
for (label in unique(ct$model)) for (fam in c("focal", "secondary")) {
  sub <- ct[ct$model == label & ct$family == fam & ct$scale == "log_odds", ]
  n_family <- ifelse(fam == "focal", 12L, 6L)
  if (label == "primary") check(paste(label, fam, "complete_family"), nrow(sub) == n_family)
  p <- 2 * pnorm(-abs(sub$estimate / sub$SE))
  check(paste(label, fam, "Holm"), max_diff(p.adjust(p, "holm", n = n_family), sub$p_holm) < 1e-12)
}
write.csv(do.call(rbind, refit_comparison), file.path(run, "independent_refit_comparison.csv"), row.names = FALSE)
writeLines(c("VERIFIED: independent counts, direct linear algebra contrasts, intervals, Holm correction and exact refits.",
  "This verifies computation, not psychological interpretation or all model assumptions.",
  paste("Checks passed:", length(tests)), capture.output(sessionInfo())), file.path(run, "VERIFICATION.txt"))
message("All ", length(tests), " independent checks passed.")
