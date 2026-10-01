# Unmodified extracts from the live manuscript analysis; see engine_provenance.json.
custom_control_lmer  <- lmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 200000),
                                     calc.derivs = FALSE)
custom_control_lmer_final <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000),
  calc.derivs = TRUE
)
custom_control_glmmTMB <- glmmTMB::glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS"),
  parallel = 1
)
custom_control_glmmTMB_nlminb <- glmmTMB::glmmTMBControl(
  optimizer = nlminb,
  parallel = 1
)

# Null-coalescing helper for safe attribute fallback.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Final-model acceptance gate. Model selection may use derivative-free LMM
# fits for speed, but no model is reported or saved until the exact final
# formula has been refitted with derivatives and passes these checks.
mixed_model_diagnostics <- function(model, require_derivatives = TRUE) {
  is_glmm <- inherits(model, "glmmTMB")
  is_lmm <- inherits(model, "merMod")
  if (!is_glmm && !is_lmm) stop("Unsupported mixed-model class.")

  singular <- if (is_glmm) {
    vc <- VarCorr(model)$cond
    if (is.null(vc)) FALSE else {
      vars <- unlist(lapply(vc, function(x) diag(as.matrix(x))))
      any(!is.finite(vars)) || any(vars < 1e-6)
    }
  } else {
    lme4::isSingular(model, tol = 1e-4)
  }

  if (is_glmm) {
    optimizer_code <- model$fit$convergence %||% NA_integer_
    convergence_messages <- model$fit$message %||% character(0)
    hessian_available <- !is.null(model$sdr$pdHess)
    hessian_positive_definite <- isTRUE(model$sdr$pdHess)
    fixed_se <- suppressWarnings(sqrt(diag(vcov(model)$cond)))
  } else {
    optimizer_code <- model@optinfo$conv$opt %||% NA_integer_
    convergence_messages <- model@optinfo$conv$lme4$messages %||% character(0)
    hessian <- model@optinfo$derivs$Hessian
    hessian_available <- !is.null(hessian)
    hessian_positive_definite <- if (hessian_available) {
      eigenvalues <- tryCatch(
        eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NA_real_
      )
      all(is.finite(eigenvalues)) && min(eigenvalues) > 0
    } else {
      FALSE
    }
    fixed_se <- sqrt(diag(vcov(model)))
  }

  finite_fixed_se <- length(fixed_se) > 0 && all(is.finite(fixed_se))
  finite_likelihood <- is.finite(AIC(model)) && is.finite(as.numeric(logLik(model)))
  optimizer_ok <- isTRUE(optimizer_code == 0)
  # glmmTMB/nlminb stores a textual success message even when convergence = 0;
  # for glmmTMB, the code plus pdHess are the operative checks. For lme4,
  # distinguish hard convergence/Hessian failures from recorded soft scaling
  # warnings; the latter remain visible in REPRODUCIBILITY.txt but do not by
  # themselves invalidate a nonsingular fit with a positive-definite Hessian.
  hard_message_pattern <- paste(
    c("failed to converge", "unable to evaluate scaled gradient",
      "degenerate.*Hessian", "negative eigenvalue"),
    collapse = "|"
  )
  messages_ok <- if (is_glmm) TRUE else {
    !any(grepl(hard_message_pattern, convergence_messages,
               ignore.case = TRUE))
  }
  derivatives_ok <- !require_derivatives ||
    (hessian_available && hessian_positive_definite)
  valid <- optimizer_ok && messages_ok && !singular && finite_fixed_se &&
    finite_likelihood && derivatives_ok

  list(
    valid = valid,
    optimizer_code = optimizer_code,
    convergence_messages = paste(convergence_messages, collapse = " | "),
    singular = singular,
    hessian_available = hessian_available,
    hessian_positive_definite = hessian_positive_definite,
    finite_fixed_se = finite_fixed_se,
    finite_likelihood = finite_likelihood
  )
}

format_model_diagnostics <- function(d) {
  sprintf(
    paste0("valid=%s; optimizer_code=%s; messages=%s; singular=%s; ",
           "hessian_available=%s; hessian_positive_definite=%s; ",
           "finite_fixed_se=%s; finite_likelihood=%s"),
    d$valid, d$optimizer_code,
    ifelse(nchar(d$convergence_messages) == 0, "none", d$convergence_messages),
    d$singular, d$hessian_available, d$hessian_positive_definite,
    d$finite_fixed_se, d$finite_likelihood
  )
}

refit_final_lmm <- function(model, data) {
  attrs_to_keep <- attributes(model)[c("step4_selected", "step4_lrt_table")]
  fit_and_check <- function(candidate_formula, random_tier) {
    candidate <- lme4::lmer(
      candidate_formula, data = data, REML = TRUE,
      control = custom_control_lmer_final
    )
    diagnostics <- mixed_model_diagnostics(
      candidate, require_derivatives = TRUE
    )
    for (an in names(attrs_to_keep)) {
      if (!is.null(attrs_to_keep[[an]])) attr(candidate, an) <- attrs_to_keep[[an]]
    }
    attr(candidate, "final_validation") <- diagnostics
    attr(candidate, "final_refit_optimizer") <- "bobyqa"
    attr(candidate, "final_random_tier") <- random_tier
    candidate
  }

  exact <- fit_and_check(formula(model), "selected_formula")
  if (isTRUE(attr(exact, "final_validation")$valid)) return(exact)

  stop("The exact buildmer-selected LMM failed the acceptance gate: ",
       format_model_diagnostics(attr(exact, "final_validation")))
}

fit_glmm_candidate <- function(formula, data, control, optimizer_label) {
  model <- tryCatch(
    glmmTMB::glmmTMB(
      formula, data = data, family = binomial(link = "logit"),
      control = control, REML = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  attr(model, "final_refit_optimizer") <- optimizer_label
  attr(model, "final_validation") <- mixed_model_diagnostics(
    model, require_derivatives = TRUE
  )
  model
}

# Build a true zero-correlation-parameter representation for factor slopes.
# Separate terms such as (0 + factor | id) are not diagonal for a factor: they
# create a correlated dummy-level block and can be redundant with (1 | id).
# Here each sum-contrast column receives its own one-dimensional random term.
add_random_contrast_columns <- function(data, random_vars, group_name) {
  # Include the full factorial of the specified within-participant variables.
  # Numeric sum-contrast columns make both main and interaction ZCP blocks diagonal.
  stopifnot(length(random_vars) > 0L, !anyDuplicated(random_vars))
  if (!all(c(random_vars, group_name) %in% names(data))) {
    stop("Missing random-design or grouping variable.")
  }
  design_data <- data[random_vars]
  contrast_args <- list()
  for (v in random_vars) {
    x <- design_data[[v]]
    if (anyNA(x)) stop("Missing values in random-design variable: ", v)
    if (is.factor(x)) {
      if (nlevels(x) < 2L) stop("Random-design factor has fewer than two levels: ", v)
      contrast_args[[v]] <- contr.sum(nlevels(x))
    } else if (is.numeric(x)) {
      if (any(!is.finite(x))) stop("Non-finite random-design variable: ", v)
      design_data[[v]] <- x - mean(x)
    } else stop("Unsupported random-slope variable type: ", v)
  }
  design_formula <- reformulate(paste(random_vars, collapse = " * "))
  mm <- model.matrix(design_formula, data = design_data,
                     contrasts.arg = if (length(contrast_args)) contrast_args else NULL)
  assignments <- attr(mm, "assign")
  labels <- attr(terms(design_formula), "term.labels")
  random_terms <- random_columns <- random_blocks <- character(0)
  for (k in seq_along(labels)) {
    indices <- which(assignments == k)
    term_name <- gsub(":", "_x_", labels[k], fixed = TRUE)
    for (j in seq_along(indices)) {
      column_name <- sprintf("RE_%s_%d", term_name, j)
      data[[column_name]] <- mm[, indices[j]]
      random_columns <- c(random_columns, column_name)
      random_blocks[column_name] <- paste0(term_name, "_slope_block")
      random_terms <- c(random_terms, sprintf("(0 + %s | %s)", column_name, group_name))
    }
  }
  list(data = data, random_terms = random_terms,
       random_columns = random_columns, random_blocks = random_blocks)
}

# Buildmer's documented terms-list interface allows the k - 1 contrast columns
# belonging to one factor to be evaluated as a single block. This preserves the
# preregistered backward-LRT selection while preventing basis-dependent removal
# of individual contrasts.
make_blocked_buildmer_terms <- function(formula, random_blocks) {
  terms <- buildmer::tabulate.formula(formula)
  for (column_name in names(random_blocks)) {
    is_random_column <- !is.na(terms$grouping) & terms$term == column_name
    terms$block[is_random_column] <- unname(random_blocks[[column_name]])
  }
  terms
}

# Refit the exact buildmer-selected GLMM formula with two optimizers. No
# alternative random-effects structure is substituted at this stage: if the
# selected formula remains invalid, the analysis stops for explicit review.
stabilize_final_glmm <- function(model, data) {
  selected_formula <- formula(model)

  fit_with_optimizers <- function(candidate_formula, tier_name) {
    candidates <- list(
      BFGS = fit_glmm_candidate(candidate_formula, data,
                                custom_control_glmmTMB, "BFGS"),
      nlminb = fit_glmm_candidate(candidate_formula, data,
                                  custom_control_glmmTMB_nlminb, "nlminb")
    )
    valid <- Filter(function(x) {
      !is.null(x) && isTRUE(attr(x, "final_validation")$valid)
    }, candidates)
    if (length(valid) == 0) return(NULL)
    aics <- vapply(valid, AIC, numeric(1))
    chosen <- valid[[which.min(aics)]]
    attr(chosen, "final_random_tier") <- tier_name
    chosen
  }

  exact <- fit_with_optimizers(selected_formula, "selected_formula")
  if (!is.null(exact)) return(exact)

  stop("The exact buildmer-selected GLMM failed the acceptance gate under ",
       "both BFGS and nlminb.")
}

# ------------------------------------------------------------------------------
# Buildmer engineering helpers.
# 
# Two recurring failure modes for buildmer on large datasets:
#   (1) lmerTest, when loaded, intercepts ALL anova() calls and forces
#       Satterthwaite df computation. Inside buildmer's iterative LRT loop,
#       this triggers Hessian inversions for every candidate model. On certain
#       random structures (e.g., (0 + offer_type | id)) at large N (> 20K),
#       these Hessian inversions hang indefinitely, manifesting as CPU 0%
#       deadlock. SOLUTION: detach lmerTest before buildmer, reload after.
#   (2) Even without lmerTest, buildmer can occasionally hang on a single
#       optimization iteration. SOLUTION: wrap buildmer in withTimeout
#       (R.utils package); on timeout, stop for explicit review.
# ------------------------------------------------------------------------------

# Safe lmerTest detach helper. Returns TRUE if lmerTest was attached and we
# detached it (caller should reload after buildmer). Returns FALSE if it
# wasn't attached.
detach_lmerTest_safely <- function() {
  if ("package:lmerTest" %in% search()) {
    detach("package:lmerTest", unload = FALSE, character.only = FALSE)
    return(TRUE)
  }
  FALSE
}

# Re-attach lmerTest (only if it was detached by us).
reattach_lmerTest_safely <- function(was_detached) {
  if (isTRUE(was_detached)) {
    suppressMessages(library(lmerTest))
  }
  invisible(NULL)
}

# Wrap buildmer call in withTimeout; on timeout, return NULL so caller falls
# back to ZCP. timeout_sec defaults to 1800 (30 minutes) but is overridden
# per-stage via buildmer_timeout_map.
run_buildmer_protected <- function(buildmer_fn, ..., timeout_sec = 1800) {
  # IMPORTANT: R.utils::withTimeout interrupts R-level loops but cannot interrupt
  # blocking C-level calls (e.g., the optimizer's inner loop). For deeply stuck
  # buildmer iterations, you may still need to manually interrupt R (Ctrl+C / 
  # ESC). Timeout offers protection against MOST stuck states (R loop hangs)
  # but is not 100% guaranteed.
  was_lmerTest <- detach_lmerTest_safely()
  on.exit(reattach_lmerTest_safely(was_lmerTest), add = TRUE)
  
  result <- tryCatch({
    R.utils::withTimeout({
      buildmer_fn(...)
    }, timeout = timeout_sec, onTimeout = "error")
  }, TimeoutException = function(ex) {
    progress_log(sprintf("[TIMEOUT] buildmer exceeded %d seconds; aborting.",
                         timeout_sec), indent = 2)
    NULL
  }, error = function(e) {
    if (grepl("reached.*time.*limit|timeout", conditionMessage(e), ignore.case = TRUE)) {
      progress_log(sprintf("[TIMEOUT] buildmer exceeded %d seconds; aborting.",
                           timeout_sec), indent = 2)
      NULL
    } else {
      progress_log(sprintf("[WARN] buildmer error: %s", conditionMessage(e)),
                   indent = 2)
      NULL
    }
  })
  
  result
}

# ------------------------------------------------------------------------------
# Progress reporting helpers.
# ------------------------------------------------------------------------------
flush_progress <- function() {
  flush.console()
  if (interactive()) Sys.sleep(0)
}

progress_log <- function(msg, indent = 0) {
  prefix <- strrep("  ", indent)
  cat(sprintf("%s%s\n", prefix, msg))
  flush_progress()
}

format_elapsed <- function(t_start) {
  secs <- as.numeric(Sys.time() - t_start, units = "secs")
  if (secs < 60) sprintf("%.1fs", secs)
  else if (secs < 3600) sprintf("%dm%02ds", floor(secs / 60), round(secs %% 60))
  else sprintf("%dh%02dm", floor(secs / 3600), floor((secs %% 3600) / 60))
}
# ==============================================================================
# 4. Bates 5-step random-structure selection
# ==============================================================================
#
# Generalizes the original two-factor engine to accept an arbitrary set of
# within-subject random-slope candidates. This is necessary because:
#   E1/E2 GLMM_rejection / LMM_RT_main: random_max_vars = c("emotion", "offer_type")
#   E1/E2 LMM_RT_unfair:                random_max_vars = c("emotion", "reaction")
#   Integrative GLMM_rejection / LMM_RT_main: same as single-exp (Exp is
#                                              between-subjects, NOT in random)
#   Integrative LMM_RT_unfair:          random_max_vars = c("emotion", "reaction")
#
# CRITICAL: Exp is between-subjects in Integrative stage. Each subject has
# exactly one Exp value, so a random slope (0 + Exp | participant_id_internal)
# is mathematically unidentifiable and must NOT appear in random_max_vars.
#
# Random-effect grouping uses participant_id_internal throughout (the
# stage-appropriate uniquified id). This is set by the loader functions.
#
# RANDOM-EFFECT STRUCTURE (Bates 2015 / Matuschek 2017):
#   Step 1: Full factorial random slopes for random_max_vars, with correlations.
#           Interaction slopes are candidates; diagnostics and selection determine retention.
#   Step 2: Zero-correlation-parameter (ZCP) model
#   Step 3: LRT-based backward elimination (alpha = 0.20, Matuschek et al. 2017)
#   Step 4: Re-test correlation parameters on surviving slopes
#   Step 5: rePCA dimensionality diagnostic
# ==============================================================================

# Helper: parse a fitted model's formula to identify which random slopes
# survived Step 3.
extract_surviving_slopes <- function(model) {
  re_terms <- lme4::findbars(formula(model))
  surviving <- character(0)
  for (term in re_terms) {
    rhs <- deparse(term[[2]])
    rhs <- gsub("^\\s*0\\s*\\+\\s*", "", rhs)
    rhs <- gsub("^\\s*1\\s*\\+\\s*", "", rhs)
    rhs <- gsub("^\\s*1\\s*$",       "", rhs)
    rhs <- trimws(rhs)
    if (rhs == "" || rhs == "1") next
    parts <- trimws(strsplit(rhs, "\\+")[[1]])
    surviving <- c(surviving, parts)
  }
  unique(surviving)
}

# rePCA dimensionality summary as tidy data frame. Uses the model's actual
# random grouping name; pass group_name explicitly.
extract_repca_summary <- function(model, group_name = "participant_id_internal") {
  empty_df <- data.frame(component = NA_integer_, sd = NA_real_,
                         proportion_variance = NA_real_,
                         cumulative_variance = NA_real_)
  
  if (inherits(model, "merMod")) {
    pca_obj <- tryCatch(lme4::rePCA(model), error = function(e) NULL)
    if (is.null(pca_obj)) return(empty_df)
    # rePCA returns a list named by grouping factor. If the expected name is
    # not present, fall back to the first available grouping.
    if (!group_name %in% names(pca_obj)) {
      if (length(pca_obj) >= 1) {
        group_name <- names(pca_obj)[1]
      } else {
        return(empty_df)
      }
    }
    imp <- summary(pca_obj[[group_name]])$importance
    return(data.frame(
      component = seq_len(ncol(imp)),
      sd = imp["Standard deviation", ],
      proportion_variance = imp["Proportion of Variance", ],
      cumulative_variance = imp["Cumulative Proportion", ]
    ))
  }
  
  if (inherits(model, "glmmTMB")) {
    # glmmTMB path: extract VarCorr -> compute SVD on combined cov matrix.
    vc_list <- VarCorr(model)$cond
    if (is.null(vc_list)) return(empty_df)
    
    relevant <- vc_list[names(vc_list) == group_name |
                        grepl(paste0("^", group_name), names(vc_list))]
    if (length(relevant) == 0) {
      relevant <- vc_list
    }
    
    cov_mats <- lapply(relevant, function(x) as.matrix(x))
    total_dim <- sum(sapply(cov_mats, function(x) nrow(x)))
    big_cov <- matrix(0, total_dim, total_dim)
    pos <- 0
    for (m in cov_mats) {
      d <- nrow(m)
      big_cov[(pos + 1):(pos + d), (pos + 1):(pos + d)] <- m
      pos <- pos + d
    }
    
    sv <- tryCatch(svd(big_cov)$d, error = function(e) NULL)
    if (is.null(sv) || length(sv) == 0) return(empty_df)
    
    sv <- pmax(sv, 0)
    sds <- sqrt(sv)
    if (sum(sv) <= 0) return(empty_df)
    prop <- sv / sum(sv)
    cum_prop <- cumsum(prop)
    
    return(data.frame(
      component = seq_along(sds),
      sd = sds,
      proportion_variance = prop,
      cumulative_variance = cum_prop
    ))
  }
  
  empty_df
}

count_effective_dims <- function(repca_df, threshold = 0.999) {
  if (all(is.na(repca_df$cumulative_variance))) return(NA_integer_)
  hit <- which(repca_df$cumulative_variance >= threshold)
  if (length(hit) == 0) nrow(repca_df) else hit[1]
}

# Step 4 helper: enumerate candidate correlated structures.
# Generalized to accept arbitrary random grouping variable name.
add_correlations_step4 <- function(model_zcp_reduced, fixed_str, data,
                                   group_name = "participant_id_internal",
                                   model_family = "gaussian", alpha_lrt = 0.20) {
  if (!inherits(model_zcp_reduced, c("merMod", "glmmTMB"))) {
    progress_log("[Step 4] Reduced model is not a mixed model; skipping correlation tests.",
                 indent = 1)
    return(model_zcp_reduced)
  }
  
  surviving <- extract_surviving_slopes(model_zcp_reduced)
  
  if (length(surviving) == 0) {
    progress_log("[Step 4] Reduced model has only random intercept; no correlations to test.",
                 indent = 1)
    # SAFETY: explicit refit to clear buildmer's @call data placeholder.
    f_int_only <- formula(model_zcp_reduced)
    m_int_only <- if (model_family == "binomial") {
      glmmTMB::glmmTMB(f_int_only, data = data,
                       family = binomial(link = "logit"),
                       control = custom_control_glmmTMB,
                       REML = FALSE)
    } else {
      lme4::lmer(f_int_only, data = data, REML = TRUE,
                 control = custom_control_lmer)
    }
    attr(m_int_only, "step4_selected") <- "RandomIntercept_only"
    return(m_int_only)
  }
  
  candidates <- list()
  candidates[["ZCP_reference"]] <- formula(model_zcp_reduced)
  
  all_str <- paste(surviving, collapse = " + ")
  candidates[["All_correlated"]] <- as.formula(
    paste(fixed_str, "+ (1 +", all_str, "|", group_name, ")")
  )
  
  if (length(surviving) > 1) {
    for (i in seq_along(surviving)) {
      others <- setdiff(surviving, surviving[i])
      if (length(others) > 0) {
        others_terms <- paste(
          sprintf("(0 + %s | %s)", others, group_name),
          collapse = " + "
        )
        cand_str <- paste(fixed_str,
                          sprintf("+ (1 + %s | %s)", surviving[i], group_name),
                          "+", others_terms)
      } else {
        cand_str <- paste(fixed_str,
                          sprintf("+ (1 + %s | %s)", surviving[i], group_name))
      }
      candidates[[paste0("Int_with_", surviving[i])]] <- as.formula(cand_str)
    }
  }
  
  total_cands <- length(candidates)
  progress_log(sprintf("[Step 4] Fitting %d candidate correlated structures (%s)...",
                       total_cands, model_family), indent = 1)
  
  fits <- list()
  cand_idx <- 0
  for (nm in names(candidates)) {
    cand_idx <- cand_idx + 1
    progress_log(sprintf("[Step 4 %d/%d] Trying '%s'...",
                         cand_idx, total_cands, nm), indent = 2)
    t_cand <- Sys.time()
    
    fit <- tryCatch({
      if (model_family == "binomial") {
        glmmTMB::glmmTMB(candidates[[nm]], data = data,
                         family = binomial(link = "logit"),
                         control = custom_control_glmmTMB,
                         REML = FALSE)
      } else {
        lme4::lmer(candidates[[nm]], data = data, REML = FALSE,
                   control = custom_control_lmer_final)
      }
    }, error = function(e) NULL)
    
    fit_diag <- if (is.null(fit)) NULL else mixed_model_diagnostics(
      fit, require_derivatives = TRUE
    )
    if (!is.null(fit) && isTRUE(fit_diag$valid)) {
      fits[[nm]] <- fit
      progress_log(sprintf("[Step 4 %d/%d] '%s' converged (%s).",
                           cand_idx, total_cands, nm, format_elapsed(t_cand)),
                   indent = 3)
    } else if (!is.null(fit)) {
      progress_log(sprintf("[Step 4 %d/%d] '%s' invalid; excluded (%s): %s",
                           cand_idx, total_cands, nm, format_elapsed(t_cand),
                           format_model_diagnostics(fit_diag)),
                   indent = 3)
    } else {
      progress_log(sprintf("[Step 4 %d/%d] '%s' failed to converge (%s).",
                           cand_idx, total_cands, nm, format_elapsed(t_cand)),
                   indent = 3)
    }
  }
  
  if (is.null(fits[["ZCP_reference"]])) {
    stop(
      "[Step 4] The buildmer-selected ZCP reference failed the mandatory ",
      "acceptance gate. Correlated candidates cannot be selected without the ",
      "prespecified LRT reference; stopping for explicit review."
    )
  }
  
  if (length(fits) == 1) {
    progress_log("[Step 4] Only ZCP candidate viable; correlations not added.", indent = 1)
    f_zcp_only <- formula(model_zcp_reduced)
    if (model_family == "binomial") {
      return(glmmTMB::glmmTMB(f_zcp_only, data = data,
                              family = binomial(link = "logit"),
                              control = custom_control_glmmTMB,
                              REML = FALSE))
    } else {
      return(lme4::lmer(f_zcp_only, data = data, REML = TRUE,
                        control = custom_control_lmer_final))
    }
  }
  
  reference_fit <- fits[["ZCP_reference"]]
  lrt_table <- data.frame(model = character(), df_diff = integer(),
                          chisq = numeric(), p_value = numeric(),
                          AIC = numeric(), stringsAsFactors = FALSE)
  
  best_name <- "ZCP_reference"
  best_fit  <- reference_fit

  for (nm in setdiff(names(fits), "ZCP_reference")) {
    cand <- fits[[nm]]
    lrt <- tryCatch(anova(reference_fit, cand), error = function(e) NULL)
    if (is.null(lrt) || nrow(lrt) < 2) next
    
    p_val <- lrt$"Pr(>Chisq)"[2]
    if (is.na(p_val)) next
    
    lrt_table <- rbind(lrt_table, data.frame(
      model = nm,
      df_diff = lrt$Df[2],
      chisq = lrt$Chisq[2],
      p_value = p_val,
      AIC = AIC(cand),
      stringsAsFactors = FALSE
    ))

    if (p_val < alpha_lrt && AIC(cand) < AIC(best_fit)) {
      best_name <- nm
      best_fit  <- cand
    }
  }
  
  attr(best_fit, "step4_lrt_table") <- lrt_table
  attr(best_fit, "step4_selected") <- best_name
  
  progress_log(sprintf("[Step 4] Selected: %s (alpha_LRT = %.2f).",
                       best_name, alpha_lrt), indent = 1)
  
  if (model_family == "binomial") {
    best_fit
  } else {
    f_best <- formula(best_fit)
    m_reml <- lme4::lmer(f_best, data = data, REML = TRUE,
                         control = custom_control_lmer_final)
    for (an in c("step4_lrt_table", "step4_selected")) {
      a <- attr(best_fit, an)
      if (!is.null(a)) attr(m_reml, an) <- a
    }
    m_reml
  }
}

# ==============================================================================
# Main Bates 5-step engine. Branches on model_family. Generalized for any
# random-slope variables and arbitrary random grouping (default
# participant_id_internal).
#
# Inputs:
#   data:               clean trial-level data frame
#   outcome_var:        name of dependent variable
#   fixed_terms:        right-hand side of fixed effects (without intercept)
#                       e.g., "(emotion + offer_type)^2",
#                             "(Exp + emotion + offer_type)^3"
#   include_terms:      character vector of fixed effects to force-keep in Step 3
#   model_family:       "gaussian" (LMM) or "binomial" (GLMM)
#   random_max_vars:    within-subject factors to enter as random slopes
#   group_name:         random grouping variable (default participant_id_internal)
#   timeout_sec:        per-stage buildmer timeout (default 1800s)
# ==============================================================================
fit_stage1_actor <- function(data, outcome_var, fixed_terms, include_terms,
                             model_family = "gaussian",
                             random_max_vars = c("emotion", "offer_type"),
                             group_name = "participant_id_internal",
                             timeout_sec = 1800, actor_audit_dir = NULL) {
  
  stopifnot(model_family == "binomial", "actor_id" %in% names(data))
  fixed_str <- sprintf("%s ~ %s + (1 | actor_id)", outcome_var, fixed_terms)

  # Expand factor slopes into their k - 1 sum-contrast columns once, then use
  # those numeric columns throughout all five Bates steps. This makes the ZCP
  # model genuinely diagonal and avoids the redundant dummy-level block from
  # (1 | id) + (0 + factor | id).
  random_design <- add_random_contrast_columns(
    data, random_vars = random_max_vars, group_name = group_name
  )
  data <- random_design$data
  random_slope_vars <- random_design$random_columns
  
  # Force ALL design fixed-effect terms (plus the random intercept) into
  # buildmer's `include` so only RANDOM slopes are eligible for LRT elimination;
  # the full factorial fixed structure is retained and reported (Bates protocol
  # simplifies variance components, NOT fixed effects). NB: collapse with " + "
  # into a SINGLE formula -- paste() without collapse + as.formula() on a vector
  # silently keeps only the first term, which let buildmer drop non-significant
  # fixed effects (e.g., a non-significant 3-way interaction).
  fixed_include_str <- paste(c(include_terms, sprintf("(1 | %s)", group_name), "(1 | actor_id)"),
                             collapse = " + ")
  fixed_include <- as.formula(sprintf("~ %s", fixed_include_str))
  
  diag_log <- list(
    outcome = outcome_var,
    family = model_family,
    n_obs = nrow(data),
    n_subjects = dplyr::n_distinct(data[[group_name]]),
    random_max_vars = random_max_vars,
    random_contrast_columns = random_slope_vars,
    random_slope_blocks = unique(unname(random_design$random_blocks)),
    group_name = group_name,
    timeout_sec = timeout_sec,
    step1_status = NA_character_,
    step1_pca = NULL,
    step2_status = NA_character_,
    step3_final_formula = NA_character_,
    step4_status = NA_character_,
    step5_pca = NULL
  )
  
  random_max_str <- paste("(1 +", paste(random_slope_vars, collapse = " + "),
                          "|", group_name, ")")
  
  fit_func <- function(formula, data, REML_flag = FALSE) {
    if (model_family == "binomial") {
      glmmTMB::glmmTMB(formula, data = data,
                       family = binomial(link = "logit"),
                       control = custom_control_glmmTMB,
                       REML = REML_flag)
    } else {
      lme4::lmer(formula, data = data, REML = REML_flag,
                 control = custom_control_lmer)
    }
  }
  
  is_singular_check <- function(m) {
    if (inherits(m, "glmmTMB")) {
      vc <- VarCorr(m)$cond
      if (is.null(vc)) return(FALSE)
      vars <- unlist(lapply(vc, function(x) diag(as.matrix(x))))
      any(vars < 1e-6)
    } else {
      lme4::isSingular(m, tol = 1e-4)
    }
  }
  
  # ---- STEP 1: Barr-style maximal model ----
  progress_log(sprintf("[%s | %s] Stage 1: Bates 5-step random-structure selection",
                       outcome_var, model_family))
  progress_log(sprintf("Data: N = %d trials, %d subjects (group = %s).",
                       nrow(data), dplyr::n_distinct(data[[group_name]]),
                       group_name), indent = 1)
  progress_log(sprintf("Fixed terms: %s", fixed_terms), indent = 1)
  progress_log(sprintf("Forced include: %s",
                       paste(include_terms, collapse = ", ")), indent = 1)
  progress_log(sprintf("Random max vars: %s",
                       paste(random_max_vars, collapse = ", ")), indent = 1)
  progress_log(sprintf("Random contrast columns: %s",
                       paste(random_slope_vars, collapse = ", ")), indent = 1)
  progress_log(sprintf("buildmer timeout: %d sec", timeout_sec), indent = 1)
  maximal_diagnostic_timeout_sec <- min(timeout_sec, 300)
  progress_log(sprintf("maximal diagnostic timeout: %d sec",
                       maximal_diagnostic_timeout_sec), indent = 1)
  
  f_max <- as.formula(paste(fixed_str, "+", random_max_str))
  diag_log$step1_formula <- paste(deparse(f_max), collapse = " ")
  progress_log(diag_log$step1_formula, indent = 2)
  
  progress_log("[Step 1] Fitting design-justified maximal model with correlations...",
               indent = 1)
  t_step1 <- Sys.time()
  m_max <- tryCatch(
    R.utils::withTimeout(
      fit_func(f_max, data, REML_flag = FALSE),
      timeout = maximal_diagnostic_timeout_sec,
      onTimeout = "error"
    ),
    error = function(e) {
      progress_log(sprintf("[WARN] Maximal diagnostic failed or timed out: %s",
                           conditionMessage(e)),
                   indent = 2)
      NULL
    }
  )
  progress_log(sprintf("[Step 1] elapsed: %s", format_elapsed(t_step1)),
               indent = 2)
  
  if (!is.null(m_max)) {
    is_sing <- is_singular_check(m_max)
    step1_diag <- mixed_model_diagnostics(
      m_max, require_derivatives = model_family == "binomial"
    )
    diag_log$step1_status <- if (!isTRUE(step1_diag$valid)) {
      paste0("invalid: ", format_model_diagnostics(step1_diag))
    } else if (is_sing) {
      "converged_singular"
    } else {
      "converged_clean"
    }
    diag_log$step1_pca <- extract_repca_summary(m_max, group_name = group_name)
    progress_log(sprintf("[Step 1] Status: %s. Effective dims = %d / %d.",
                         diag_log$step1_status,
                         count_effective_dims(diag_log$step1_pca),
                         nrow(diag_log$step1_pca)),
                 indent = 2)
  } else {
    diag_log$step1_status <- "diagnostic_failed_or_timed_out"
  }
  
  # ---- STEP 2: ZCP model ----
  progress_log("[Step 2] Fitting zero-correlation-parameter model...", indent = 1)
  t_step2 <- Sys.time()
  
  zcp_components <- c(
    sprintf("(1 | %s)", group_name),
    random_design$random_terms
  )
  random_zcp <- paste(zcp_components, collapse = " + ")
  f_zcp <- as.formula(paste(fixed_str, "+", random_zcp))
  diag_log$step2_formula <- paste(deparse(f_zcp), collapse = " ")
  progress_log(diag_log$step2_formula, indent = 2)
  blocked_terms <- make_blocked_buildmer_terms(
    f_zcp, random_design$random_blocks
  )
  
  m_zcp <- tryCatch(fit_func(f_zcp, data, REML_flag = FALSE),
                    error = function(e) NULL)
  progress_log(sprintf("[Step 2] elapsed: %s", format_elapsed(t_step2)), indent = 2)

  step2_diag <- if (is.null(m_zcp)) NULL else mixed_model_diagnostics(
    m_zcp, require_derivatives = model_family == "binomial"
  )
  if (model_family == "binomial" &&
      (is.null(step2_diag) || !isTRUE(step2_diag$valid))) {
    progress_log("[Step 2] BFGS ZCP invalid; retrying exact formula with nlminb.",
                 indent = 2)
    m_zcp <- tryCatch(
      glmmTMB::glmmTMB(
        f_zcp, data = data, family = binomial(link = "logit"),
        control = custom_control_glmmTMB_nlminb, REML = FALSE
      ),
      error = function(e) NULL
    )
    step2_diag <- if (is.null(m_zcp)) NULL else mixed_model_diagnostics(
      m_zcp, require_derivatives = TRUE
    )
  }
  if (!is.null(actor_audit_dir) && !is.null(m_zcp)) saveRDS(m_zcp, file.path(actor_audit_dir, "step2_zcp.rds"))
  diag_log$step2_status <- if (is.null(step2_diag)) {
    "fit_failed"
  } else if (isTRUE(step2_diag$valid)) {
    "converged_clean"
  } else {
    paste0("diagnostic_only_invalid: ",
           format_model_diagnostics(step2_diag))
  }
  progress_log(sprintf("[Step 2] Status: %s.", diag_log$step2_status), indent = 2)
  
  # ---- STEP 3: buildmer LRT-based random-structure simplification ----
  progress_log("[Step 3] LRT-based backward elimination (buildmer, alpha = 0.20)...",
               indent = 1)
  progress_log(paste("This step iteratively tests model terms via LRT and may fit",
                     "many candidate models. With large N or complex random",
                     "structures, expect several minutes."),
               indent = 2)
  t_step3 <- Sys.time()
  
  if (model_family == "binomial") {
    progress_log("[Step 3] Engine: buildglmmTMB (glmmTMB-based LRT elimination).",
                 indent = 2)
    bm_ctrl <- buildmer::buildmerControl(
      dep = outcome_var,
      include = fixed_include,
      direction = c("order", "backward"),
      crit = "LRT",
      elim = buildmer::LRTalpha(0.20),
      quiet = FALSE,
      ddf = "Wald",
      args = list(control = custom_control_glmmTMB)
    )
    
    bm_res <- run_buildmer_protected(
      buildmer::buildglmmTMB,
      formula = blocked_terms, data = data,
      family = binomial(link = "logit"),
      buildmerControl = bm_ctrl,
      timeout_sec = timeout_sec
    )
    
    diag_log$step3_method <- "buildglmmTMB_LRT"
    
  } else {
    progress_log("[Step 3] Engine: buildmer (lme4-based LRT elimination).", indent = 2)
    bm_ctrl <- buildmer::buildmerControl(
      dep = outcome_var,
      include = fixed_include,
      direction = c("order", "backward"),
      crit = "LRT",
      elim = buildmer::LRTalpha(0.20),
      quiet = FALSE,
      args = list(control = custom_control_lmer)
    )
    
    bm_res <- run_buildmer_protected(
      buildmer::buildmer,
      formula = blocked_terms, data = data,
      buildmerControl = bm_ctrl,
      timeout_sec = timeout_sec
    )
    
    diag_log$step3_method <- "buildmer_LRT"
  }
  
  if (is.null(bm_res)) {
    stop("Step 3 buildmer selection failed or timed out; no unvalidated ",
         "ZCP fallback will be reported.")
  } else {
    m_reduced <- tryCatch({
      bm_res@model
    }, error = function(e) {
      stop("Could not extract the buildmer-selected model: ",
           conditionMessage(e))
    })
    diag_log$step3_fallback <- NA_character_
  }
  
  progress_log(sprintf("[Step 3] elapsed: %s", format_elapsed(t_step3)), indent = 2)
  
  if (!inherits(m_reduced, c("merMod", "glmmTMB"))) {
    progress_log("[WARN] Reduced model lost all random effects; refitting with intercept-only.",
                 indent = 2)
    reduced_fixed_str_recovered <- paste(deparse(formula(m_reduced)), collapse = " ")
    f_recover <- as.formula(paste(reduced_fixed_str_recovered,
                                  sprintf("+ (1 | %s)", group_name)))
    m_reduced <- fit_func(f_recover, data, REML_flag = FALSE)
  }
  
  if (!is.null(actor_audit_dir)) saveRDS(m_reduced, file.path(actor_audit_dir, "step3_selected.rds"))
  diag_log$step3_final_formula <- paste(deparse(formula(m_reduced)), collapse = " ")
  progress_log(sprintf("[Step 3] Reduced formula: %s",
                       stringr::str_squish(diag_log$step3_final_formula)),
               indent = 2)
  
  # ---- STEP 4: Test correlation parameters on surviving slopes ----
  progress_log("[Step 4] Testing correlation parameters on surviving slopes...",
               indent = 1)
  t_step4 <- Sys.time()
  
  reduced_fixed <- nobars(formula(m_reduced))
  reduced_fixed_str <- paste(paste(deparse(reduced_fixed), collapse = " "), "+ (1 | actor_id)")
  
  m_final <- add_correlations_step4(m_reduced, reduced_fixed_str, data,
                                    group_name = group_name,
                                    model_family = model_family,
                                    alpha_lrt = 0.20)
  if (!is.null(actor_audit_dir)) saveRDS(m_final, file.path(actor_audit_dir, "step4_selected.rds"))
  diag_log$step4_status <- attr(m_final, "step4_selected") %||% "ZCP_kept"
  progress_log(sprintf("[Step 4] elapsed: %s. Selected: %s.",
                       format_elapsed(t_step4), diag_log$step4_status),
               indent = 2)
  
  # ---- STEP 5: rePCA diagnostic ----
  diag_log$step5_pca <- extract_repca_summary(m_final, group_name = group_name)
  step5_eff <- count_effective_dims(diag_log$step5_pca)
  step5_nom <- nrow(diag_log$step5_pca)
  progress_log(sprintf("[Step 5] Final model dimensionality: %s effective / %s nominal.",
                       if (is.na(step5_eff)) "NA" else as.character(step5_eff),
                       if (is.na(step5_nom)) "NA" else as.character(step5_nom)),
               indent = 1)
  
  # ---- Defensive REML guard for LMM ----
  if (model_family == "gaussian" && inherits(m_final, "merMod")) {
    is_reml <- tryCatch(lme4::isREML(m_final), error = function(e) FALSE)
    if (!isTRUE(is_reml)) {
      progress_log("[Final] Refitting LMM with REML=TRUE for Satterthwaite compatibility.",
                   indent = 1)
      m_final_reml <- tryCatch(
        update(m_final, REML = TRUE, control = custom_control_lmer),
        error = function(e) {
          progress_log(sprintf("[Final] REML refit failed: %s; retaining ML fit.",
                               conditionMessage(e)), indent = 2)
          NULL
        }
      )
      if (!is.null(m_final_reml)) {
        existing_attrs <- attributes(m_final)[c("step4_selected",
                                                "step4_lrt_table")]
        for (an in names(existing_attrs)) {
          if (!is.null(existing_attrs[[an]])) {
            attr(m_final_reml, an) <- existing_attrs[[an]]
          }
        }
        m_final <- m_final_reml
      }
    }
  }

  # ---- Mandatory final-model validation ----
  # Selection fits above remain derivative-free for LMM speed. Here the exact
  # selected LMM is refitted with derivatives. GLMMs are re-estimated under a
  # bounded, fixed-effect-preserving optimizer/random-structure protocol.
  if (model_family == "gaussian" && inherits(m_final, "merMod")) {
    progress_log("[Final] Refitting exact LMM formula with derivatives enabled.",
                 indent = 1)
    m_final <- refit_final_lmm(m_final, data)
  } else if (model_family == "binomial" && inherits(m_final, "glmmTMB")) {
    progress_log("[Final] Validating the exact buildmer-selected GLMM formula.",
                 indent = 1)
    existing_attrs <- attributes(m_final)[c("step4_selected", "step4_lrt_table")]
    m_final <- stabilize_final_glmm(m_final, data)
    for (an in names(existing_attrs)) {
      if (!is.null(existing_attrs[[an]])) attr(m_final, an) <- existing_attrs[[an]]
    }
  }

  final_validation <- attr(m_final, "final_validation")
  diag_log$final_validation <- final_validation
  diag_log$final_refit_optimizer <- attr(m_final, "final_refit_optimizer") %||% NA_character_
  diag_log$final_random_tier <- attr(m_final, "final_random_tier") %||%
    if (model_family == "gaussian") "selected_formula" else NA_character_
  diag_log$step5_pca <- extract_repca_summary(m_final, group_name = group_name)
  progress_log(sprintf("[Final] Acceptance gate: %s",
                       format_model_diagnostics(final_validation)), indent = 1)
  progress_log(sprintf("[Final] Optimizer: %s; random tier: %s.",
                       diag_log$final_refit_optimizer,
                       diag_log$final_random_tier), indent = 1)
  
  # ---- Defensive lmerTest class upgrade for LMM ----
  if (model_family == "gaussian" && inherits(m_final, "merMod") &&
      !inherits(m_final, "lmerModLmerTest")) {
    m_final_lt <- tryCatch(
      lmerTest::as_lmerModLmerTest(m_final),
      error = function(e) {
        progress_log(sprintf("[Final] lmerTest class upgrade failed: %s",
                             conditionMessage(e)), indent = 2)
        NULL
      }
    )
    if (!is.null(m_final_lt)) {
      existing_attrs <- attributes(m_final)[c(
        "step4_selected", "step4_lrt_table", "final_validation",
        "final_refit_optimizer", "final_random_tier"
      )]
      for (an in names(existing_attrs)) {
        if (!is.null(existing_attrs[[an]])) {
          attr(m_final_lt, an) <- existing_attrs[[an]]
        }
      }
      m_final <- m_final_lt
      progress_log("[Final] LMM upgraded to lmerModLmerTest for Satt/KR dispatch.",
                   indent = 1)
    }
  }
  
  attr(m_final, "bates_diag") <- diag_log
  m_final
}
