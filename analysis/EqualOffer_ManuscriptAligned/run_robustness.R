#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply primary run and NEW output directory")
run <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
if (dir.exists(out) || file.exists(out)) stop("Output already exists; preserve it")
root <- run
while (!file.exists(file.path(root, "UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if (parent == root) stop("Project root not found"); root <- parent
}
if (!startsWith(out, paste0(root, "/results/"))) stop("Output must be in project results")
.libPaths(c(file.path(root, "renv/library/R-4.3/x86_64-w64-mingw32"), .libPaths()))
suppressPackageStartupMessages({library(lme4); library(glmmTMB); library(emmeans)})
code <- file.path(root, "analysis/EqualOffer_ManuscriptAligned")
source(file.path(code, "manuscript_engine.R"))
options(contrasts = c("contr.sum", "contr.poly"), digits = 17)
if (!identical(readLines(file.path(run, "STATUS.txt"))[1], "COMPLETE")) stop("Primary run not complete")
provenance <- jsonlite::read_json(file.path(code, "engine_provenance.json"))
stopifnot(digest::digest(file=file.path(code,"manuscript_engine.R"),algo="sha256") == provenance$snapshot_sha256,
          digest::digest(file=file.path(root,"Sta_Behaviour_E1_E2_Integrative.Rmd"),algo="sha256") == provenance$source_sha256)
dir.create(out)
dir.create(file.path(out,"code"))
file.copy(file.path(code,c("run_robustness.R","ROBUSTNESS_PLAN.md","manuscript_engine.R","engine_provenance.json")),file.path(out,"code"))
writeLines("RUNNING",file.path(out,"STATUS.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
save_csv <- function(x,name) readr::write_csv(x,file.path(out,name),na="NA")
read <- function(name) read.csv(file.path(run,name),check.names=FALSE,stringsAsFactors=FALSE)
source_paths <- c(list.files(run,full.names=TRUE,recursive=TRUE),
  file.path(root,paste0("data/02_Pipeline_Output_",c("E1","E2"),"/Method_Regression/Stimulus_Locked/trials.csv")))
source_paths <- source_paths[!dir.exists(source_paths)]
hashes <- data.frame(path=source_paths,md5=unname(tools::md5sum(source_paths)))
save_csv(hashes,"input_manifest.csv")
attempts <- covariances <- effects <- omnibus <- sims <- checks <- list()
check <- function(label,value) {
  checks[[length(checks)+1L]] <<- data.frame(check=label,passed=isTRUE(value))
  save_csv(do.call(rbind,checks),"verification_checks.csv")
  if (!isTRUE(value)) stop("Check failed: ",label)
}
covariance_rows <- function(m,ex,label,opt) {
  vv <- VarCorr(m)$cond
  do.call(rbind,lapply(seq_along(vv),function(j) {
    ev <- eigen(vv[[j]],symmetric=TRUE,only.values=TRUE)$values
    cr <- attr(vv[[j]],"correlation"); cr <- cr[upper.tri(cr)]
    data.frame(experiment=ex,model=label,optimizer=opt,block=j,group=names(vv)[j],
      min_variance=min(diag(vv[[j]])),min_eigenvalue=min(ev),max_eigenvalue=max(ev),
      relative_min_eigenvalue=min(ev)/max(ev),max_abs_correlation=if(length(cr)) max(abs(cr)) else 0)
  }))
}
fit_exact <- function(f,d,ex,label) {
  candidates <- list()
  for (opt in c("BFGS","nlminb")) {
    warnings <- character(); err <- ""
    start <- proc.time()[[3]]
    m <- tryCatch(withCallingHandlers(R.utils::withTimeout(glmmTMB(f,data=d,
      family=binomial(),REML=FALSE,control=if(opt=="BFGS") custom_control_glmmTMB else custom_control_glmmTMB_nlminb),
      timeout=300,onTimeout="error"),warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")}),
      error=function(e) {err <<- conditionMessage(e); NULL})
    dg <- if(is.null(m)) NULL else mixed_model_diagnostics(m)
    valid <- !is.null(dg) && isTRUE(dg$valid)
    attempts[[length(attempts)+1L]] <<- data.frame(experiment=ex,model=label,optimizer=opt,
      n_trials=nrow(d),n_participants=length(unique(d$participant_id_internal)),valid=valid,
      AIC=if(is.null(m)) NA_real_ else AIC(m),elapsed_seconds=proc.time()[[3]]-start,
      error=err,warnings=paste(unique(warnings),collapse=" | "),
      diagnostics=if(is.null(dg)) "fit_failed" else format_model_diagnostics(dg))
    save_csv(do.call(rbind,attempts),"fit_attempts.csv")
    if (!is.null(m)) {
      covariances[[length(covariances)+1L]] <<- covariance_rows(m,ex,label,opt)
      save_csv(do.call(rbind,covariances),"covariance_diagnostics.csv")
    }
    if(valid) candidates[[opt]] <- m
  }
  if(!length(candidates)) return(NULL)
  opt <- names(which.min(vapply(candidates,AIC,numeric(1))))
  m <- candidates[[opt]]; attr(m,"robustness_optimizer") <- opt
  saveRDS(m,file.path(out,paste0(ex,"_",label,".rds")))
  m
}
extract <- function(m,ex,label) {
  ct <- as.data.frame(summary(pairs(emmeans(m,~emotion,weights="equal")),infer=c(TRUE,TRUE),adjust="none"))
  ct$p_BH <- p.adjust(ct$p.value,"BH"); ct$family_n <- 10L; ct$effect_type <- "emotion_marginal"
  al <- as.data.frame(summary(pairs(emmeans(m,~allocation,weights="equal")),infer=c(TRUE,TRUE),adjust="none"))
  al$p_BH <- al$p.value; al$family_n <- 1L; al$effect_type <- "allocation_marginal"
  z <- rbind(ct,al); z$experiment <- ex; z$model <- label
  z$OR <- exp(z$estimate); z$OR_lower <- exp(z$asymp.LCL); z$OR_upper <- exp(z$asymp.UCL)
  effects[[length(effects)+1L]] <<- z
  save_csv(do.call(rbind,effects),"marginal_effects.csv")
  om <- as.data.frame(car::Anova(m,type=3)); om$term <- rownames(om); rownames(om) <- NULL
  om$experiment <- ex; om$model <- label
  omnibus[[length(omnibus)+1L]] <<- om
  save_csv(do.call(rbind,omnibus),"omnibus.csv")
  # Independent fixed-design arithmetic for every extracted marginal contrast.
  grid <- expand.grid(emotion=factor(c("neu","aff","dis","dom","enj"),levels=c("neu","aff","dis","dom","enj")),
    allocation=factor(c("5:5","6:4"),levels=c("5:5","6:4")))
  b <- fixef(m)$cond; V <- vcov(m)$cond
  X <- model.matrix(~emotion*allocation,grid)[,names(b),drop=FALSE]
  for(i in seq_len(nrow(z))) {
    pair <- strsplit(z$contrast[i]," - ",fixed=TRUE)[[1]]
    v <- if(z$effect_type[i]=="emotion_marginal") "emotion" else "allocation"
    w <- (as.numeric(grid[[v]]==pair[1])-as.numeric(grid[[v]]==pair[2])) / if(v=="emotion") 2 else 5
    g <- as.numeric(w %*% X); est <- sum(g*b); se <- sqrt(as.numeric(t(g)%*%V%*%g))
    check(paste(ex,label,i,"direct_estimate_SE"),max(abs(c(est,se)-c(z$estimate[i],z$SE[i])))<1e-10)
  }
}
for(ex in c("E1","E2")) {
  message("Starting ",ex)
  d <- read(paste0(ex,"_analysis_trials.csv"))
  check(paste(ex,"required_columns"),all(c("participant_id","participant_id_internal","index","actor_id","reaction","reject_binary","emotion","allocation","RT") %in% names(d)))
  check(paste(ex,"complete_valid_input"),!anyNA(d) && all(d$reject_binary %in% 0:1) &&
    all(d$reject_binary==as.integer(d$reaction==2)) && all(d$RT>=300 & d$RT<=3000) && !anyDuplicated(paste(d$participant_id,d$index)))
  raw <- read.csv(file.path(root,paste0("data/02_Pipeline_Output_",ex,"/Method_Regression/Stimulus_Locked/trials.csv")))
  keep <- raw$Offers_Other %in% c(5,6) & raw$reaction %in% c(1,2) & raw$RT>=300 & raw$RT<=3000
  raw <- raw[keep,]
  check(paste(ex,"canonical_rows_outcomes_actor"),identical(paste(raw$participant_id,raw$index),paste(d$participant_id,d$index)) &&
    identical(as.integer(raw$reaction==2),d$reject_binary) && identical(sub("^[LR]_(neu|aff|dis|dom|enj)","",raw$stim),d$actor_id))
  d$emotion <- factor(d$emotion,levels=c("neu","aff","dis","dom","enj"))
  d$allocation <- factor(d$allocation,levels=c("5:5","6:4"))
  d$participant_id_internal <- factor(d$participant_id_internal); d$actor_id <- factor(d$actor_id)
  d <- add_random_contrast_columns(d,c("emotion","allocation"),"participant_id_internal")$data
  primary <- readRDS(file.path(run,paste0(ex,"_fair_choice.rds")))
  frame <- model.frame(primary)
  check(paste(ex,"primary_frame"),all(vapply(names(frame),function(n) isTRUE(all.equal(frame[[n]],d[[n]],check.attributes=FALSE)),logical(1))))
  extract(primary,ex,"primary")
  baseline <- effects[[length(effects)]]
  ref <- read("emotion_main_contrasts.csv"); ref <- ref[ref$experiment==ex,]
  check(paste(ex,"saved_primary_effects"),identical(ref$contrast,baseline$contrast[1:10]) &&
    max(abs(ref$estimate-baseline$estimate[1:10]))<1e-12 && max(abs(ref$p_BH-baseline$p_BH[1:10]))<1e-12)
  covariances[[length(covariances)+1L]] <- covariance_rows(primary,ex,"primary",attr(primary,"final_refit_optimizer"))
  save_csv(do.call(rbind,covariances),"covariance_diagnostics.csv")
  participants <- levels(d$participant_id_internal)
  save_csv(data.frame(participant_id_internal=participants),paste0(ex,"_deletion_roster.csv"))
  for(id in participants) {
    label <- paste0("delete_",id); message(ex," ",label)
    dd <- droplevels(d[d$participant_id_internal!=id,])
    m <- fit_exact(formula(primary),dd,ex,label)
    if(!is.null(m)) extract(m,ex,label)
    rm(m,dd); gc(verbose=FALSE)
  }
  message(ex," actor intercept sensitivity")
  actor <- fit_exact(update(formula(primary),. ~ . + (1|actor_id)),d,ex,"actor_intercept")
  if(!is.null(actor)) extract(actor,ex,"actor_intercept")
  # Fixed-parameter, new-random-effects simulation; no model refitting.
  nsim <- 1000L; seed <- if(ex=="E1") 20260928L else 20260929L
  draws <- simulate(primary,nsim=nsim,seed=seed)
  grouping <- interaction(d$allocation,d$emotion,drop=TRUE)
  summarize_draw <- function(y) {
    if(is.matrix(y)) {stopifnot(ncol(y)==2L,all(rowSums(y)==1)); y <- y[,1]}
    stopifnot(length(y)==nrow(d),all(y %in% 0:1))
    totals <- tapply(y,d$participant_id_internal,sum)
    c(setNames(as.numeric(tapply(y,grouping,sum)),paste0("cell_",levels(grouping))),
      total_rejections=sum(y),zero_rejection_participants=sum(totals==0),
      max_participant_share=if(sum(y)>0) max(totals)/sum(y) else NA_real_)
  }
  observed <- summarize_draw(d$reject_binary)
  simulated <- t(vapply(draws,summarize_draw,numeric(length(observed))))
  colnames(simulated) <- names(observed)
  save_csv(data.frame(simulation=seq_len(nsim),simulated,check.names=FALSE),paste0(ex,"_simulation_draws.csv"))
  summary <- data.frame(experiment=ex,metric=names(observed),observed=as.numeric(observed),
    lower=apply(simulated,2,quantile,probs=.025,na.rm=TRUE),median=apply(simulated,2,median,na.rm=TRUE),
    upper=apply(simulated,2,quantile,probs=.975,na.rm=TRUE),n_simulations=nsim,seed=seed)
  summary$outside_pointwise_envelope <- summary$observed<summary$lower | summary$observed>summary$upper
  sims[[ex]] <- summary
  save_csv(do.call(rbind,sims),"simulation_summary.csv")
  check(paste(ex,"simulation_cell_totals"),all(abs(rowSums(simulated[,seq_len(10)])-simulated[,"total_rejections"])==0))
  rm(primary,actor,draws,simulated,d,raw); gc(verbose=FALSE)
}
check("input_hashes_unchanged",identical(unname(tools::md5sum(hashes$path)),hashes$md5))
writeLines("COMPLETE: numerical checks executed; interpretation in ROBUSTNESS_REPORT.md",file.path(out,"STATUS.txt"))
message("Completed: ",out)
