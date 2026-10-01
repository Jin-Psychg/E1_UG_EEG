#!/usr/bin/env Rscript
# Independent reconstruction and exact full-sample refits; does not reselect models.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Supply primary run and actor run with both selections resolved")
primary <- normalizePath(args[1],winslash="/",mustWork=TRUE)
out <- normalizePath(args[2],winslash="/",mustWork=TRUE)
if(file.exists(file.path(out,"REFIT_VERIFICATION.txt"))) stop("Preserve verification")
root <- dirname(dirname(primary))
.libPaths(c(file.path(root,"renv/library/R-4.3/x86_64-w64-mingw32"),.libPaths()))
suppressPackageStartupMessages({library(lme4);library(glmmTMB)})
source(file.path(out,"code/actor_engine.R"))
options(contrasts=c("contr.sum","contr.poly"),digits=17)
controls <- list(BFGS=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),optCtrl=list(maxit=2000),parallel=1),
  nlminb=glmmTMBControl(optimizer=nlminb,optCtrl=list(iter.max=2000,eval.max=4000),parallel=1))
status <- read.csv(file.path(out,"selection_status.csv"))
stopifnot(nrow(status)==2L,setequal(status$experiment,c("E1","E2")))
results <- stages <- correlation_tests <- list()
for(ex in status$experiment) {
  for(path in list.files(file.path(out,paste0(ex,"_selection_stages")),pattern="[.]rds$",full.names=TRUE)) {
    stage <- readRDS(path)
    f <- formula(stage)
    random_terms <- vapply(lme4::findbars(f),function(b)paste(deparse(b),collapse=" "),character(1))
    stopifnot("1 | actor_id" %in% random_terms,
      setequal(attr(terms(lme4::nobars(f)),"term.labels"),c("emotion","allocation","emotion:allocation")))
    stages[[paste(ex,basename(path))]] <- data.frame(experiment=ex,stage=basename(path),
      formula=paste(deparse(f),collapse=" "),diagnostics=format_model_diagnostics(mixed_model_diagnostics(stage)))
    tests <- attr(stage,"step4_lrt_table")
    if(!is.null(tests)) {tests$experiment <- ex;correlation_tests[[ex]] <- tests}
  }
  raw <- read.csv(file.path(root,paste0("data/02_Pipeline_Output_",ex,"/Method_Regression/Stimulus_Locked/trials.csv")))
  d <- raw[raw$Offers_Other %in% c(5,6)&raw$reaction %in% c(1,2)&raw$RT>=300&raw$RT<=3000,]
  d$participant_id_internal <- factor(sprintf("Vp%04d",as.integer(sub("^[^0-9]*","",d$participant_id))))
  d$actor_id <- factor(sub("^[LR]_(neu|aff|dis|dom|enj)","",d$stim))
  d$emotion <- factor(d$emotion,levels=c("neu","aff","dis","dom","enj"))
  d$allocation <- factor(d$Offers_Other,levels=c(5,6),labels=c("5:5","6:4"))
  d$reject_binary <- as.integer(d$reaction==2)
  exported <- read.csv(file.path(primary,paste0(ex,"_analysis_trials.csv")))
  for(n in names(exported)[names(exported) %in% names(d)]) {
    actual <- if(is.factor(d[[n]])) as.character(d[[n]]) else d[[n]]
    stopifnot(isTRUE(all.equal(actual,exported[[n]],check.attributes=FALSE)))
  }
  d <- add_random_contrast_columns(d,c("emotion","allocation"),"participant_id_internal")$data
  s <- status[status$experiment==ex,]
  labels <- c(if(s$actor_valid) "actor_selected",if(s$matched_pair_valid) "matched_no_actor")
  for(label in labels) {
    message("Independent full-sample refit: ",ex," ",label)
    saved <- readRDS(file.path(out,paste0(ex,"_",label,".rds")))
    opt <- attr(saved,if(label=="actor_selected") "final_refit_optimizer" else "robustness_optimizer")
    stopifnot(opt %in% names(controls))
    m <- glmmTMB(formula(saved),data=d,family=binomial(),REML=FALSE,control=controls[[opt]])
    delta <- c(coefficients=max(abs(fixef(saved)$cond-fixef(m)$cond)),
      covariance=max(abs(vcov(saved)$cond-vcov(m)$cond)),logLik=abs(as.numeric(logLik(saved)-logLik(m))),
      all_optimization_parameters=max(abs(saved$fit$par-m$fit$par)))
    results[[paste(ex,label)]] <- data.frame(experiment=ex,model=label,optimizer=opt,
      metric=names(delta),maximum_absolute_difference=unname(delta),valid=isTRUE(mixed_model_diagnostics(m)$valid))
    readr::write_csv(do.call(rbind,results),file.path(out,"independent_refit_comparison.csv"))
    stopifnot(all(delta==0),isTRUE(mixed_model_diagnostics(m)$valid),nobs(m)==nrow(d))
    rm(saved,m);gc(verbose=FALSE)
  }
}
readr::write_csv(do.call(rbind,stages),file.path(out,"independent_stage_audit.csv"))
if(length(correlation_tests)) readr::write_csv(do.call(rbind,correlation_tests),file.path(out,"correlation_selection.csv"))
writeLines(c(paste("VERIFIED: canonical data reconstructed for both experiments; exact full-sample refits:",length(results)),
  "Selection, invalid attempts and deletion fits were not independently repeated. Numerical reproduction does not establish scientific robustness.",
  capture.output(sessionInfo())),file.path(out,"REFIT_VERIFICATION.txt"))
