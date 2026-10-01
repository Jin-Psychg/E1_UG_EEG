#!/usr/bin/env Rscript
# Independently rebuild canonical input and exactly refit every recovered candidate.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Supply primary run and numerical retry run")
primary_run <- normalizePath(args[1],winslash="/",mustWork=TRUE)
out <- normalizePath(args[2],winslash="/",mustWork=TRUE)
if(file.exists(file.path(out,"REFIT_VERIFICATION.txt"))) stop("Verification already exists")
root <- primary_run
while(!file.exists(file.path(root,"UG_ERP_Project.Rproj"))) {
  parent <- dirname(root);if(parent==root) stop("Project root not found");root <- parent
}
.libPaths(c(file.path(root,"renv/library/R-4.3/x86_64-w64-mingw32"),.libPaths()))
suppressPackageStartupMessages({library(lme4);library(glmmTMB)})
source(file.path(out,"code/manuscript_engine.R"))
options(contrasts=c("contr.sum","contr.poly"),digits=17)
controls <- readRDS(file.path(out,"optimizer_controls.rds"))
status <- read.csv(file.path(out,"retry_status.csv"))
status <- status[status$recovered_original_gate,]
results <- list()
for(ex in unique(status$experiment)) {
  raw <- read.csv(file.path(root,paste0("data/02_Pipeline_Output_",ex,"/Method_Regression/Stimulus_Locked/trials.csv")))
  d <- raw[raw$Offers_Other %in% c(5,6)&raw$reaction %in% c(1,2)&raw$RT>=300&raw$RT<=3000,]
  d$participant_id_internal <- factor(sprintf("Vp%04d",as.integer(sub("^[^0-9]*","",d$participant_id))))
  d$actor_id <- factor(sub("^[LR]_(neu|aff|dis|dom|enj)","",d$stim))
  d$emotion <- factor(d$emotion,levels=c("neu","aff","dis","dom","enj"))
  d$allocation <- factor(d$Offers_Other,levels=c(5,6),labels=c("5:5","6:4"))
  d$reject_binary <- as.integer(d$reaction==2)
  d <- add_random_contrast_columns(d,c("emotion","allocation"),"participant_id_internal")$data
  for(label in status$model[status$experiment==ex]) {
    message("Independent recovered-model refit: ",ex," ",label)
    dd <- if(startsWith(label,"delete_")) droplevels(d[d$participant_id_internal!=sub("^delete_","",label),]) else d
    saved <- readRDS(file.path(out,paste0(ex,"_",label,".rds")))
    opt <- attr(saved,"retry_optimizer")
    m <- glmmTMB(formula(saved),data=dd,family=binomial(),REML=FALSE,control=controls[[opt]])
    delta <- c(coefficients=max(abs(fixef(saved)$cond-fixef(m)$cond)),
      covariance=max(abs(vcov(saved)$cond-vcov(m)$cond)),logLik=abs(as.numeric(logLik(saved)-logLik(m))),
      all_optimization_parameters=max(abs(saved$fit$par-m$fit$par)))
    results[[paste(ex,label)]] <- data.frame(experiment=ex,model=label,optimizer=opt,
      metric=names(delta),maximum_absolute_difference=unname(delta),valid=isTRUE(mixed_model_diagnostics(m)$valid))
    readr::write_csv(do.call(rbind,results),file.path(out,"independent_refit_comparison.csv"))
    stopifnot(all(delta==0),isTRUE(mixed_model_diagnostics(m)$valid),nobs(m)==nrow(dd))
    rm(saved,m,dd);gc(verbose=FALSE)
  }
}
if(!nrow(status)) {
  writeLines("ANALYZED: no recovered candidate available for independent valid-fit reproduction.",file.path(out,"REFIT_VERIFICATION.txt"))
} else {
  writeLines(c(paste("VERIFIED: exact independent reproduction of all",nrow(status),"recovered original-gate-valid models."),
    "Canonical input reconstruction; coefficients, covariance, log-likelihood and all optimization parameters matched exactly.",
    "Invalid attempts were not independently refitted. This does not establish full covariance rank or scientific validity.",capture.output(sessionInfo())),
    file.path(out,"REFIT_VERIFICATION.txt"))
}
