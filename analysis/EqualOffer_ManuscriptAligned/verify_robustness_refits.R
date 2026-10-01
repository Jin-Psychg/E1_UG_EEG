#!/usr/bin/env Rscript
# Reproduce the largest standardized-change deletion and valid actor fit in each experiment.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2L) stop("Supply primary run and completed robustness run")
primary_run <- normalizePath(a[1],winslash="/",mustWork=TRUE)
out <- normalizePath(a[2],winslash="/",mustWork=TRUE)
dest <- file.path(out,"REFIT_VERIFICATION.txt")
if(file.exists(dest)) stop("Refit verification exists; preserve it")
root <- primary_run
while(!file.exists(file.path(root,"UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if(parent==root) stop("Project root not found"); root <- parent
}
.libPaths(c(file.path(root,"renv/library/R-4.3/x86_64-w64-mingw32"),.libPaths()))
suppressPackageStartupMessages({library(lme4);library(glmmTMB)})
source(file.path(out,"code/manuscript_engine.R"))
options(contrasts=c("contr.sum","contr.poly"),digits=17)
selection <- read.csv(file.path(out,"selected_optimizers.csv"))
influence <- read.csv(file.path(out,"influence_summary.csv"))
comparison <- list()
for(ex in c("E1","E2")) {
  # Rebuild from canonical trial inputs, not the exported analysis data.
  raw <- read.csv(file.path(root,paste0("data/02_Pipeline_Output_",ex,"/Method_Regression/Stimulus_Locked/trials.csv")))
  d <- raw[raw$Offers_Other %in% c(5,6) & raw$reaction %in% c(1,2) & raw$RT>=300 & raw$RT<=3000,]
  d$participant_id_internal <- factor(sprintf("Vp%04d",as.integer(sub("^[^0-9]*","",d$participant_id))))
  d$actor_id <- factor(sub("^[LR]_(neu|aff|dis|dom|enj)","",d$stim))
  d$emotion <- factor(d$emotion,levels=c("neu","aff","dis","dom","enj"))
  d$allocation <- factor(d$Offers_Other,levels=c(5,6),labels=c("5:5","6:4"))
  d$reject_binary <- as.integer(d$reaction==2)
  d <- add_random_contrast_columns(d,c("emotion","allocation"),"participant_id_internal")$data
  z <- influence[influence$experiment==ex,]
  labels <- if(nrow(z)) z$most_influential_model[which.max(z$max_abs_change_in_primary_SE)] else character()
  if(any(selection$experiment==ex & selection$model=="actor_intercept")) labels <- c(labels,"actor_intercept")
  for(label in labels) {
    message("Verification refit: ",ex," ",label)
    dd <- if(startsWith(label,"delete_")) droplevels(d[d$participant_id_internal!=sub("^delete_","",label),]) else d
    m <- readRDS(file.path(out,paste0(ex,"_",label,".rds")))
    opt <- attr(m,"robustness_optimizer")
    expected <- selection[selection$experiment==ex & selection$model==label,]
    stopifnot(nrow(expected)==1L,opt==expected$optimizer,abs(AIC(m)-expected$AIC)<1e-10)
    fitted <- glmmTMB(formula(m),data=dd,family=binomial(),REML=FALSE,
      control=if(opt=="BFGS") custom_control_glmmTMB else custom_control_glmmTMB_nlminb)
    delta <- c(beta=max(abs(fixef(m)$cond-fixef(fitted)$cond)),
      vcov=max(abs(vcov(m)$cond-vcov(fitted)$cond)),logLik=abs(as.numeric(logLik(m)-logLik(fitted))))
    comparison[[paste(ex,label)]] <- data.frame(experiment=ex,model=label,optimizer=opt,
      metric=names(delta),maximum_absolute_difference=unname(delta),valid=isTRUE(mixed_model_diagnostics(fitted)$valid))
    readr::write_csv(do.call(rbind,comparison),file.path(out,"independent_refit_comparison.csv"))
    stopifnot(all(delta==0),isTRUE(mixed_model_diagnostics(fitted)$valid))
    rm(m,fitted,dd);gc(verbose=FALSE)
  }
}
writeLines(c("VERIFIED: exact deterministic numerical reproduction for the valid actor model (if any) and the largest standardized-change deletion per experiment.",
  "Canonical input reconstruction; identical beta, covariance and log-likelihood. Full deletion sequence and simulations were not rerun.",
  "This checks computation, not all model assumptions, theory identification or human verification.",capture.output(sessionInfo())),dest)
