#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Supply primary run, prior robustness run and NEW output")
run <- normalizePath(args[1],winslash="/",mustWork=TRUE)
previous <- normalizePath(args[2],winslash="/",mustWork=TRUE)
out <- normalizePath(args[3],winslash="/",mustWork=FALSE)
if(file.exists(out)||dir.exists(out)) stop("Preserve existing output; choose a new destination")
root <- run
while(!file.exists(file.path(root,"UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if(parent==root) stop("Project root not found"); root <- parent
}
if(!startsWith(out,paste0(root,"/results/"))) stop("Output must be under project results")
.libPaths(c(file.path(root,"renv/library/R-4.3/x86_64-w64-mingw32"),.libPaths()))
suppressPackageStartupMessages({library(lme4);library(glmmTMB);library(emmeans)})
options(contrasts=c("contr.sum","contr.poly"),digits=17)
helper <- file.path(previous,"delivery_code/run_robustness.R")
manifest <- jsonlite::read_json(file.path(previous,"delivery_code/code_manifest_robustness.json"))
stopifnot(digest::digest(file=helper,algo="sha256")==manifest[["run_robustness.R"]])
engine <- file.path(previous,"delivery_code/manuscript_engine.R")
stopifnot(digest::digest(file=engine,algo="sha256")==manifest[["manuscript_engine.R"]])
source(engine)
old_inputs <- read.csv(file.path(previous,"input_manifest.csv"))
stopifnot(identical(unname(tools::md5sum(old_inputs$path)),old_inputs$md5))
old_attempts <- read.csv(file.path(previous,"fit_attempts.csv"))
coverage <- aggregate(valid~experiment+model,old_attempts,any)
targets <- coverage[!coverage$valid,c("experiment","model")]
stopifnot(nrow(targets)>0,!anyDuplicated(targets))
dir.create(out);dir.create(file.path(out,"code"));dir.create(file.path(out,"attempt_models"))
code <- file.path(root,"analysis/EqualOffer_ManuscriptAligned")
file.copy(file.path(code,c("run_numeric_retry.R","NUMERIC_RETRY_PLAN.md")),file.path(out,"code"))
file.copy(c(helper,engine),file.path(out,"code"))
writeLines("RUNNING",file.path(out,"STATUS.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
save_csv <- function(x,name) readr::write_csv(x,file.path(out,name),na="NA")
save_csv(targets,"retry_targets.csv")
paths <- unique(c(old_inputs$path,helper,engine,file.path(previous,"fit_attempts.csv")))
hashes <- data.frame(path=paths,md5=unname(tools::md5sum(paths)))
save_csv(hashes,"input_manifest.csv")
# Reuse only these named function definitions; no top-level fit or I/O is evaluated.
names_needed <- c("check","covariance_rows","extract"); loaded <- character()
for(expr in parse(helper)) {
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&is.symbol(expr[[2]])&&
     as.character(expr[[2]]) %in% names_needed) {
    eval(expr,envir=.GlobalEnv);loaded <- c(loaded,as.character(expr[[2]]))
  }
}
stopifnot(setequal(loaded,names_needed))
checks <- effects <- omnibus <- attempts <- covariances <- comparisons <- selections <- formula_checks <- list()
controls <- list(BFGS=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),
  optCtrl=list(maxit=2000),parallel=1),
  nlminb=glmmTMBControl(optimizer=nlminb,optCtrl=list(iter.max=2000,eval.max=4000),parallel=1))
saveRDS(controls,file.path(out,"optimizer_controls.rds"))
writeLines(capture.output(lapply(controls,function(x)x[c("optCtrl","optArgs","parallel") ])),file.path(out,"optimizer_controls.txt"))
for(ex in unique(targets$experiment)) {
  d <- read.csv(file.path(run,paste0(ex,"_analysis_trials.csv")))
  d$emotion <- factor(d$emotion,levels=c("neu","aff","dis","dom","enj"))
  d$allocation <- factor(d$allocation,levels=c("5:5","6:4"))
  d$participant_id_internal <- factor(d$participant_id_internal);d$actor_id <- factor(d$actor_id)
  check(paste(ex,"input"),!anyNA(d)&&all(d$reject_binary==as.integer(d$reaction==2)))
  d <- add_random_contrast_columns(d,c("emotion","allocation"),"participant_id_internal")$data
  primary <- readRDS(file.path(run,paste0(ex,"_fair_choice.rds")))
  frame <- model.frame(primary)
  check(paste(ex,"primary_frame"),all(vapply(names(frame),function(n)isTRUE(all.equal(frame[[n]],d[[n]],check.attributes=FALSE)),logical(1))))
  extract(primary,ex,"primary")
  for(label in targets$model[targets$experiment==ex]) {
    message("Numerical retry: ",ex," ",label)
    dd <- if(startsWith(label,"delete_")) droplevels(d[d$participant_id_internal!=sub("^delete_","",label),]) else d
    f <- if(label=="actor_intercept") update(formula(primary),.~.+(1|actor_id)) else formula(primary)
    candidates <- list()
    for(opt in names(controls)) {
      warnings <- character();err <- "";start <- proc.time()[[3]]
      m <- tryCatch(withCallingHandlers(R.utils::withTimeout(
        glmmTMB(f,data=dd,family=binomial(),REML=FALSE,control=controls[[opt]]),
        timeout=300,onTimeout="error"),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
        error=function(e){err <<- conditionMessage(e);NULL})
      dg <- if(is.null(m)) NULL else mixed_model_diagnostics(m)
      valid <- !is.null(dg)&&isTRUE(dg$valid)
      attempts[[length(attempts)+1L]] <- data.frame(experiment=ex,model=label,optimizer=opt,
        n_trials=nrow(dd),n_participants=length(unique(dd$participant_id_internal)),valid=valid,
        AIC=if(is.null(m)) NA_real_ else AIC(m),objective=if(is.null(m)) NA_real_ else m$fit$objective,
        elapsed_seconds=proc.time()[[3]]-start,error=err,warnings=paste(unique(warnings),collapse=" | "),
        diagnostics=if(is.null(dg)) "fit_failed" else format_model_diagnostics(dg))
      save_csv(do.call(rbind,attempts),"fit_attempts.csv")
      if(!is.null(m)) {
        # Compare canonical formula text, excluding retained call/environment attributes.
        formula_checks[[length(formula_checks)+1L]] <- data.frame(experiment=ex,model=label,optimizer=opt,
          expected_formula=paste(deparse(f),collapse=" "),actual_formula=paste(deparse(formula(m)),collapse=" "),
          expression_equal=identical(paste(deparse(formula(m)),collapse=" "),paste(deparse(f),collapse=" ")),
          expected_rows=nrow(dd),actual_rows=nobs(m))
        save_csv(do.call(rbind,formula_checks),"formula_checks.csv")
        check(paste(ex,label,opt,"formula_and_rows"),isTRUE(tail(formula_checks,1)[[1]]$expression_equal)&&nobs(m)==nrow(dd))
        attr(m,"retry_optimizer") <- opt
        # Invalid objects are explicitly diagnostic-only; extract() is never called on them.
        saveRDS(m,file.path(out,"attempt_models",paste0(ex,"_",label,"_",opt,".rds")))
        covariances[[length(covariances)+1L]] <- covariance_rows(m,ex,label,opt)
        save_csv(do.call(rbind,covariances),"covariance_diagnostics.csv")
        candidates[[opt]] <- m
      }
    }
    if(length(candidates)==2L) {
      a <- candidates[[1]];b <- candidates[[2]]
      comparisons[[length(comparisons)+1L]] <- data.frame(experiment=ex,model=label,
        absolute_objective_difference=abs(a$fit$objective-b$fit$objective),
        maximum_fixed_beta_difference=max(abs(fixef(a)$cond-fixef(b)$cond)),
        both_pass_original_gate=isTRUE(mixed_model_diagnostics(a)$valid)&&isTRUE(mixed_model_diagnostics(b)$valid))
      save_csv(do.call(rbind,comparisons),"optimizer_comparison.csv")
    }
    ok <- vapply(candidates,function(m)isTRUE(mixed_model_diagnostics(m)$valid),logical(1))
    accepted <- candidates[ok]
    opt <- if(length(accepted)) names(which.min(vapply(accepted,AIC,numeric(1)))) else NA_character_
    selections[[length(selections)+1L]] <- data.frame(experiment=ex,model=label,
      recovered_original_gate=length(accepted)>0,selected_optimizer=opt,formula=paste(deparse(f),collapse=" "))
    save_csv(do.call(rbind,selections),"retry_status.csv")
    if(length(accepted)) {
      chosen <- accepted[[opt]]
      saveRDS(chosen,file.path(out,paste0(ex,"_",label,".rds")))
      extract(chosen,ex,label)
    }
    rm(candidates,accepted);gc(verbose=FALSE)
  }
}
check("all_sources_unchanged",identical(unname(tools::md5sum(hashes$path)),hashes$md5))
check("complete_attempt_coverage",length(attempts)==2L*nrow(targets))
writeLines("COMPLETE: bounded unchanged-formula retries; independent report verification pending",file.path(out,"STATUS.txt"))
message("Completed numerical retry stage: ",out)
