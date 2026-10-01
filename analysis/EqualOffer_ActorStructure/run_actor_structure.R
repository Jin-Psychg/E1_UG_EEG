#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Supply primary run, first robustness run and NEW output directory")
run <- normalizePath(args[1],winslash="/",mustWork=TRUE)
previous <- normalizePath(args[2],winslash="/",mustWork=TRUE)
out <- normalizePath(args[3],winslash="/",mustWork=FALSE)
if(dir.exists(out)||file.exists(out)) stop("Preserve existing output")
root <- run
while(!file.exists(file.path(root,"UG_ERP_Project.Rproj"))) {
  parent <- dirname(root);if(parent==root) stop("Project root not found");root <- parent
}
if(!startsWith(out,paste0(root,"/results/"))) stop("Output must be under project results")
.libPaths(c(file.path(root,"renv/library/R-4.3/x86_64-w64-mingw32"),.libPaths()))
required <- c("lme4","lmerTest","glmmTMB","buildmer","dplyr","stringr","emmeans","car","R.utils","readr","jsonlite","digest")
if(any(!vapply(required,requireNamespace,logical(1),quietly=TRUE))) stop("Existing dependencies missing; do not auto-install")
suppressPackageStartupMessages(pacman::p_load(char=required,install=FALSE))
options(contrasts=c("contr.sum","contr.poly"),digits=17)
code <- file.path(root,"analysis/EqualOffer_ActorStructure")
provenance <- jsonlite::read_json(file.path(code,"engine_provenance.json"))
stopifnot(digest::digest(file=provenance$source,algo="sha256")==provenance$source_sha256,
          digest::digest(file=file.path(code,"actor_engine.R"),algo="sha256")==provenance$adapted_sha256)
source(file.path(code,"actor_engine.R"))
custom_control_glmmTMB <- glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),optCtrl=list(maxit=2000),parallel=1)
custom_control_glmmTMB_nlminb <- glmmTMBControl(optimizer=nlminb,optCtrl=list(iter.max=2000,eval.max=4000),parallel=1)
helper <- file.path(previous,"delivery_code/run_robustness.R")
helper_manifest <- jsonlite::read_json(file.path(previous,"delivery_code/code_manifest_robustness.json"))
stopifnot(digest::digest(file=helper,algo="sha256")==helper_manifest[["run_robustness.R"]])
inputs <- read.csv(file.path(previous,"input_manifest.csv"))
stopifnot(identical(unname(tools::md5sum(inputs$path)),inputs$md5))
dir.create(out);dir.create(file.path(out,"code"))
file.copy(list.files(code,full.names=TRUE),file.path(out,"code"))
file.copy(helper,file.path(out,"code/run_robustness_helpers_source.R"))
save_csv <- function(x,name) readr::write_csv(x,file.path(out,name),na="NA")
writeLines("RUNNING",file.path(out,"STATUS.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
paths <- unique(c(inputs$path,helper,file.path(code,"actor_engine.R")))
hashes <- data.frame(path=paths,md5=unname(tools::md5sum(paths)))
save_csv(hashes,"input_manifest.csv")
# Reuse validated fit/contrast helpers without executing their original pipeline.
needed <- c("check","covariance_rows","fit_exact","extract");loaded <- character()
for(expr in parse(helper)) if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&
  is.symbol(expr[[2]])&&as.character(expr[[2]]) %in% needed) {
    eval(expr,envir=.GlobalEnv);loaded <- c(loaded,as.character(expr[[2]]))
}
stopifnot(setequal(needed,loaded))
checks <- covariances <- attempts <- effects <- omnibus <- selection_status <- structure_audit <- list()
formula_text <- function(f) paste(deparse(f),collapse=" ")
bars <- function(f) vapply(lme4::findbars(f),function(b)paste(deparse(b),collapse=" "),character(1))
for(ex in c("E1","E2")) {
  message("Actor-preserving selection: ",ex)
  d <- read.csv(file.path(run,paste0(ex,"_analysis_trials.csv")))
  d$emotion <- factor(d$emotion,levels=c("neu","aff","dis","dom","enj"))
  d$allocation <- factor(d$allocation,levels=c("5:5","6:4"))
  d$participant_id_internal <- factor(d$participant_id_internal);d$actor_id <- factor(d$actor_id)
  check(paste(ex,"data_integrity"),!anyNA(d)&&all(d$reject_binary==as.integer(d$reaction==2))&&
    all(d$RT>=300&d$RT<=3000)&&!anyDuplicated(paste(d$participant_id,d$index)))
  d <- add_random_contrast_columns(d,c("emotion","allocation"),"participant_id_internal")$data
  primary <- readRDS(file.path(run,paste0(ex,"_fair_choice.rds")))
  frame <- model.frame(primary)
  check(paste(ex,"primary_frame"),all(vapply(names(frame),function(n)isTRUE(all.equal(frame[[n]],d[[n]],check.attributes=FALSE)),logical(1))))
  extract(primary,ex,"original_primary")
  stage_dir <- file.path(out,paste0(ex,"_selection_stages"));dir.create(stage_dir)
  log <- file(file.path(out,paste0(ex,"_selection.log")),"wt");sink(log,split=TRUE)
  selected <- tryCatch(fit_stage1_actor(d,"reject_binary","emotion * allocation",
    c("emotion","allocation","emotion:allocation"),model_family="binomial",
    random_max_vars=c("emotion","allocation"),timeout_sec=1800,actor_audit_dir=stage_dir),error=identity)
  sink();close(log)
  if(inherits(selected,"error")) {
    selection_status[[ex]] <- data.frame(experiment=ex,actor_valid=FALSE,matched_pair_valid=FALSE,
      reason=conditionMessage(selected),actor_formula=NA_character_,no_actor_formula=NA_character_)
    save_csv(dplyr::bind_rows(selection_status),"selection_status.csv")
    message(ex," stopped after selection failure: ",conditionMessage(selected))
    next
  }
  check(paste(ex,"selected_actor_gate"),isTRUE(mixed_model_diagnostics(selected)$valid))
  fa <- formula(selected); participant_bars <- bars(fa)[grepl("participant_id_internal",bars(fa),fixed=TRUE)]
  actor_bars <- bars(fa)[grepl("actor_id",bars(fa),fixed=TRUE)]
  check(paste(ex,"actor_retained"),length(actor_bars)==1L && actor_bars=="1 | actor_id")
  check(paste(ex,"all_fixed_design_terms"),setequal(attr(terms(lme4::nobars(fa)),"term.labels"),c("emotion","allocation","emotion:allocation")))
  present <- unique(unlist(lapply(lme4::findbars(fa),function(b) all.vars(b[[2]]))))
  blocks <- list(emotion=paste0("RE_emotion_",1:4),allocation="RE_allocation_1",interaction=paste0("RE_emotion_x_allocation_",1:4))
  for(n in names(blocks)) check(paste(ex,"whole_block",n),sum(blocks[[n]] %in% present) %in% c(0,length(blocks[[n]])))
  saveRDS(selected,file.path(out,paste0(ex,"_actor_selected.rds")))
  saveRDS(attr(selected,"bates_diag"),file.path(out,paste0(ex,"_selection_diagnostics.rds")))
  jsonlite::write_json(attr(selected,"bates_diag"),file.path(out,paste0(ex,"_selection_diagnostics.json")),auto_unbox=TRUE,pretty=TRUE,na="null")
  covariances[[length(covariances)+1L]] <- covariance_rows(selected,ex,"actor_selected",attr(selected,"final_refit_optimizer"))
  save_csv(dplyr::bind_rows(covariances),"covariance_diagnostics.csv")
  fn <- update(fa,. ~ . - (1|actor_id))
  check(paste(ex,"matched_random_structure"),setequal(participant_bars,bars(fn))&&length(bars(fn))==length(participant_bars))
  check(paste(ex,"matched_fixed_structure"),identical(colnames(model.matrix(lme4::nobars(fa),d)),colnames(model.matrix(lme4::nobars(fn),d))))
  matched <- fit_exact(fn,d,ex,"matched_no_actor")
  pair_valid <- !is.null(matched)
  selection_status[[ex]] <- data.frame(experiment=ex,actor_valid=TRUE,matched_pair_valid=pair_valid,
    reason=if(pair_valid) "accepted_matched_pair" else "matched_no_actor_failed_original_gate",
    actor_formula=formula_text(fa),no_actor_formula=formula_text(fn))
  save_csv(dplyr::bind_rows(selection_status),"selection_status.csv")
  if(!pair_valid) next
  extract(selected,ex,"actor_selected");extract(matched,ex,"matched_no_actor")
  ids <- levels(d$participant_id_internal)
  save_csv(data.frame(participant_id_internal=ids),paste0(ex,"_deletion_roster.csv"))
  for(id in ids) {
    label <- paste0("actor_delete_",id);message(ex," ",label)
    dd <- droplevels(d[d$participant_id_internal!=id,])
    m <- fit_exact(fa,dd,ex,label)
    if(!is.null(m)) {
      check(paste(ex,label,"fixed_selected_formula"),formula_text(formula(m))==formula_text(fa)&&nobs(m)==nrow(dd))
      extract(m,ex,label)
    }
    rm(m,dd);gc(verbose=FALSE)
  }
  rm(primary,selected,matched,d);gc(verbose=FALSE)
}
check("source_hashes_unchanged",identical(unname(tools::md5sum(hashes$path)),hashes$md5))
writeLines(c("Retrospective exploratory actor-intercept sensitivity. Fixed emotion * allocation retained.",
  "Original grouped random-slope LRT alpha=.20 and correlation comparison rules; actor intercept forced throughout.",
  "Extended optimizer budgets recorded in ANALYSIS_PLAN.md. Marginal emotion pairs use BH; single allocation comparison; unadjusted Wald CIs.",
  "Original primary models/reporting routes are unchanged. Failed selection blocks dependent stages; no fallback.",
  "Actor intercepts address baseline identity variation only. Valid deletion fits do not cover failed ones."),file.path(out,"METHODS.md"))
writeLines("COMPLETE: inspect selection_status.csv for scientific coverage; independent verification pending",file.path(out,"STATUS.txt"))
message("Actor structure stage complete: ",out)
