#!/usr/bin/env Rscript
# Rscript --vanilla analysis/EqualOffer_Conditional/run_analysis.R results/EqualOffer_Conditional_NEW
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Supply one NEW output directory')
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
code <- normalizePath(dirname(self),winslash='/',mustWork=TRUE)
root <- normalizePath(file.path(code,'../..'),winslash='/',mustWork=TRUE)
stopifnot(file.exists(file.path(root,'UG_ERP_Project.Rproj')))
out <- normalizePath(args[1],winslash='/',mustWork=FALSE)
if(dir.exists(out)||file.exists(out)) stop('Preserve existing output: choose a new directory')
if(!startsWith(out,paste0(root,'/results/'))) stop('Output must be in project results')
.libPaths(c(file.path(root,'renv/library/R-4.3/x86_64-w64-mingw32'),.libPaths()))
required <- c('glmmTMB','lme4','emmeans','jsonlite','digest','readr')
if(any(!vapply(required,requireNamespace,logical(1),quietly=TRUE))) stop('Restore existing dependencies; no automatic install')
suppressPackageStartupMessages(pacman::p_load(char=required,install=FALSE))
options(contrasts=c('contr.sum','contr.poly'),digits=17)
primary <- file.path(root,'results/EqualOffer_Aligned_20260928_v2')
prior <- file.path(root,'results/EqualOffer_Robustness_20260928_v1')
retry <- file.path(root,'results/EqualOffer_NumericRetry_20260929_v3')
actor <- file.path(root,'results/EqualOffer_ActorStructure_20260929_v1')
engine <- file.path(root,'analysis/EqualOffer_ManuscriptAligned/manuscript_engine.R')
ep <- jsonlite::read_json(file.path(dirname(engine),'engine_provenance.json'))
stopifnot(digest::digest(file=engine,algo='sha256')==ep$snapshot_sha256)
source(engine)
dir.create(out);dir.create(file.path(out,'code'))
file.copy(list.files(code,full.names=TRUE),file.path(out,'code'))
writeLines('RUNNING',file.path(out,'STATUS.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'SESSION_INFO.txt'))
save_csv <- function(x,n) readr::write_csv(x,file.path(out,n),na='NA')
roster <- data_rows <- list()
for(ex in c('E1','E2')) {
 d <- read.csv(file.path(primary,paste0(ex,'_analysis_trials.csv')))
 stopifnot(all(c('participant_id','participant_id_internal','emotion','allocation','reject_binary','index','reaction','RT')%in%names(d)),!anyNA(d),all(d$reject_binary==as.integer(d$reaction==2)),all(d$RT>=300&d$RT<=3000),!anyDuplicated(paste(d$participant_id,d$index)))
 data_rows[[ex]] <- d
 add <- function(label,path,id='') roster[[length(roster)+1L]] <<- data.frame(experiment=ex,model=label,deleted_id=id,path=path,exists=file.exists(path))
 add('primary',file.path(primary,paste0(ex,'_fair_choice.rds')))
 add('actor_selected',file.path(actor,paste0(ex,'_actor_selected.rds')))
 add('matched_no_actor',file.path(actor,paste0(ex,'_matched_no_actor.rds')))
 for(id in unique(d$participant_id_internal)) {
  old <- file.path(prior,paste0(ex,'_delete_',id,'.rds'))
  new <- file.path(retry,paste0(ex,'_delete_',id,'.rds'))
  if(file.exists(old)&&file.exists(new)) stop('Ambiguous accepted deletion source')
  add('primary_delete',if(file.exists(new)) new else old,id)
  add('actor_delete',file.path(actor,paste0(ex,'_actor_delete_',id,'.rds')),id)
 }
}
roster <- do.call(rbind,roster);save_csv(roster,'model_roster.csv')
paths <- unique(c(roster$path[roster$exists],engine,file.path(code,c('ANALYSIS_PLAN.md','run_analysis.R')),file.path(primary,c('E1_analysis_trials.csv','E2_analysis_trials.csv','descriptive_cells.csv','participant_cells.csv')),file.path(retry,'retry_status.csv'),file.path(actor,'deletion_coverage.csv')))
manifest <- data.frame(path=paths,sha256=vapply(paths,function(p)digest::digest(file=p,algo='sha256'),character(1)))
save_csv(manifest,'input_manifest.csv')
results <- estimates <- diagnostics <- verification <- list()
for(i in seq_len(nrow(roster))) {
 row <- roster[i,];if(!row$exists) next
 m <- readRDS(row$path);dg <- mixed_model_diagnostics(m)
 diagnostics[[length(diagnostics)+1L]] <- data.frame(experiment=row$experiment,model=row$model,deleted_id=row$deleted_id,path=row$path,valid=dg$valid,reason=format_model_diagnostics(dg),formula=paste(deparse(formula(m)),collapse=' '))
 if(!isTRUE(dg$valid)) next
 d <- data_rows[[row$experiment]]
 if(nchar(row$deleted_id)) d <- d[d$participant_id_internal!=row$deleted_id,]
 mf <- model.frame(m)
 stopifnot(nrow(mf)==nrow(d),identical(as.integer(mf$reject_binary),as.integer(d$reject_binary)),identical(as.character(mf$emotion),as.character(d$emotion)),identical(as.character(mf$allocation),as.character(d$allocation)),identical(as.character(mf$participant_id_internal),as.character(d$participant_id_internal)))
 stopifnot(setequal(attr(terms(lme4::nobars(formula(m))),'term.labels'),c('emotion','allocation','emotion:allocation')))
 eg <- emmeans::emmeans(m,~emotion,at=list(allocation='5:5'),data=mf)
 tab <- as.data.frame(summary(eg,type='response',infer=c(TRUE,FALSE)))
 stopifnot(identical(as.character(tab$emotion),c('neu','aff','dis','dom','enj')))
 weights <- lapply(2:5,function(j){x<-numeric(5);x[1]<- -1;x[j]<-1;x})
 names(weights) <- paste(c('aff','dis','dom','enj'),'minus neu')
 ct <- emmeans::contrast(eg,method=weights)
 joint <- as.data.frame(emmeans::test(ct,joint=TRUE))
 stopifnot(joint$df1==4,is.infinite(joint$df2))
 b <- fixef(m)$cond;V <- vcov(m)$cond
 grid <- data.frame(emotion=factor(c('neu','aff','dis','dom','enj'),levels=c('neu','aff','dis','dom','enj')),allocation=factor(rep('5:5',5),levels=c('5:5','6:4')))
 X <- model.matrix(~emotion*allocation,grid)[,names(b),drop=FALSE]
 L <- X[2:5,,drop=FALSE]-matrix(X[1,],nrow=4,ncol=ncol(X),byrow=TRUE)
 theta <- as.vector(L%*%b);S <- L%*%V%*%t(L)
 ev <- eigen(S,symmetric=TRUE,only.values=TRUE)$values
 stopifnot(all(is.finite(ev)),min(ev)>0)
 W <- as.numeric(crossprod(theta,solve(S,theta)))
 pval <- pchisq(W,df=4,lower.tail=FALSE)
 eta <- as.vector(X%*%b);se <- sqrt(diag(X%*%V%*%t(X)))
 stopifnot(abs(round(W,3)-joint$Chisq)<1e-10,abs(pval-joint$p.value)<1e-10,max(abs(plogis(eta)-tab$prob))<1e-10)
 verification[[length(verification)+1L]] <- data.frame(experiment=row$experiment,model=row$model,deleted_id=row$deleted_id,frame_matches=TRUE,joint_direct_matches=TRUE,probability_direct_matches=TRUE)
 results[[length(results)+1L]] <- data.frame(experiment=row$experiment,model=row$model,deleted_id=row$deleted_id,Chisq=W,df=4,p_raw=pval,contrast_cov_condition=max(ev)/min(ev),n_participants=length(unique(d$participant_id_internal)),n_trials=nrow(d),n_equal_trials=sum(d$allocation=='5:5'),equal_rejections=sum(d$reject_binary[d$allocation=='5:5']),source=row$path)
 tab$experiment <- row$experiment;tab$model <- row$model;tab$deleted_id <- row$deleted_id;tab$logit <- eta;tab$logit_SE <- se
 estimates[[length(estimates)+1L]] <- tab
 if(i%%20==0) message('Processed ',i,'/',nrow(roster))
 rm(m);gc(verbose=FALSE)
}
r <- do.call(rbind,results);r$p_BH_two_experiments <- NA_real_
for(label in c('primary','actor_selected','matched_no_actor')) {
 k <- which(r$model==label);stopifnot(length(k)==2)
 r$p_BH_two_experiments[k] <- p.adjust(r$p_raw[k],method='BH')
}
save_csv(r,'conditional_omnibus.csv');save_csv(do.call(rbind,estimates),'conditional_probabilities.csv')
save_csv(do.call(rbind,diagnostics),'diagnostics.csv');save_csv(do.call(rbind,verification),'verification_checks.csv')
stopifnot(identical(vapply(manifest$path,function(p)digest::digest(file=p,algo='sha256'),character(1)),setNames(manifest$sha256,manifest$path)))
for(ex in c('E1','E2')) writeLines(c('Exploratory 5:5 conditional expression omnibus from the existing full emotion * allocation Bernoulli-logit GLMM.',
'Four independent expression contrasts at 5:5 jointly tested by asymptotic Wald chi-square; BH across the two primary experiment tests.',
'Participant structure and fixed interaction unchanged. Conditional probabilities set random effects to zero; pointwise 95% Wald CIs.',
'No pairwise significance tests. Existing actor/matched and accepted deletion fits provide sensitivity estimates; failed fits are retained in the roster as missing coverage.',
'No refits or new model selection. Sparse-data finite-sample calibration and model-selection uncertainty remain unresolved.'),file.path(out,paste0('Methods_paragraph_',ex,'.md')))
writeLines(c('See input_manifest.csv for exact source identities; all hashes checked before and after extraction.',
'No raw inputs, historical models or historical outputs modified. Independently verified conditional estimates and joint tests.',
'Primary and sensitivity results must remain separate. Missing deletion models are not zero effects.'),file.path(out,'REPRODUCIBILITY.txt'))
writeLines('COMPLETE: independent summary verification pending',file.path(out,'STATUS.txt'))
print(r[r$model%in%c('primary','actor_selected','matched_no_actor'),c('experiment','model','Chisq','p_raw','p_BH_two_experiments')])


