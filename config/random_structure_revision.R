# Fixed formulas and original model identities for the September 2026 random-slope revision.
# Extracted from original saved fits; no fitted estimates are embedded.
revision_models <- list(
  "Behavior__E1__GLMM_rejection" = list(formula = reject_binary ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "Behavior/E1_TwoStage_Bates/GLMM_rejection/final_model_GLMM_rejection.rds", md5 = "bc80c821cef3bfde140bd2e49d4645f4"),
  "Behavior__E1__LMM_RT_main" = list(formula = logRT ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "Behavior/E1_TwoStage_Bates/LMM_RT_main/final_model_LMM_RT_main.rds", md5 = "ab119570c85a09e8b9c83f9a24290661"),
  "Behavior__E1__LMM_RT_unfair" = list(formula = logRT ~ 1 + emotion + reaction + emotion:reaction,
    path = "Behavior/E1_TwoStage_Bates/LMM_RT_unfair/final_model_LMM_RT_unfair.rds", md5 = "55570ee48c46ac52fcd10b5036548ef6"),
  "Behavior__E2__GLMM_rejection" = list(formula = reject_binary ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "Behavior/E2_TwoStage_Bates/GLMM_rejection/final_model_GLMM_rejection.rds", md5 = "8fb80b53ff9928cd2b6a5307a4edafe2"),
  "Behavior__E2__LMM_RT_main" = list(formula = logRT ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "Behavior/E2_TwoStage_Bates/LMM_RT_main/final_model_LMM_RT_main.rds", md5 = "8a59c87e8cb9e2d25e93b775224d8582"),
  "Behavior__E2__LMM_RT_unfair" = list(formula = logRT ~ 1 + emotion + reaction + emotion:reaction,
    path = "Behavior/E2_TwoStage_Bates/LMM_RT_unfair/final_model_LMM_RT_unfair.rds", md5 = "65bfc5a8bea064d8b4aff0129d086acb"),
  "Behavior__Integrative__GLMM_rejection" = list(formula = reject_binary ~ 1 + Exp + emotion + offer_type + Exp:emotion + Exp:offer_type + emotion:offer_type + Exp:emotion:offer_type,
    path = "Behavior/Integrative_TwoStage_Bates/GLMM_rejection/final_model_GLMM_rejection.rds", md5 = "832666ce86c5eaf2a0be0acda85bbf15"),
  "Behavior__Integrative__LMM_RT_main" = list(formula = logRT ~ 1 + Exp + emotion + offer_type + Exp:emotion + Exp:offer_type + emotion:offer_type + Exp:emotion:offer_type,
    path = "Behavior/Integrative_TwoStage_Bates/LMM_RT_main/final_model_LMM_RT_main.rds", md5 = "037b2ed548d415f9bc68f7fbc71eb83d"),
  "Behavior__Integrative__LMM_RT_unfair" = list(formula = logRT ~ 1 + Exp + emotion + reaction + Exp:emotion + Exp:reaction + emotion:reaction + Exp:emotion:reaction,
    path = "Behavior/Integrative_TwoStage_Bates/LMM_RT_unfair/final_model_LMM_RT_unfair.rds", md5 = "4e5ee820679a32dc3f67fb9b327b8036"),
  "Alday__E1__FRN_pre" = list(formula = FRN_pre ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + offer_type:Baseline_c + emotion:Baseline_c,
    path = "EEG/E1_TwoStage_Bates_Alday/Stage1_TrialLevel/FRN_pre/final_model_FRN_pre.rds", md5 = "00cc0549896a9d89a1fde11e0f8054b3"),
  "Alday__E1__LPP_pre" = list(formula = LPP_pre ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + offer_type:Baseline_c + emotion:Baseline_c,
    path = "EEG/E1_TwoStage_Bates_Alday/Stage1_TrialLevel/LPP_pre/final_model_LPP_pre.rds", md5 = "75b667551329677b2e5ebd0b8542b273"),
  "Alday__E1__FRN_explor" = list(formula = FRN_explor ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + offer_type:Baseline_c + emotion:Baseline_c,
    path = "EEG/E1_TwoStage_Bates_Alday/Stage1_TrialLevel/FRN_explor/final_model_FRN_explor.rds", md5 = "ae33a650c66d73f2bec8fd9c3b17dbdb"),
  "Alday__E1__P3_explor" = list(formula = P3_explor ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + offer_type:Baseline_c + emotion:Baseline_c,
    path = "EEG/E1_TwoStage_Bates_Alday/Stage1_TrialLevel/P3_explor/final_model_P3_explor.rds", md5 = "9f1bcee1d057e6df19e9438f2d49a6a6"),
  "Alday__E2__FRN_pre" = list(formula = FRN_pre ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type,
    path = "EEG/E2_TwoStage_Bates_Alday/Stage1_TrialLevel/FRN_pre/final_model_FRN_pre.rds", md5 = "762ec77db7d52bf78f7d8dad05846119"),
  "Alday__E2__LPP_pre" = list(formula = LPP_pre ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + emotion:Baseline_c + offer_type:Baseline_c,
    path = "EEG/E2_TwoStage_Bates_Alday/Stage1_TrialLevel/LPP_pre/final_model_LPP_pre.rds", md5 = "b0771c52054b77171140e102e62b6def"),
  "Alday__E2__FRN_explor" = list(formula = FRN_explor ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + offer_type:Baseline_c + emotion:Baseline_c,
    path = "EEG/E2_TwoStage_Bates_Alday/Stage1_TrialLevel/FRN_explor/final_model_FRN_explor.rds", md5 = "c02c6f86877a6cb66e9d9ef150f0170f"),
  "Alday__E2__P3_explor" = list(formula = P3_explor ~ 1 + emotion + offer_type + Baseline_c + emotion:offer_type + emotion:Baseline_c + offer_type:Baseline_c,
    path = "EEG/E2_TwoStage_Bates_Alday/Stage1_TrialLevel/P3_explor/final_model_P3_explor.rds", md5 = "70214e73eedc2f5006989a2aa4e8e0ce"),
  "Traditional__E1__FRN_pre" = list(formula = FRN_pre ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E1_TwoStage_Bates_Traditional/Stage1_TrialLevel/FRN_pre/final_model_FRN_pre.rds", md5 = "3eb7da9d8796d0602ba445ce4499f7bf"),
  "Traditional__E1__LPP_pre" = list(formula = LPP_pre ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E1_TwoStage_Bates_Traditional/Stage1_TrialLevel/LPP_pre/final_model_LPP_pre.rds", md5 = "efa03b27c6a63cc354d9a0056d64dec8"),
  "Traditional__E1__FRN_explor" = list(formula = FRN_explor ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E1_TwoStage_Bates_Traditional/Stage1_TrialLevel/FRN_explor/final_model_FRN_explor.rds", md5 = "bc224119d503b1829a235f34729c5e9c"),
  "Traditional__E1__P3_explor" = list(formula = P3_explor ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E1_TwoStage_Bates_Traditional/Stage1_TrialLevel/P3_explor/final_model_P3_explor.rds", md5 = "bbe2e3b6ca309bd00ec233419496c530"),
  "Traditional__E2__FRN_pre" = list(formula = FRN_pre ~ emotion * offer_type,
    path = "EEG/E2_TwoStage_Bates_Traditional/Stage1_TrialLevel/FRN_pre/final_model_FRN_pre.rds", md5 = "3494f49bf4cbe02a5847d0010cab5b69"),
  "Traditional__E2__LPP_pre" = list(formula = LPP_pre ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E2_TwoStage_Bates_Traditional/Stage1_TrialLevel/LPP_pre/final_model_LPP_pre.rds", md5 = "7cbcf454711f094319ee48f30023a339"),
  "Traditional__E2__FRN_explor" = list(formula = FRN_explor ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E2_TwoStage_Bates_Traditional/Stage1_TrialLevel/FRN_explor/final_model_FRN_explor.rds", md5 = "5bca932e28b7cc1aef94a1d44624db4d"),
  "Traditional__E2__P3_explor" = list(formula = P3_explor ~ 1 + emotion + offer_type + emotion:offer_type,
    path = "EEG/E2_TwoStage_Bates_Traditional/Stage1_TrialLevel/P3_explor/final_model_P3_explor.rds", md5 = "ce54ce2bc21f604b39c816b7b2f9843c")
)

# Validate loaded input without replacing it with saved model data.
# Enable ug.verify_reference for the controlled original-versus-revised check.
revision_reference <- function(data, model_id, output_dir) {
  ref <- revision_models[[model_id]]
  if (is.null(ref)) stop("No reference specification: ", model_id)
  status <- "Fixed formula preserved; original-data comparison not requested."
  if (isTRUE(getOption("ug.verify_reference", FALSE))) {
    path <- here::here("results", ref$path)
    if (!file.exists(path) || unname(tools::md5sum(path)) != ref$md5)
      stop("Original reference model missing or changed: ", path)
    old <- readRDS(path)
    frame <- model.frame(old)
    columns <- names(frame)[!startsWith(names(frame), "RE_")]
    if (nrow(data) != nrow(frame) || !all(columns %in% names(data)))
      stop("Analysis rows or columns differ: ", model_id)
    group_name <- if (startsWith(model_id, "Behavior__")) "participant_id_internal" else "participant_id"
    for (v in columns) {
      # lme4 stores character grouping IDs as factors without changing membership.
      if (v == group_name) {
        if (!identical(as.character(data[[v]]), as.character(frame[[v]])))
          stop("Participant identities or order differ: ", model_id)
        next
      }
      if (!isTRUE(all.equal(data[[v]], frame[[v]], tolerance = 0, check.attributes = FALSE)) ||
          !identical(levels(data[[v]]), levels(frame[[v]])))
        stop("Loaded data differ from original fit: ", model_id, " / ", v)
    }
    expected <- if (inherits(old, "glmmTMB")) model.matrix(old, component = "cond") else lme4::getME(old, "X")
    actual <- model.matrix(ref$formula, data)
    if (!identical(colnames(actual), colnames(expected)) ||
        !isTRUE(all.equal(unname(actual), unname(expected), tolerance = 0, check.attributes = FALSE)))
      stop("Fixed design matrix differs: ", model_id)
    status <- "PASS: ordered analysis values, factor levels and fixed design matrix exactly match original fit."
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(model_id, paste(deparse(ref$formula), collapse = " "), status),
             file.path(output_dir, "INPUT_FIXED_VERIFICATION.txt"))
  ref$formula
}
