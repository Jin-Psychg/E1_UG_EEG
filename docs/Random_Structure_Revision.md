# Random structure revision and reproduction

The three `Sta_*.Rmd` statistics scripts now include the interaction of their within-participant experimental factors among the participant random-slope candidates. Fixed effects, eligibility rules, contrast coding, baseline data, optimizer policies, selection thresholds and final acceptance checks are retained.

## Run the actual Rmd

Open the project directory in PowerShell, using the recorded R 4.3.3 environment:

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla Run_Statistical_Revision.R Behavior E1 all preflight
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla Run_Statistical_Revision.R Behavior E1 all fit
```

Arguments are method (`Behavior`, `Alday`, `Traditional`), experiment (`E1`, `E2`; `Integrative` for behavior), target (`all`, a behavioral analysis name, or an ERP component name), and mode (`preflight`, `fit`). For example, `Alday E2 FRN_pre fit` runs that component through its formal Rmd. `Behavior E1 LMM_RT_unfair fit` runs the unfair-only RT model.

This entrypoint calls `rmarkdown::render()` on the project Rmd. It does not substitute a saved fit or replace the Rmd's loaders, fitting functions, or reporting functions. It requires the recorded software versions and stops if dependencies are missing instead of installing them. Inputs are existing preprocessing tables; it does not run preprocessing or HDDM. The bounded revision run omits trait screens, rating controls and already-completed baseline balance diagnostics. Direct Knit retains their existing defaults.

Every invocation writes to a new `results/RandomStructure_<...>/` directory. Inspect `RUN_RECEIPT.txt`, the HTML in `Reports/`, and model-level `INPUT_FIXED_VERIFICATION.txt`, `REPRODUCIBILITY.txt`, `final_model_*.rds`, and `STATS_*.xlsx`. Older result trees remain inputs and are not overwritten. Downstream figures must use the chosen new output folder explicitly; their old default paths do not automatically identify this revision.

## Fixed formulas and input checks

`config/random_structure_revision.R` records the fixed formulas and original-model identities extracted for this revision. In particular, EEG baseline interactions remain as selected in the original fits; they are not reselected when evaluating the enlarged random candidate set. This configuration is study-specific, not a rule for a future dataset.

The entrypoint compares freshly loaded analysis values, participant identities in row order, condition-factor levels, and the fixed design matrix against the original models. lme4's conversion of character grouping IDs to factors is distinguished from a change in participant membership. No values are replaced during validation. Original saved RDS files are required for this exact historical comparison and are not included in the public source repository. Direct Knit uses the fixed-formula configuration but does not require this historical comparison unless `options(ug.verify_reference = TRUE)` is set.

R 4.3.3 and lme4 1.1-38, lmerTest 3.1-3, glmmTMB 1.1.14, buildmer 2.12, emmeans 2.0.0 and Matrix 1.6-5 are checked. Full session information is included in the Rmd outputs. Check actual versions after restoring renv; the existing lockfile is not silently rewritten by this revision.

## Validation status

All 25 model-specific input/fixed-design checks passed through the formal scripts on 2026-09-15. The formal E1 choice fit reproduced the earlier controlled refit exactly in its final formula, fixed estimates, fixed covariance matrix and log likelihood, and passed the Hessian check. Full remaining fit/output validation is in progress; input checks alone do not establish full numerical reproduction.

The pooled choice model has a documented non-positive-definite Hessian in the Step 4 ZCP reference under the existing procedure. No optimizer-policy extension is included here. A failed fit must not be reported as an accepted result.

The Traditional script also contains a previously local, uncommitted dependency: when buildmer 2.12 reports that no removable terms remain, it fits the predefined participant-intercept endpoint and applies the unchanged final gate. This dependency and the existing session-log writing correction are included so the repository reproduces the local code used for the comparison; they were not invented in this revision.
