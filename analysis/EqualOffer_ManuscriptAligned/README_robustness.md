# Reproduce the fair-only robustness checks

These scripts supplement the selected-model analysis in
`results/EqualOffer_Aligned_20260928_v2`; they do not alter that run or its report.
Read ROBUSTNESS_PLAN.md before interpreting or rerunning. All checks are exploratory.

## Inputs and environment

- The complete local UG_ERP_Project tree, including its Rproj anchor and original
  `Sta_Behaviour_E1_E2_Integrative.Rmd` matching engine_provenance.json.
- The completed aligned run, including saved models, trial exports, primary
  contrast tables and STATUS.txt; no private participant data are in the ZIP.
- Canonical E1/E2 `data/02_Pipeline_Output_*/Method_Regression/Stimulus_Locked/trials.csv`.
- R 4.3.3 and the project's existing renv library. R dependencies: lme4, glmmTMB,
  emmeans, car, R.utils, readr, jsonlite and digest. The engine file retains some
  unused original selection helpers; this runner does not execute them.
- Python with numpy, pandas and scipy. Exact local versions are recorded in each
  completed run's provenance files. No script installs packages.

Place the scripts alongside manuscript_engine.R in
`analysis/EqualOffer_ManuscriptAligned/`, and run from the project root:

```powershell
$env:UG_PYTHON = 'C:/Users/neuro-lab/AppData/Local/Programs/Python/Python311/python.exe'
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/run_robustness.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_NEW
& $env:UG_PYTHON analysis/EqualOffer_ManuscriptAligned/summarize_robustness.py results/EqualOffer_Robustness_NEW
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/verify_robustness_refits.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_NEW
& $env:UG_PYTHON analysis/EqualOffer_ManuscriptAligned/verify_robustness_interpretation.py results/EqualOffer_Robustness_NEW
```

Use your own interpreter paths on another computer. The existing project R
library path is intentionally tied to R 4.3 Windows; review it when porting.
Each command must succeed before the next. NEW must name a nonexistent directory;
never rerun into an existing output. Running again means an intentional new run,
not a reason to replace previously unfavorable diagnostics.

## Outputs

- ROBUSTNESS_REPORT.md: fit coverage, all marginal contrasts and simulation checks.
- fit_attempts.csv: both optimizers for each deletion/actor model, including warnings.
- sensitivity_status.csv: complete valid/failed deletion counts and identifiers.
- marginal_effects.csv: estimates, SE, unadjusted 95% Wald CI, OR and BH p-values.
- influence_summary.csv: full-sample versus valid deletion estimates, with source IDs.
- actor_comparison.csv: valid actor fits only; an empty table means no valid fit.
- covariance_diagnostics.csv: all attempted covariance blocks, not only valid fits.
- simulation_summary.csv and E*_simulation_draws.csv: all fixed-parameter draws and
  pointwise descriptive envelopes, with seed and observed values.
- unchanged_primary_routes.csv: primary omnibus-based decomposition routing.
- stat_source_audit.csv: estimates and explicit claim scope/source row keys.
- verification_checks.csv / python_verification_checks.csv: independent arithmetic
  and input verification; input_manifest.csv checks source files remain unchanged.
- independent_refit_comparison.csv / REFIT_VERIFICATION.txt: exact targeted
  selected-model refits from independently reconstructed canonical inputs.
- validation_claim_audit.csv: source-linked verification of the delivered
  interpretation. This last check asserts the findings of this specific primary
  run; changed data or primary models may legitimately fail those assertions.

For the delivered run, start with INTERPRETATION.md, then ROBUSTNESS_REPORT.md.
The interpretation note is an AI-authored assessment checked by the final audit;
it is not automatically rewritten by the numerical scripts.

Failed deletion fits limit influence coverage and are not null effects. Numerical
replication is not proof of distributional validity, practical importance or a
psychological mechanism. No participant is removed from the primary analysis.
