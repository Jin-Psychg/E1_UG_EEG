# Actor-preserving sensitivity analysis

This folder implements Jin's approved second-stage assessment of actor identity
in fair-only choice data. Read `ANALYSIS_PLAN.md` before running. The analysis is
retrospective and exploratory. It does not replace the primary manuscript model.

## Required local project

- Run from `C:/Code/UG_ERP_Project` with R 4.3.3 and the existing project renv library.
- Python 3 with numpy, pandas and scipy is needed for independent table checks.
- Primary input: `results/EqualOffer_Aligned_20260928_v2` (trial exports and models).
- Audited helper source: `results/EqualOffer_Robustness_20260928_v1/delivery_code`.
- Canonical source trials in `data/02_Pipeline_Output_E1` and `E2` are read-only.
- Original engine in `analysis/EqualOffer_ManuscriptAligned/manuscript_engine.R`.

No dependencies are automatically installed. Missing dependencies or changed
source hashes stop execution. A fresh machine needs the project data, existing
input runs, and the matching R package environment. The code ZIP contains no
participant data and is not a standalone data release.
Recorded input manifests use absolute paths. These commands were tested in the
stated Windows project layout. Relocation requires an explicit, documented
manifest re-verification; do not bypass source checks or edit hashes to force a run.

## Reproduce in a NEW output directory

PowerShell commands, from the project root (change the final run name if it exists):

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ActorStructure/run_actor_structure.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_20260928_v1 results/EqualOffer_ActorStructure_replication
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ActorStructure/verify_actor_structure.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_ActorStructure_replication
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ActorStructure/verify_correlation_table.R results/EqualOffer_ActorStructure_replication
& 'C:/Users/neuro-lab/AppData/Local/Programs/Python/Python311/python.exe' analysis/EqualOffer_ActorStructure/summarize_actor_structure.py results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_ActorStructure_replication
```

Execute commands in order and stop if one fails. All output paths must be new;
do not delete old results to bypass this protection. The checked-in adapted
`actor_engine.R` is ready to use. `build_actor_engine.py` documents its generation
from the original engine; do not regenerate it over the supplied copy.

## Statistical contract

Trial outcome is rejection (1) versus acceptance (0). Keep 5:5 and 6:4 trials,
the original 300-3000 ms inclusion bounds, sum coding, and all fixed terms in
`emotion * allocation`. The actor random intercept remains in every structure
candidate. Original grouped participant slope selection uses LRT alpha=.20,
followed by the original correlation selection rules and numerical acceptance
gate. Extended BFGS/nlminb budgets are stated in the plan and code.

After full-sample selection, remove only actor to obtain a matched no-actor
model. Only if both fit acceptably, compare estimates and delete each participant
in turn from the selected actor model, keeping that formula fixed. Selection or
fit failures are retained; no alternative structure is silently substituted.

Type III Wald omnibus tests retain the interaction. Emotion contrasts average
equally across allocations on the log-odds scale; ten contrasts per experiment
and model receive BH correction. The sole allocation contrast averages equally
across emotions. Its OR compares rejection odds at 5:5 with those at 6:4.
Confidence intervals are unadjusted Wald intervals. Do not use these marginal
contrasts as evidence of simple effects at 5:5.

Original primary decomposition routing is unchanged: interaction p<.05 permits
main decomposition; .05<=p<.10 supplementary exploratory decomposition; p>=.10
audit-only. This is a project convention, not a universal statistical rule.
Sensitivity results cannot unlock that gate. All new tests remain exploratory.

## Files to inspect

- `REPORT.md`: statistical tables and explicit interpretation limits.
- `INTERPRETATION.md`, when supplied: reviewed current-run summary, not generated
  automatically for a replication run.
- `selection_status.csv`: whether each experiment reached the matched comparison.
- `E*_selection.log`, `E*_selection_stages/`: structure selection and saved stages.
- `fit_attempts.csv`: both optimizers, including failed attempts and diagnostics.
- `marginal_effects.csv`, `omnibus.csv`: valid-model estimates and tests.
- `matched_comparisons.csv`: actor/no-actor comparison with participant structure
  held constant; original primary shown separately.
- `deletion_coverage.csv`, `deletion_influence.csv`: coverage and valid-only ranges.
- `covariance_diagnostics.csv`, `selected_covariance_summary.csv`: variance and
  eigenvalue diagnostics, beyond the original per-variance acceptance gate.
- `stat_audit_*.csv`: exact result values and row-key provenance.
- `verification_checks.csv`, `python_verification_checks.csv`: structure and
  arithmetic checks; `REFIT_VERIFICATION.txt` records independent full-sample
  reproduction from canonical input. Deletions are not independently repeated.
- `actor_engine.diff`, `engine_provenance.json`: all adaptations from the original.
- `correlation_selection_verified.csv`: checked correlation LRT table. The raw
  legacy `df_diff` column contains candidate parameter count, not the parameter
  count difference. The verified export labels it correctly and computes the
  actual difference; p-values and model selections are unchanged.

An actor intercept accounts for baseline actor variation only. It does not model
actor-specific expression effects, prove generalization to all faces, or identify
psychological mechanisms. Passing numerical checks is not a scientific validity
certificate. Read failed-fit coverage and uncertainty alongside effect estimates.
