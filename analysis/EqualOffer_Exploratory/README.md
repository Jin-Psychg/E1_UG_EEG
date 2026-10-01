# Equal and mildly unequal offers: exploratory reproducibility package

## Purpose and status

This analysis was requested on 2026-09-28 after reviewing the manuscript results.
It is **retrospective and exploratory**, not preregistered. It asks whether facial
expression differences in rejection are present at 5:5 and whether they differ
from 6:4. It does not measure perceived WTR, association value, status threat, or
intentions, and cannot identify a psychological mechanism from a statistical interaction.

The script reads the live project's canonical preprocessing outputs, preserves
the original analyses, and writes to a **new** results folder. It does not run EEG
preprocessing, HDDM, SVO analyses, or the original statistical pipeline.

## Run locally

Open `UG_ERP_Project.Rproj`, then use the RStudio Terminal from the project root:

```sh
Rscript --vanilla analysis/EqualOffer_Exploratory/run_all.R
```

Or run `source("analysis/EqualOffer_Exploratory/run_all.R")` in the RStudio console.
This creates a fresh timestamped results directory and runs analysis, independent
verification, and report/figure rendering in sequence. For individual stages:

```sh
Rscript --vanilla analysis/EqualOffer_Exploratory/run_analysis.R --stage=audit --out=results/EqualOffer_audit_new
Rscript --vanilla analysis/EqualOffer_Exploratory/run_analysis.R --stage=all --out=results/EqualOffer_full_new
Rscript --vanilla analysis/EqualOffer_Exploratory/verify_results.R results/EqualOffer_full_new
Rscript --vanilla analysis/EqualOffer_Exploratory/render_report.R results/EqualOffer_full_new
```

On this Windows machine, if Rscript is not on PATH, PowerShell accepts:

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_Exploratory/run_all.R
```

Choose a fresh output name on every run. An existing directory is rejected, not
overwritten. Paths derive from the project root; moving the entire project does
not require editing script paths. `--root=PATH` is available for unusual launch contexts.
The output must be inside this project's `results/` directory.

Use R 4.3.3 and the existing renv project library. On a fresh machine, inspect the
project's `renv.lock` and run `renv::restore()` from R as needed. The script never
installs or upgrades packages. `SESSION_INFO.txt` records the actual versions;
it is authoritative if the historical lock file lags the local library.
Participant-level data are needed locally and are not included in this code folder.

## Visualization-only update (2026-09-28)

Jin requested manuscript-consistent figures and explicitly asked to retain the
current statistics for this turn. No model was refitted and no p-value adjustment
was changed. Read [METHODS_COMPARISON.md](METHODS_COMPARISON.md) for the source-checked
comparison of the original five-step / BH-FDR workflow and this exploratory
diagonal-model / Holm workflow. Numerical reproducibility does not make the two
methodologies identical.

The report renderer now uses the manuscript's Python paired-facet figure design,
via `plot_paired.py` and the self-contained `manuscript_style.json` palette snapshot.
Python needs numpy, pandas and matplotlib (locally tested with 1.26.4, 2.2.2 and
3.9.0). R still performs all statistics. If Python is not on PATH, set `UG_PYTHON`
to its executable before running the R report renderer. No package installation
or environment update is performed automatically.

To **update figures only**, without invoking `run_all.R` or fitting anything:

```sh
python analysis/EqualOffer_Exploratory/plot_paired.py --run results/EqualOffer_replication_20260928_224854 --out results/EqualOffer_replication_20260928_224854/paired_figures_new
```

The output folder must be new. The delivered revised figure is in
[`paired_figures/`](../../results/EqualOffer_replication_20260928_224854/paired_figures/FIGURE_NOTES.md).
The alternative report view is
[`REPORT_manuscript_style.md`](../../results/EqualOffer_replication_20260928_224854/REPORT_manuscript_style.md).
Original figures and the original report are preserved. Future report generation
uses the paired figure by default. Thin lines connect the same participant across
5:5 and 6:4, diamonds show conditional model predictions with unadjusted 95% CIs,
and the exact manuscript palette/order is retained. The figure does not merge
observed participant means and model predictions into the same statistic.

## Inputs and provenance gates

1. `data/02_Pipeline_Output_<E1|E2>/Method_Regression/Stimulus_Locked/trials.csv`.
2. The two explicit manuscript reference runs in `config.R`, selected from
   `E:/MyPaper_withARS/draft/round180_formal_script_revision/Stat_Source_Map.md`
   and `check_E2_all.R`. They are not selected by newest timestamp.
3. The saved choice and all-trial RT models from those runs.

The input CSV MD5 must match the value recorded in the original run's
`REPRODUCIBILITY.txt`. The cleaned choice outcome, emotion, pooled fairness,
participant IDs, and row count must match the saved model frame in order.
The manifest records MD5 and SHA256 for source data, saved models, current
original script, new script/config, and environment lock file. A changed input
causes a stop for provenance review; do not simply replace expected hashes.

Original inclusion is preserved: offers 5, 6, 8, 9; reaction 1/2; RT 300--3000 ms,
inclusive. The new models use only 5 and 6. Reaction 1 means accept; 2 means reject.
`Offers_Other` is the proposer share: 5 is 5:5 and 6 is 6:4. Original IDs, numeric
responses and trial indices remain in the audit export. The original loader's
internal participant-ID mapping is a separate, checked, exported column.
**These are preprocessing-output trials**; exclusions before that export are
not reconstructed from raw logs. This analysis does not generalize its counts
to all behavioral trials originally collected.

## Statistical specification (fixed before new focal results)

- Separate E1 and E2 Bernoulli-logit mixed models:
  `reject_binary ~ emotion * allocation`.
- All five expressions remain separate; allocation is categorical (5:5, 6:4).
  All participants, including zero-rejection participants, remain in the data.
- Participant intercept plus all nine sum-coded emotion, allocation and interaction
  slopes, initially uncorrelated. This transparent diagonal structure reduces
  covariance parameters without presuming that within-person effects are fixed.
- Fit by glmmTMB Laplace maximum likelihood. Try nlminb, then BFGS on numerical
  failure. Remove only slope SDs below 0.001, smallest first, then refit. The
  intercept is not removed. No fixed-effect p-value enters this procedure.
- A fit must have optimizer code 0, positive-definite Hessian, finite coefficients,
  SEs and likelihood, and no flagged extreme coefficient (absolute value >15) or
  SE (>10). Thresholds are numerical flags, not proof of no separation.
  If no valid fit is obtained, stop new inference. A regularized Bayesian analysis
  would then need its own explicit priors and checks; it is not silently substituted.
- Each fit has an evaluation/iteration cap and an R-level timeout. Native compiled
  optimizers may not respond immediately to an R timeout; the iteration cap still applies.
- Primary contrast family: Dominance-Affiliative and Dominance-Reward, each within
  5:5, within 6:4 and the difference **6:4 minus 5:5**. Holm adjustment across all
  12 log-odds tests (both experiments).
- Secondary family: Disgust-Reward, the same three contrasts across two experiments;
  Holm across six tests. The secondary label limits multiplicity scope, not importance.
- Model-scale estimates are log-odds differences (OR after exponentiation).
  For interaction contrasts, exponentiation gives a **ratio of ORs**, not a simple OR.
  Probability differences and differences-of-differences are also exported but
  are not interchangeable with the logistic interaction. Probability-scale p-values
  in the mechanical export are not the designated tests; `p_holm` is intentionally NA.
- All CIs are **unadjusted 95% Wald intervals**. Holm-adjusted p-values can disagree
  with whether these unadjusted intervals include zero. Intervals are conditional
  on the selected random structure and omit model-selection uncertainty.
- Predicted probabilities set random effects to zero. They are not integrated
  population-average probabilities; show them alongside observed participant and
  trial-weighted rates. Do not label them as observed rejection percentages.
- Sensitivity: add a crossed actor intercept to the selected participant structure.
  Actor IDs are parsed using the existing project's validated stimulus convention.
- Influence: refit the same selected formula after omitting each participant in
  turn. Invalid deletion fits are reported, not interpreted. No participant is
  removed from the primary analysis in response to the result.
- Parametric model check: 200 Bernoulli datasets simulated with new random effects;
  compare cell rejection counts and the number of zero-rejection participants.
  These are descriptive checks at the fitted parameters, not posterior predictive
  intervals or formal tests of overall model adequacy.

The existing-model audit retains original choice Wald chi-square tests and RT
Satterthwaite F tests. Its new targeted contrasts use the original asymptotic
emmeans convention and a separate nine-test Holm family per experiment/outcome.
These audit contrasts are not renamed as the original preregistered tests.

## Outputs

| File | Meaning |
|---|---|
| `STATUS.txt` | COMPLETE only if the requested stage finished |
| `input_manifest.csv`, `SESSION_INFO.txt` | Source fingerprints, code versions, actual environment |
| `manuscript_frame_checks.csv` | Trial-by-trial consistency with original saved choice models |
| `exclusions.csv` | Mutually exclusive exclusion reasons within preprocessing outputs |
| `participant_cells.csv`, `participant_rejection_totals.csv` | Event counts, zero rejecters, missing cells |
| `descriptive_cells.csv` | Condition counts, rejection prevalence and concentration |
| `existing_omnibus.csv`, `existing_targeted_contrasts.csv` | Existing choice/RT evidence, no original refits |
| `E*_analysis_trials.csv`, `E*_id_mapping.csv` | Reproducible new-model inputs, original keys preserved |
| `E*_selection_log.json`, `E*_model_summary.txt`, `E*_fair_choice.rds` | Full fitting decisions and saved models |
| `new_contrasts.csv`, `new_probabilities.csv` | New effect estimates, uncertainty and multiplicity adjustments |
| `participant_influence.csv`, `sensitivity_diagnostics.csv` | Robustness and invalid-fit disclosure |
| `Methods_paragraph_E*.md`, `REPRODUCIBILITY.txt` | Methods and reporting boundaries |
| `verification_checks.csv`, `independent_refit_comparison.csv` | Independent counts, contrast algebra and deterministic refits |
| `REPORT.md`, `paired_figures/equal_offers_paired.*` | English report and manuscript-style participant pairs plus conditional predictions; older runs retain the earlier plots |
| `participant_influence_summary.csv`, `stat_source_audit.csv` | Deletion stability and exact statistical source mapping |

Code snapshots are saved with each run. A checkpoint CSV during fitting is
unadjusted and incomplete: use the final CSV only when `STATUS.txt` says COMPLETE.
Outputs include participant-level data; review data-sharing permissions before
sharing output folders. The scripts and README can be read independently.

## Local validation record

The delivered, end-to-end tested run is
[`EqualOffer_replication_20260928_224854`](../../results/EqualOffer_replication_20260928_224854/REPORT.md).
It passed 812 independent checks and exact refits of the two selected primary
models. Ten key result CSVs are byte-identical to the independently verified
development run `EqualOffer_full_20260928_v4`; see `end_to_end_replication.json`.
Earlier development outputs are preserved for audit:
`EqualOffer_audit_20260928_v1` contained an aggregation bug subsequently caught
and fixed; `EqualOffer_full_20260928_v1` stopped at a package print-method issue;
`v2` completed before additional boundary/simulation checks; `v3` stopped at a
binomial-simulation shape mismatch. **Use the delivered replication or a newly
completed and verified replication, not those development runs.** No source data or original results
were modified to fix these issues. Console logs are in this script folder.

The independent verifier reconstructs counts with base R, recalculates contrast
estimates/SEs with model matrices rather than emmeans, checks Holm families, and
refits selected models from the canonical CSVs. It does not repeat the entire
random-structure search or all sensitivity fits; `run_all.R` repeats those too.
Numerical reproducibility does not prove that the model assumptions or theoretical
interpretation are correct. Local exact equality may not hold across software versions.

## Interpretation boundaries

An effect at 5:5 constrains an explanation based solely on monetary inequality.
It does not establish high perceived WTR, an AV effect, humiliation, or inferred
hostile intent. A nonsignificant interaction is not evidence of additivity or
equivalence. Significance at one allocation and not another does not establish
moderation. No theory ranking or manuscript revision is automated.
