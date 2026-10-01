# Reproduce the conditional equal-offer checks

This folder addresses one retrospective question: whether the five expressions
differ in rejection probability at exactly 5:5. Read ANALYSIS_PLAN.md first.
It extracts new conditional tests from existing, validated model objects; it does
not refit models or select a new random structure. No new pairwise significance
tests are performed. The primary Wald p-values use BH-FDR across E1 and E2.

## Run locally

From the UG_ERP_Project root, using R 4.3.3 and the existing renv library:

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_Conditional/run_analysis.R results/EqualOffer_Conditional_NEW
python analysis/EqualOffer_Conditional/verify_and_summarize.py results/EqualOffer_Conditional_NEW
```

Use your local Python executable with numpy and pandas installed. No scipy is
needed: the 4-df chi-square survival probability is verified by its closed form.
The output directory must not exist. Stop if either command fails. Raw data and
historical model/output folders remain read-only.

## Required existing inputs

- EqualOffer_Aligned_20260928_v2: primary models and analysis trial exports.
- EqualOffer_Robustness_20260928_v1: accepted original-structure deletion models.
- EqualOffer_NumericRetry_20260929_v3: accepted additional deletion models.
- EqualOffer_ActorStructure_20260929_v1: actor-selected, matched-no-actor, and
  accepted actor-structure deletion models.
- Canonical trial CSVs and the manuscript_engine.R snapshot: hash/frame checks.

Exact paths, hashes and missing deletion models are recorded in each run.
The scripts are readable and reusable locally; independent researchers need
access to these study inputs to reproduce the numerical results. No private
trial data or participant model objects are copied into the code folder.

## Read the outputs

- conditional_omnibus.csv: 5:5 joint chi-square, raw p and two-experiment BH p.
  Primary, alternative structures and deletions are identified separately.
- conditional_probabilities.csv: five rejection probabilities and pointwise
  95% Wald intervals, with random effects set to zero. These differ from pooled
  observed percentages and from population-averaged predictions.
- equal_offer_counts.csv: events, denominators, contributing participants and
  the participant with the largest event count in each expression condition.
- deletion_summary.csv: accepted/attempted coverage and nominal p ranges.
- model_roster.csv and diagnostics.csv: missing models and acceptance checks.
- stat_source_audit.csv, input_manifest.csv, verification.json: provenance and
  independent checks. SESSION_INFO.txt records package versions.

Do not infer specific expression-pair differences from the omnibus test.
Do not treat deletion fits as independent replications. Finite-sample calibration
of asymptotic Wald tests and random-structure selection uncertainty are not
resolved by these diagnostics. Historical failures remain missing coverage.

## Current completed run

results/EqualOffer_Conditional_20260930_v3 passed extraction and independent
summary verification. v1 stopped because the verification compared full-precision
arithmetic with emmeans' internally rounded chi-square display. v2 stopped because
a saved deletion object's original dataset name was unavailable. v3 checks rounded
chi-square plus exact p, and supplies each validated stored model frame explicitly.
Neither repair changed fitted models, the hypothesis, factor coding or multiplicity.
The two debug scripts document these interface checks and are not pipeline steps.
