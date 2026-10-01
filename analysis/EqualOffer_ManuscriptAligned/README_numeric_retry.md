# Reproduce the bounded numerical retry stage

This stage retries only targets that failed both optimizers in the first
robustness run. Read NUMERIC_RETRY_PLAN.md. Original formulas and acceptance
criteria are unchanged. No random-structure selection is performed.

Use the same project environment documented in README_robustness.md: R 4.3.3 with
the existing project renv library; Python with numpy, pandas and scipy. Inputs
include the canonical trial files, completed primary run, and the earlier
robustness run's diagnostic tables and hash-verified delivery_code snapshot.
The runner loads only named helper definitions from that snapshot, not its fitting
loop. No private data or model fits are distributed in the code ZIP.

From the UG_ERP_Project root, using a NEW nonexistent output directory:

```powershell
$env:UG_PYTHON = 'C:/Users/neuro-lab/AppData/Local/Programs/Python/Python311/python.exe'
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/run_numeric_retry.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_20260928_v1 results/EqualOffer_NumericRetry_NEW
& $env:UG_PYTHON analysis/EqualOffer_ManuscriptAligned/summarize_numeric_retry.py results/EqualOffer_Robustness_20260928_v1 results/EqualOffer_NumericRetry_NEW
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/verify_numeric_retry.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_NumericRetry_NEW
& $env:UG_PYTHON analysis/EqualOffer_ManuscriptAligned/combine_numeric_retry.py results/EqualOffer_Robustness_20260928_v1 results/EqualOffer_NumericRetry_NEW
```

Use local interpreter paths on another machine. Run commands sequentially; stop
if a command fails. Never overwrite a prior result directory. Each model has two
bounded attempts using default starting values. Invalid fit objects in
attempt_models/ are for numerical diagnosis ONLY, not fixed-effect inference.

Delivered run: results/EqualOffer_NumericRetry_20260929_v3. Read INTERPRETATION.md,
REPORT.md, retry_status.csv and selected_covariance_summary.csv together. The
interpretation note is a source-checked AI assessment, not generated prose from
the numerical runner. Exact refit verification covers all recovered valid models,
not all failed attempts. The combined tables retain the original valid fits and
append recovered fits; the combiner verifies an existing table rather than
overwriting it. Primary models and reporting routes are unchanged.

Development v1/v2 stopped during formula object/list attribute comparisons;
v2 formula_checks.csv shows matching formula strings and trial counts. v3 compares
canonical formula text, excluding irrelevant R object attributes. All versions
are retained; v1/v2 are not inferential sources. No numerical model specification
was changed to fix that verification-code issue.
