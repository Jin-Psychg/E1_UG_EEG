# Manuscript-aligned equal-offer follow-up

This version implements Jin's approved consistency revision after the earlier
Holm analysis was examined. It remains retrospective and exploratory. It retains
rejection as the binary outcome: complementary acceptance coding supplies no new
information for the same logit model. No manuscript or theory record is edited.

## Current reporting policy (supersedes the earlier exception)

Jin explicitly requires the original manuscript reporting gate. The current report
is `results/EqualOffer_Aligned_20260928_v2/standard_reporting/REPORT.md`, with a
separate SUPPLEMENT.md. The prior REPORT.md and code archive remain historical
audit artifacts and no longer define the reporting hierarchy. Main emotion
comparisons average equally across offer-ratio levels; the two levels remain
separate in the fitted model. E1 decomposition is supplementary; E2 decomposition
is audit-only. This decision was made after inspecting results and is not preregistered.

The default run_all.R now creates the paired figure and standard report. To apply
only the reporting change to an already completed run without refitting:

```sh
python analysis/EqualOffer_ManuscriptAligned/render_standard_report.py results/EqualOffer_Aligned_NEW
```

This refuses an existing standard_reporting directory. Stored estimates, p-values,
CIs, fits and figures remain unchanged. The old render_report.py is retained solely
for reproducing the earlier explicitly superseded report, not as the default.

## Local execution

From the UG_ERP_Project root, using the existing R 4.3.3 renv library:

```powershell
$env:UG_PYTHON = 'C:/Users/neuro-lab/AppData/Local/Programs/Python/Python311/python.exe'
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/run_all.R
```

The Python path above is this machine's tested interpreter; use your own interpreter
on another machine. Required Python packages are numpy, pandas and matplotlib.
No package installation, updates, raw-data edits or original-model overwrites occur.

For explicit stage execution, replace NEW with a new run name:

```sh
Rscript --vanilla analysis/EqualOffer_ManuscriptAligned/run_analysis.R --out=results/EqualOffer_Aligned_NEW
Rscript --vanilla analysis/EqualOffer_ManuscriptAligned/verify_results.R results/EqualOffer_Aligned_NEW
python analysis/EqualOffer_ManuscriptAligned/plot_paired.py --run results/EqualOffer_Aligned_NEW
python analysis/EqualOffer_ManuscriptAligned/render_standard_report.py results/EqualOffer_Aligned_NEW
```

The first stage rejects an existing directory. Verification rejects prior check
outputs. Rendering rejects an existing REPORT.md. A later intentional rerun must
use a new directory. Keep earlier versions for audit.

## Provenance and method contract

Read ANALYSIS_PLAN.md for the analysis decisions recorded before inspecting this
version's contrasts. The source engine snapshot is extracted without substantive
changes from the current manuscript Rmd; engine_provenance.json records both
hashes. The runner checks them before fitting. If the original manuscript script
changes, review the changes rather than merely updating its expected hash.

The unchanged sibling EqualOffer_Exploratory/run_analysis.R is run in audit mode
to rebuild canonical cleaned inputs, verify their MD5s against the original E1/E2
formal runs, and check saved original model frames. The new runner then uses the
original five-step selection with allocation replacing pooled offer_type.
The original RT thresholds, response mapping, factors and IDs are preserved.

This requires the complete local project: canonical trial CSVs, the two formal
manuscript model runs, the preserved Holm replication output used for comparison,
the original Rmd, and the sibling audit script/config. These private inputs are not
included in the code archive. See the sibling config.R for explicit run identities.

BH correction families include all ten marginal expression pairs per experiment,
all ten expression pairs within each experiment and allocation, and a separate
ten-test interaction family per experiment. Only the original omnibus reporting
route determines the location of decomposition results. The earlier exception is
superseded, as specified above. No effect becomes preregistered because it passes
a statistical threshold. CI coverage is unadjusted.

## Selected-model robustness checks

Read ROBUSTNESS_PLAN.md first. These checks preserve the selected primary formula,
BH families and reporting gate. They do not rerun selection, delete participants
from the primary analysis or promote simple contrasts based on sensitivity results.
Use a new destination for each intentional rerun:

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/run_robustness.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_NEW
& $env:UG_PYTHON analysis/EqualOffer_ManuscriptAligned/summarize_robustness.py results/EqualOffer_Robustness_NEW
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/verify_robustness_refits.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_NEW
```

Set UG_PYTHON as in Local execution above. The summarizer also requires scipy and
pandas, present in the tested local environment. The runner saves every fitting
attempt and warning. It interprets only models passing the original gate; a failed
deletion fit is not a null effect and limits the coverage of the influence check.
All estimates and CIs are in marginal_effects.csv, with ROBUSTNESS_REPORT.md providing
the overview. The final verification reconstructs canonical inputs and repeats the
valid actor fit and the largest standardized-change deletion per experiment; it
does not repeat the entire deletion sequence or simulation.

## Outputs and how to read them

- REPORT.md: methods, all previously specified comparisons, omnibus tests, diagnostic
  limits and manuscript-style figure.
- new_contrasts.csv: all 60 simple/interaction contrasts; primary p column is p_BH.
- new_omnibus.csv, emotion_main_contrasts.csv, E*_allocation_main.csv: omnibus and
  main-effect outputs. Interpret main effects in the context of interactions.
- version_comparison.csv: old model/old Holm, new model/old Holm, new model/new BH.
  The middle column is a diagnostic bridge, not a competing primary test.
- stat_source_audit.csv: exact-value comparison and source locators.
- E*_aligned_selection.log/json/rds: every selection step and final diagnostics.
- covariance_diagnostics.csv: covariance eigenvalues and correlations, supplementing
  the original acceptance gate without changing selection after results were seen.
- verification_checks.csv and independent_refit_comparison.csv: independent trial
  counts, model-matrix contrasts, interval/probability checks, BH arithmetic and
  deterministic exact-formula refits. Full selection is not repeated by verification.
- paired_figures/: PNG/PDF/SVG/TIFF, caption and source coordinates. Participant points
  match the earlier approved figure; model estimates use this version's fits.

The existing_* files, run_analysis.R/config.R at the output root and original
REPRODUCIBILITY.txt belong to the legacy audit stage. They are not the new aligned
model results. Use ALIGNED_REPRODUCIBILITY.txt and aligned_code/ for this revision.
The delivered run additionally contains delivery_code/, with the final verification,
rendering and wrapper files completed after fitting began. aligned_code/ preserves
the snapshot at fit start; these two timestamps are not interchangeable.

## Delivered run and checks

The delivered run is results/EqualOffer_Aligned_20260928_v2/REPORT.md. Both final
models passed the original acceptance gate. Eight hundred independent checks
passed; coefficients, covariance and log likelihood exactly reproduced for both
selected models on this machine. Figure coordinates exactly match the approved
prior figure. The run_all.R wrapper's component stages were executed separately;
the wrapper itself has been syntax-checked, not end-to-end rerun.

Development v1 stopped after the E1 fit because the contrast exporter read factor
levels from the wrong environment. v2 fixes that interface; both runs and logs
are retained. A first report rendering was retained separately after improving
floating-point/literal-value preservation; the model result files were unchanged.

Important limitation: the original actor-intercept, participant-deletion and
simulation diagnostics apply to the older fits. They were not rerun for the new
structures. Some covariance directions are weak despite passing the original
gate. Neither exact reproducibility nor BH significance proves full robustness or
identifies WTR, AV, status threat or inferred intentions.
