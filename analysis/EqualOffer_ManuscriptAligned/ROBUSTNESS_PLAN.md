# Selected-model robustness checks

Recorded before running these checks, following Jin's authorization on 2026-09-28.
This is retrospective exploratory analysis, not a preregistration. Primary source:
`results/EqualOffer_Aligned_20260928_v2`. No primary selection is repeated.

1. Delete each participant in turn, retaining the exact selected fixed/random
   formula. Retain all participants in the primary analysis. Fit with the original
   BFGS and nlminb controls; select the lowest-AIC candidate passing the original
   numerical gate. Retain failed attempts and warnings; do not simplify formulas
   to obtain a successful deletion fit. Compare all ten marginal emotion pairs,
   the allocation contrast, and the three omnibus tests with the primary model.
   Report coefficient changes in units of the primary SE, directions and intervals.
   These descriptive influence measures have no universal pass/fail threshold.
2. Add `(1 | actor_id)` to the selected participant structure as a separate
   sensitivity fit. Apply the same numerical gate and dual-optimizer rule.
   This checks an actor intercept only, not actor-specific emotion slopes or full
   generalization across all possible faces. Do not interpret invalid candidates.
3. Generate 1,000 Bernoulli data sets from each saved primary glmmTMB model using
   fixed seeds (20260928 for E1; 20260929 for E2). Keep the observed trial design,
   sample new random effects from fitted distributions, and hold fitted parameters
   fixed. Compare ten allocation-by-emotion rejection totals, the overall rejection
   total, the number of zero-rejection participants, and the largest participant
   share of rejections. Retain all simulated summaries. Pointwise 2.5/97.5 percentiles
   are descriptive simulation envelopes, not simultaneous confidence intervals,
   posterior predictive intervals, or formal multiplicity-adjusted tests. Model
   parameter/selection uncertainty and trial-order dependence are not assessed.

Preserve the original reporting gate: interaction p<.05 permits main decomposition;
.05<=p<.10 permits supplementary exploratory decomposition; p>=.10 leaves it audit-only.
Sensitivity results cannot unlock a primary decomposition gate. Marginal emotion
contrasts average equally over the two allocations on the log-odds scale. BH families
remain the ten emotion pairs separately by experiment/model. The sole allocation
comparison needs no multi-pair adjustment. All Wald intervals remain unadjusted.
Deletion fits are descriptive diagnostics, not new opportunities for discovery.

Verification: hash source files before/after; verify exported trial keys, outcomes
and model frames; reproduce primary marginal contrasts from model-matrix algebra;
independently verify correction/interval arithmetic, complete deletion coverage,
and simulation counts. Save versions and sessionInfo. No manuscripts, figures,
raw data, primary models or existing results are overwritten.

## Local execution

From the project root (existing R 4.3.3 project library; no installation):

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/run_robustness.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_NEW
```

The destination must not exist. Both experiment loops and optimizer attempts are
finite. Each attempt has an R elapsed-time guard of 300 seconds (native routines
may return control late); failed/timeout attempts remain in the diagnostics table.
The original optimizer iteration settings are unchanged. Expect several minutes.
