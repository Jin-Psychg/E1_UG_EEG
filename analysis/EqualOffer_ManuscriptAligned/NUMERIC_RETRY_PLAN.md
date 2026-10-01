# Bounded numerical retries, unchanged formulas

Recorded 2026-09-29 before new retry outcomes. This implements the first step of
ACTOR_REPAIR_PROPOSAL.md only. No structural simplification is authorized here.

Sources: primary EqualOffer_Aligned_20260928_v2; diagnostics from
EqualOffer_Robustness_20260928_v1. Retry each experiment/model for which neither
previous optimizer passed the original gate, including actor extensions. Preserve
all primary fits and already valid sensitivity fits; do not choose targets by p.

- Fixed/random formulas, all factor coding, exclusions and participants are unchanged
  for each target. A deletion target still omits only its original designated ID.
- BFGS: maxit=2000. nlminb: iter.max=2000, eval.max=4000. Both start from the package's
  default initialization, as in the original check. No warm starts or tolerance
  changes. Save both complete fit objects, including failed diagnostic-only objects.
- Use the original optimizer/Hessian/finite-SE/likelihood/variance gate. Select the
  lowest-AIC valid candidate as before. Do not change the 1e-6 variance rule to
  obtain a pass. Save objective/parameter differences and covariance eigenvalues
  for both attempts; agreement is diagnostic, not a new significance test.
- Each attempt has a 300-second R elapsed-time guard (native routines may return
  late), finite iteration caps, and one CPU thread. There are two attempts per
  target, with no open-ended retry loop.
- Recovering an original-gate pass does not prove nonsingularity or full robustness.
  The original diagonal-variance check may miss near-zero covariance directions.
  Report these separately. Never interpret invalid-fit fixed effects as findings.
- Export all marginal emotion/offer contrasts and omnibus tests for valid recovered
  fits using unchanged correction families. They cannot alter primary reporting
  gates or replace the primary models. No new simple-effect decomposition is run.
- If retries remain unsuccessful or expose structural degeneracy, stop at diagnosis
  and report the next methodological decision. Do not automatically simplify models.

Verification: source hashes, exact target coverage/formulas/data, independent
matrix calculations for contrasts/SEs, separate interval/BH verification and
reconstruction of at least one recovered fit if available. Retain failures.

```powershell
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla analysis/EqualOffer_ManuscriptAligned/run_numeric_retry.R results/EqualOffer_Aligned_20260928_v2 results/EqualOffer_Robustness_20260928_v1 results/EqualOffer_NumericRetry_NEW
```

Run from the project root using the existing environment. Output must not exist.
The runner reuses hash-verified helper function definitions from the earlier
delivery_code/run_robustness.R snapshot; it does not execute that script's fits.
