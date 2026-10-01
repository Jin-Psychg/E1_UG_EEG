# Manuscript-aligned equal-offer analysis

## Subsequent author-approved reporting amendment

The author withdrew the targeted-follow-up exception after the results discussion.
The original manuscript gate now governs reporting: p<.05 main decomposition;
.05<=p<.10 supplementary exploratory decomposition; p>=.10 audit only.
Main emotion comparisons remain equally weighted across offer-ratio levels.
All previously computed results remain preserved; no models are refitted for
this amendment. The earlier plan below is retained as historical provenance.

## Original execution plan (reporting exception superseded)

Status: retrospective exploratory; authorized by Jin after the methods comparison.
The previous Holm results were already seen. This is a consistency-motivated revision,
not a preregistered analysis or independent replication. Preserve the old version.

- Inputs, exclusions, response coding and participant identities: reuse the verified
  canonical audit in EqualOffer_Exploratory/run_analysis.R, audit stage only.
- Separate E1/E2 Bernoulli-logit GLMMs: reject_binary ~ emotion * allocation (5:5/6:4).
- Reuse the unmodified manuscript five-step function, blocked random-slope LRT
  alpha=.20, correlation testing, rePCA, and exact-formula two-optimizer final gate.
- Preserve the full fixed design and all participants, including zero rejecters.
- New omnibus Type III Wald tests; asymptotic Wald contrasts; unadjusted 95% CIs.
- BH families: all 10 expression pairs per allocation per experiment; all 10
  expression-pair allocation interactions per experiment; 10 expression-main
  pairs per experiment. Allocation main comparison uses Tukey (one comparison).
- Export all pairs, highlighting the previously specified dom-aff, dom-enj and
  dis-enj comparisons. No selection of significant pairs for reporting.
- Targeted-follow-up exception: compute these comparisons regardless of omnibus
  interaction p. Record the original p<.05/.10 reporting gate separately. All
  results remain exploratory regardless of that gate.
- If the original acceptance gate fails, retain diagnostics and stop inference
  for that experiment. Do not silently substitute the prior diagonal model.
- Independently reconstruct contrasts/SE/BH values; verify counts, source hashes,
  fitted design, and refit the exact selected model before declaring verification.
- Compare old vs new estimates and p values, separating model change from
  correction/family changes. Current scope is choice only, not RT/SVO/HDDM.
- Keep the manuscript paired-line style. Updated model diamonds must use the new
  estimates; participant points must match the prior version exactly.
