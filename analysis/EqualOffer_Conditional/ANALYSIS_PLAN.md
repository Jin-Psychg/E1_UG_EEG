# Equal-offer conditional analysis plan

Recorded before computing the new conditional tests, 2026-09-30. Exploratory;
previous descriptive and pairwise outputs have already been seen. This is not
preregistration and does not remove selection or multiplicity from earlier work.

Question: Are the five expression-specific log odds of rejection equal at 5:5?
Use the existing full emotion * allocation Bernoulli-logit GLMM, retaining its
participant random structure, all trials at 5:5 and 6:4, and original exclusions.
Test four independent conditional contrasts jointly (4-df Wald chi-square).
Do not replace this question with the interaction test or refit/select a
5:5-only model. Do not interpret absence of interaction as equivalence.

Primary family: two conditional omnibus tests, E1 and E2; BH-FDR across these
two p-values, alpha .05. Record raw p-values as well. No within-ratio pairwise
significance tests in this stage. Export five conditional rejection probabilities
with pointwise 95% Wald intervals, with random effects set to zero; these are
not population-averaged probabilities or raw participant means.

Sensitivity: extract the identical conditional estimand from saved accepted
matched-no-actor and actor-selected models, and all available participant-deletion
models (original structure plus accepted numerical retries; actor-selected structure).
Keep alternative specifications separate. No new random-structure selection,
no repair or reinterpretation of historically failed fits. Show deletion coverage,
nominal p-value ranges, changes in probabilities, and high-influence participants.
Sensitivity p-values are diagnostic estimates under alternative data/specifications,
not independent replications. BH within each full-sample two-experiment sensitivity
pair is descriptive and does not enter the primary family.

Before extraction check input hashes, selected-model diagnostic gate, retained
participants/outcomes and fixed design. Independently verify joint Wald arithmetic
against emmeans::test(joint=TRUE), and recheck input hashes after extraction.
Cross-check cell counts against trial exports and saved descriptive tables.

Limits: asymptotic Wald inference, sparse expression cells, selection uncertainty,
incomplete deletion coverage, and conditional rather than integrated probabilities.
No bootstrap calibration or Bayesian fit is performed in this bounded stage.
Do not label a finding robust solely because the full fit converged. If evidence
is specification-sensitive, retain the result with that qualification and do not
use it to choose a psychological theory. Theoretical comparisons remain in the
existing manuscript framework; no new mechanism claim is authorized.
