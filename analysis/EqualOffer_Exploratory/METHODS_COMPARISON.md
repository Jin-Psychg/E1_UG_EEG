---
verified: false
scope: Source-checked methods comparison; no new fits or revised p-values
date: 2026-09-28
---

# Comparison with the manuscript's behavioral analysis

## Author instruction for this revision

Jin requested an explanation of Holm correction and a methods comparison, and
explicitly selected **explanation plus figure revision only**. The existing Holm
analysis is preserved. No BH replacement, original-pipeline refit, additional
contrast family, exclusion change, or manuscript/theory revision is authorized
by this update. The figure script reads the existing verified CSVs only.

## What Holm does

Holm is a step-down multiple-testing procedure controlling the family-wise error
rate (FWER): the probability of at least one false rejection among the true null
hypotheses in a defined family. With valid input p-values, it controls FWER under
arbitrary dependence. Order m p-values from smallest to largest, compare them
sequentially with alpha/m, alpha/(m-1), ..., alpha, and stop rejecting at the first
failure. Equivalently, use monotone adjusted p-values and compare each with alpha.
Holm is no more conservative than ordinary Bonferroni.

The manuscript uses Benjamini-Hochberg (BH), called `fdr` in R. It targets the
expected false-discovery proportion among rejected hypotheses, rather than the
probability of any false positive. A 5% FDR target does not mean every individual
significant result has a 5% chance of being false, nor that exactly 5% of findings
in a particular paper are false. BH's standard guarantee depends on independence
or appropriate positive dependence. Neither correction fixes invalid model p-values.

The two methods address different error criteria. Holm is usually more stringent
than BH on the **same input p-values and same family**, but the present analyses
also differ in models and family definitions, so their p-values cannot be compared
as if correction alone differed.

### Source record

- Original source: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
- Inspected scope: official R documentation, Details and References, 2026-09-28.
- Verbatim excerpt supporting FWER: "strong control of the family-wise error rate".
- Verbatim excerpt supporting FDR: "expected proportion of false discoveries amongst the rejected hypotheses".
- `BH` and `fdr` are explicitly identified as aliases on that page.
- Additional source: https://arxiv.org/abs/2201.09350v3 (Ruodu Wang, original abstract HTML inspected on 2026-09-28; full paper not inspected here).
- Verbatim scope excerpt: "BH) procedure for independent or positive-regression dependent p-values". This supports the dependence qualification, not a verification of the dependence structure of this study.
- These excerpts and the source check are AI-generated records; human verification remains false.

## Detailed comparison

| Component | Manuscript / original formal pipeline | New equal-offer script | Assessment |
|---|---|---|---|
| Input provenance | Canonical Method_Regression / Stimulus_Locked trial exports; explicit E1/E2 formal runs | Same canonical inputs, MD5 checks and exact saved-choice-frame matching | Aligned |
| Inclusion | Responses 1/2, RT 300-3000 ms, pooled offers 5/6/8/9; filler 7 excluded | Same cleaning; then retain 5/6 only | Necessary narrowing for the equal/mildly unequal question |
| Outcome | Binary rejection, accept=0/reject=1 | Same | Aligned |
| Fixed effects | `emotion * offer_type`, fair={5,6}, unfair={8,9} | `emotion * allocation`, 5:5 versus 6:4 | Necessary change in scientific estimand; cannot copy the original final formula unchanged |
| Model family / engine | glmmTMB, Bernoulli-logit, ML with Laplace approximation | Same | Aligned |
| Factor coding | Sum-to-zero contrasts; explicit numeric contrast columns for diagonal random slopes | Same coding principle and expression levels | Aligned |
| Candidate random structure | Attempt maximal correlated participant intercept, factor slopes and interaction slopes; then diagonal model | Start directly with diagonal intercept and all nine numeric slopes | A new modeling decision, not a requirement imposed by splitting offers |
| Random-effect reduction | buildmer backward LRT with alpha=.20, contrast blocks for factors; evaluate retained correlations; rePCA dimensionality check | Drop individual slope SDs below .001, smallest first; no LRT reduction, correlation restoration or rePCA | Substantive difference; the new procedure is not the manuscript's five-step procedure |
| Coding dependence | Blocked LRT evaluates factor contrast columns together | Boundary deletion acts on individual sum-coded columns | The new variance restrictions depend on the chosen basis; retain coding and report this limitation |
| GLMM numerical boundary | Original GLMM gate flags variance <1e-6 | New gate flags SD <.001 | Same variance threshold expressed on different scales |
| Final optimizer | Refit exact selected formula with BFGS and nlminb; retain valid candidate with smallest AIC (same formula) | Try nlminb; use BFGS only after numerical failure | Engineering difference; the new script does not independently compare two converged optima |
| Additional numerical flags | Convergence, Hessian, finite SE/likelihood, singularity and dimensionality checks | Same core numerical checks; adds abs(beta)>15 or SE>10 flags | Additional heuristic flags are not a formal separation test |
| Omnibus tests | Type III Wald chi-square for choice; Satterthwaite F for Gaussian outcomes with fallbacks | Existing-model omnibus tests extracted; no new fair-only omnibus test used as a gate | New analysis focuses on predefined-for-this-follow-up contrasts; it is exploratory |
| When simple effects are produced | Script routes interaction p<.05 to main tables; .05<=p<.10 to supplementary exploratory tables; otherwise does not generate the same simple-effect decomposition | Computes specified within-allocation and difference-of-differences contrasts regardless of omnibus p | Targeted equal-offer questions can be estimated directly, but this is a declared exception to the original reporting gate |
| Correction method | BH-FDR (`adjust="fdr"`) for emotion comparisons and interaction contrasts | Holm for specified log-odds contrasts | Analyst-selected change, not required by rarity or GLMM family |
| Correction families | Per experiment/model: ten emotion-main-effect pairs; simple emotion pairs within each offer stratum; ten expression-pair-by-offer interaction contrasts when generated | Twelve focal tests across both experiments; six secondary Disgust-Reward tests separately | Substantive change beyond the correction method; two separate Holm families do not control FWER across the entire paper at 5% |
| Other original comparison | Offer main-effect pair uses Tukey; only two offer levels | Not an additional new focal contrast family | Avoid describing every original comparison as FDR |
| Contrast test | Asymptotic Wald z in emmeans | Same | Aligned |
| Reported contrast intervals | For asymptotic output, formatter constructs `lower.CL`/`upper.CL` from b +/- 1.96 SE and uses them for OR CIs | Explicit unadjusted 95% Wald intervals using the normal quantile | Broadly aligned; neither should be called a Holm/BH-adjusted CI |
| Interval-column caveat | Original tables can contain both emmeans `asymp.LCL/UCL` and formatter `lower.CL/upper.CL`; they need not coincide | One explicit unadjusted interval pair for contrasts | Distinguish columns before quoting; this update does not repair or reinterpret old outputs |
| Scale | Choice ORs/log-odds; back-transformed cell means in figures | Also provides probability differences; exponentiated interaction is a ratio of ORs | Additional descriptive scale; not a replacement for designated model-scale tests |
| New sensitivity checks | Original project contains actor-sensitivity analyses separately | Adds actor intercept, every-participant deletion and parametric count simulations for new fair-only models | Supplemental checks, not part of the original five-step model-selection procedure |
| Analysis status | Mix of preregistered/exploratory components; infer status from registration and deviations | Entire new analysis is retrospective exploratory | A script label `Confirmatory` triggered by p<.05 does not establish preregistration |

## Why the new workflow was different

The assistant chose a transparent diagonal random structure because fair rejection
was sparse, and Holm to constrain false-positive conclusions across a focused set
of follow-up questions. These choices were documented before inspecting the new
focal contrasts, but they are **not empirically demonstrated to be necessary or
superior** to the manuscript workflow. The broad authorization to implement the
follow-up did not make these choices identical to the original methods.

In particular, successful convergence and exact numerical reproduction establish
that the selected procedure was executed consistently. They do not establish
that the diagonal model is preferable, that Holm is uniquely appropriate, or that
the new results are an exact-method replication of the original pipeline.

Consistency would normally favor applying the manuscript's selection strategy
and correction policy to the new estimand, allowing diagnostic-driven exceptions
only when demonstrated. That would mean rerunning selection on the fair-only data,
not copying a final random-effects formula estimated on all fair/unfair data.
Family definitions and the targeted-comparison exception must also be made explicit.
Any future change must retain the current version and be motivated by consistency,
not by which version produces more significant findings. **No such change was
made in this revision, following Jin's explicit instruction.**

## Source locators

Original manuscript methods: `E:/MyPaper_withARS/draft/professor_feedback_2026-09-28/manuscript_original_view.md`,
P065-P069 (paragraph identifiers). This is the previously structurally verified
original-view extraction of the professor-commented Word file; the methods were
cross-checked against the live code below rather than inferred from prose alone.

Original executable source: `C:/Code/UG_ERP_Project/Sta_Behaviour_E1_E2_Integrative.Rmd`:
configuration lines 73-84; GLMM gates around 203-272; blocked random-column
selection around 332-386; final optimizer comparison 391-415; CI formatter
689-723; buildmer alpha=.20 around 1717/1739; post-hoc gate/families 2055-2144.

New executable sources: `run_analysis.R`, functions `select_model`, `diagnose`,
`extract_contrasts`, and the final grouped Holm adjustment; `config.R` defines
family members, thresholds and explicit input/run identity.

The original E1 formal `GLMM_Rejection/posthoc_Main_Effect_emotion.md` was inspected
to confirm the distinction between `asymp.LCL/UCL` and the reported `lower.CL/upper.CL`
and OR interval columns. No new numeric result was inferred or recalculated here.

## Figure revision

`plot_paired.py` and `manuscript_style.json` reproduce the manuscript's expression
order, color assignments, within-participant connecting lines, consistent jitter,
endpoint dots, black-bordered model diamonds, confidence whiskers and neutral
reference line. The necessary changes are allocation labels (5:5/6:4), rows
(E1/E2 rather than choice/RT), the neutral reference (6:4), and slightly larger text.
All existing model estimates, intervals and p-values remain unchanged.
