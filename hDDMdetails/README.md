> Final-analysis update: these notebooks are historical reference copies.
> The selected final notebooks, configurations and follow-up scripts are
> distributed in the [final-analysis code addendum](../FINAL_ANALYSIS_ADDENDUM.md).
> See that page for download instructions and remaining environment checks.

# hDDMdetails/ — HDDM reference copies

Reference copies of the two final Hierarchical Drift Diffusion Model (HDDM)
notebooks, one per experiment:

- `Final_vazt-Unfair_Grouplevel_InformativePr_E1.ipynb`
- `Final_vazt-Unfair_Grouplevel_InformativePr_E2.ipynb`

These are **reference copies only**. The canonical HDDM analysis — with its own
environment, dependencies, and generated artifacts (`hddm_config.py`,
`hddm_data_unfair.csv`, manifests, audit CSVs) — is maintained in a **separate
parallel project tree**. Consult that project for the runnable version and the
environment specification. Do not treat this folder as a runnable sub-project.

## What the notebooks do

Both notebooks run the same 6-step pipeline on their experiment's data:

1. Data preprocessing & exploratory behavioral analysis
2. 2a Global configuration (cryptographic lineage / SSOT) · 2b Hierarchical
   Bayesian model specification & MCMC estimation
3. Posterior trace extraction & predictive simulation (PPC)
4. Convergence diagnostics & model selection
5. Statistical inference & publication-ready visualization
6. Data-informed group-level parameter recovery

## Model

- 4-parameter drift-diffusion model: drift rate `v`, boundary separation `a`,
  non-decision time `t`, starting-point bias `z` (explicitly included).
- **No** across-trial variability parameters (`sv`, `st`, `sz`) — per Lerche &
  Voss (2016) and Boehm et al. (2018) for per-condition trial counts of ~30–40.
- Fit with `HDDMRegressor`: every DDM parameter is regressed on emotion, with the
  treatment-contrast coefficients estimated at the **group level**
  (`group_only_regressors=True`); subject-level variability is carried by the
  hierarchical DDM base parameters. Informative priors (`use_informative_priors=True`).
- Software: dockerHDDM / HDDM (Wiecki, Sofer & Frank, 2013).

## Data scope & coding

- **Unfair offers only** (the HDDM-parallel stratum; matches `LMM_RT_unfair` in
  `Sta_Behaviour_E1_E2_Integrative.Rmd`).
- Trial inclusion: `reaction != 0`, RT in [300, 3000] ms.
- Response → boundary coding: Accept (`reaction == 1`) → upper boundary (1);
  Reject (`reaction == 2`) → lower boundary (0).
- Emotion recode: `enj` → `rew` (reward smile).

---
Basis: written from the two notebooks' own cells (markdown + Step 2a/2b code),
verified for both E1 and E2. The "separate parallel project" statement follows the
project README and CLAUDE.md; that project's path is not recorded in this repo.
