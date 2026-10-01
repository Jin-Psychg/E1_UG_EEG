"""Validate and summarize the fixed-plan robustness run; never refit or overwrite."""
from pathlib import Path
import hashlib
from importlib.metadata import version
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

run = Path(sys.argv[1]).resolve()
if not (run / "STATUS.txt").read_text().startswith("COMPLETE"):
    raise RuntimeError("Fit/simulation run is incomplete")
if (run / "ROBUSTNESS_REPORT.md").exists():
    raise FileExistsError("Preserve the existing report")
checks = []


def check(name, passed):
    checks.append({"check": name, "passed": bool(passed)})
    if not passed:
        raise AssertionError(name)


def table(df):
    # No optional tabulate dependency. Only display values are rounded.
    def cell(x):
        if isinstance(x, (float, np.floating)):
            return f"{x:.6g}"
        return str(x).replace("|", "/")
    rows = [[cell(v) for v in row] for row in df.itertuples(index=False, name=None)]
    return "\n".join(["| " + " | ".join(df.columns) + " |",
                      "| " + " | ".join(["---"] * len(df.columns)) + " |"] +
                     ["| " + " | ".join(row) + " |" for row in rows])


attempts = pd.read_csv(run / "fit_attempts.csv")
effects = pd.read_csv(run / "marginal_effects.csv")
omnibus = pd.read_csv(run / "omnibus.csv")
simulation = pd.read_csv(run / "simulation_summary.csv")
routes = omnibus[(omnibus.model == "primary") & (omnibus.term == "emotion:allocation")][
    ["experiment", "Pr(>Chisq)"]].copy()
routes["decomposition_route"] = routes["Pr(>Chisq)"].map(
    lambda p: "main" if p < .05 else "supplement" if p < .10 else "audit-only")
routes.to_csv(run / "unchanged_primary_routes.csv", index=False)
check("R checks all pass", pd.read_csv(run / "verification_checks.csv").passed.all())
for row in pd.read_csv(run / "input_manifest.csv").itertuples():
    check(f"input hash: {row.path}", hashlib.md5(Path(row.path).read_bytes()).hexdigest() == row.md5)
for key, g in effects.groupby(["experiment", "model", "effect_type"]):
    n = 10 if key[2] == "emotion_marginal" else 1
    check(f"complete family: {key}", len(g) == n and not g.contrast.duplicated().any())
    p = 2 * norm.sf(abs(g.estimate.to_numpy() / g.SE.to_numpy()))
    order = np.argsort(p)
    corrected = np.empty(n)
    corrected[order] = np.minimum(1, np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1])
    check(f"p and BH: {key}", np.allclose(p, g["p.value"], atol=1e-12, rtol=1e-10) and
          np.allclose(corrected, g.p_BH, atol=1e-12, rtol=1e-10))
    check(f"Wald CI and OR: {key}",
          np.allclose(g.estimate - norm.ppf(.975) * g.SE, g["asymp.LCL"]) and
          np.allclose(g.estimate + norm.ppf(.975) * g.SE, g["asymp.UCL"]) and
          np.allclose(np.exp(g.estimate), g.OR))

summaries, status, actor_comparisons, selected = [], [], [], []
for ex in ["E1", "E2"]:
    a = attempts[attempts.experiment == ex]
    roster = pd.read_csv(run / f"{ex}_deletion_roster.csv").participant_id_internal
    expected = {f"delete_{pid}" for pid in roster} | {"actor_intercept"}
    check(f"{ex} complete attempts", set(a.model) == expected and len(a) == 2 * len(expected))
    for model, g in a.groupby("model"):
        check(f"{ex} {model} two optimizers", set(g.optimizer) == {"BFGS", "nlminb"} and len(g) == 2)
        ok = g[g.valid]
        exists = (run / f"{ex}_{model}.rds").exists()
        check(f"{ex} {model} valid-only saved", exists == (not ok.empty))
        if not ok.empty:
            best = ok.loc[ok.AIC.idxmin()]
            selected.append({"experiment": ex, "model": model, "optimizer": best.optimizer, "AIC": best.AIC})
    valid_models = set(a[a.valid].model)
    valid_deletions = {m for m in valid_models if m.startswith("delete_")}
    invalid = sorted(expected - valid_models - {"actor_intercept"})
    status.append({"experiment": ex, "valid_deletions": len(valid_deletions), "total_deletions": len(roster),
                   "invalid_deletions": ", ".join(invalid) or "none",
                   "actor_valid": "actor_intercept" in valid_models})
    e = effects[effects.experiment == ex]
    check(f"{ex} no invalid estimates", set(e.model) == valid_models | {"primary"})
    primary = e[e.model == "primary"].set_index("contrast")
    for contrast, row in primary.iterrows():
        loo = e[(e.model.isin(valid_deletions)) & (e.contrast == contrast)]
        if loo.empty:
            continue
        shift = (loo.estimate - row.estimate) / row.SE
        most = shift.abs().idxmax()
        summaries.append({"experiment": ex, "contrast": contrast, "primary_beta": row.estimate,
                          "primary_SE": row.SE, "primary_CI_low": row["asymp.LCL"],
                          "primary_CI_high": row["asymp.UCL"], "primary_p_BH": row.p_BH,
                          "valid_deletions": len(loo), "beta_min": loo.estimate.min(), "beta_max": loo.estimate.max(),
                          "max_abs_change_in_primary_SE": shift.abs().max(), "most_influential_model": loo.loc[most, "model"],
                          "direction_changes": int((np.sign(loo.estimate) != np.sign(row.estimate)).sum()),
                          "p_BH_min": loo.p_BH.min(), "p_BH_max": loo.p_BH.max()})
    actor = e[e.model == "actor_intercept"].set_index("contrast")
    if not actor.empty:
        z = actor[["estimate", "asymp.LCL", "asymp.UCL", "OR", "p_BH"]].join(
            primary[["estimate", "SE", "p_BH"]], rsuffix="_primary")
        z["change_in_primary_SE"] = (z.estimate - z.estimate_primary) / z.SE
        z["experiment"] = ex
        actor_comparisons.append(z.reset_index())
    draws = pd.read_csv(run / f"{ex}_simulation_draws.csv")
    check(f"{ex} simulation count", len(draws) == 1000 and list(draws.simulation) == list(range(1, 1001)))
    check(f"{ex} simulated total", np.array_equal(draws.filter(like="cell_").sum(axis=1), draws.total_rejections))
    for row in simulation[simulation.experiment == ex].itertuples():
        q = draws[row.metric].quantile([.025, .5, .975]).to_numpy()
        check(f"{ex} {row.metric} quantiles", np.allclose(q, [row.lower, row.median, row.upper]))

summary = pd.DataFrame(summaries)
summary.to_csv(run / "influence_summary.csv", index=False)
pd.DataFrame(status).to_csv(run / "sensitivity_status.csv", index=False)
pd.DataFrame(selected).to_csv(run / "selected_optimizers.csv", index=False)
actor = pd.concat(actor_comparisons, ignore_index=True) if actor_comparisons else pd.DataFrame(
    columns=["contrast", "estimate", "asymp.LCL", "asymp.UCL", "OR", "p_BH",
             "estimate_primary", "SE", "p_BH_primary", "change_in_primary_SE", "experiment"])
actor.to_csv(run / "actor_comparison.csv", index=False)
audit = effects.copy()
audit["source"] = "marginal_effects.csv; experiment + model + contrast identify each row"
audit["claim_scope"] = "Sensitivity only; primary reporting gate unchanged; unadjusted Wald CI"
audit.to_csv(run / "stat_source_audit.csv", index=False)
pd.DataFrame(checks).to_csv(run / "python_verification_checks.csv", index=False)

report = ["# Fair-only selected-model robustness checks", "", "## Material Passport", "",
          "- Origin Skill: academic-research-suite / experiment-agent",
          "- Origin Mode: run + validate",
          "- Origin Date: 2026-09-28",
          "- Verification Status: ANALYZED (independent refit results recorded separately)",
          "- Version Label: selected_model_robustness_v1",
          "- Overall Confidence: CAUTION (incomplete valid-fit coverage; exploratory selection)", "",
          "Source: EqualOffer_Aligned_20260928_v2, canonical trial files and recorded input hashes. "
          "Status: ANALYZED; numerical reconstruction verified. Independent sensitivity refit status, "
          "if performed, is in REFIT_VERIFICATION.txt. Human verification: false. "
          "All analyses are retrospective exploratory.", "", "## Scope and safeguards", "",
          "The selected primary models, data exclusions, BH families and primary reporting gate remain unchanged. "
          "Sensitivity p-values do not unlock a decomposition. Emotion contrasts average equally across "
          "5:5 and 6:4 on the log-odds scale, retaining their interaction in every fitted model. "
          "They are not tests restricted to 5:5. The primary reporting routes are shown below and cannot "
          "be changed by sensitivity results. See code/ROBUSTNESS_PLAN.md for the pre-execution plan.", "",
          table(routes), "",
          "## Fit coverage", "", table(pd.DataFrame(status)), "",
          "Every deletion was attempted with both original optimizers. Invalid fits contribute no effect estimates. "
          "Counts of successful deletion fits are not evidence that failed fits would agree. "
          "No participant was removed from the primary analysis. All attempt warnings and weak covariance "
          "directions are retained in fit_attempts.csv and covariance_diagnostics.csv. "
          "The original glmmTMB gate checks optimizer success, positive-definite Hessian, finite SE/likelihood "
          "and individual random variances >=1e-6. Its variance-boundary check is not a full covariance-rank "
          "test; weak eigen-directions may remain even among accepted fits. Passing this gate is not a "
          "blanket certificate of model stability.", "",
          "## Participant influence", "",
          "All contrasts are shown, including nonsignificant ones. Changes are relative to the full-sample "
          "estimate, scaled by its SE; there is no universal stability cutoff. A change in significance "
          "is not a formal test of a difference between fits. Intervals and ORs for every valid fit are "
          "in marginal_effects.csv. Display tables use six significant digits; CSV files retain numerical precision.", "",
          table(summary[["experiment", "contrast", "primary_beta", "beta_min", "beta_max",
                         "max_abs_change_in_primary_SE", "most_influential_model", "direction_changes"]]), "",
          "## Actor intercept sensitivity", "",
          table(actor.drop(columns=["SE"], errors="ignore")) if not actor.empty else "No actor fit passed the numerical gate.", "",
          "An actor random intercept addresses baseline face-identity heterogeneity only. It does not test "
          "actor-by-emotion slopes or establish stimulus-population generalizability. Invalid or boundary "
          "fits cannot be interpreted as evidence of no actor effect.", "", "## Simulation adequacy", "",
          table(simulation), "",
          "Each primary fit generated 1,000 data sets with new random effects, holding fitted parameters fixed. "
          "These are pointwise descriptive simulation envelopes, not posterior predictive intervals or formal "
          "multiplicity-adjusted tests. Being inside an envelope is not proof of model adequacy; being outside "
          "identifies a feature for scrutiny. Parameter uncertainty, selection uncertainty and serial trial "
          "dependence were not evaluated. No extra-binomial dispersion claim is made from Bernoulli residuals.", "",
          "## Interpretation audit (11/11 categories considered)", "",
          "1. Simpson's paradox: no pooling of the allocation factor; inspect the retained interaction. "
          "No comprehensive subgroup reversal scan was performed.",
          "2. Ecological fallacy: trial-level mixed models; no claim that every individual shows a group effect.",
          "3. Berkson/selection bias: original eligibility/exclusions retained; representativeness not established.",
          "4. Collider bias: no outcome-derived covariate added; this is not a causal-DAG audit.",
          "5. Base-rate neglect: sparse rejection is assessed using counts and zero-rejection participants; no diagnostic-accuracy claim.",
          "6. Regression to the mean: no extreme participant is selected for primary removal; not a pre/post efficacy analysis.",
          "7. Survivorship bias: all primary participants retained; no assessment of recruitment nonresponse.",
          "8. Look-elsewhere effect: all planned deletion attempts/contrasts recorded, including failures; no discoveries from sensitivity p-values.",
          "9. Forking paths: prior analyses and choices remain retrospective; the plan fixes only this check sequence and cannot preregister past work.",
          "10. Causal overreach: no evidence identifying status defense, WTR/AV or intention mechanisms.",
          "11. Reverse causality: these checks do not test the direction of unmeasured psychological pathways.", "",
          "## Verification and limits", "",
          f"Passed {len(checks)} Python checks plus all saved R direct-matrix/input checks. "
          "Computational verification does not validate all distributional assumptions or turn exploratory "
          "analyses into confirmatory evidence. Read sensitivity_status.csv before making robustness claims.", ""]
(run / "ROBUSTNESS_REPORT.md").write_text("\n".join(report), encoding="utf-8")
(run / "SUMMARY_QA.json").write_text(json.dumps({"checks_passed": len(checks), "human_verified": False,
    "python": sys.version, "packages": {p: version(p) for p in ["numpy", "pandas", "scipy"]}}, indent=2))
print(json.dumps({"checks_passed": len(checks), "fit_status": status}, indent=2))
