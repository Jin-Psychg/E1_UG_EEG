"""Verify completed actor sensitivity outputs and write a coverage-aware report."""
from pathlib import Path
import hashlib
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

primary, run = (Path(x).resolve() for x in sys.argv[1:3])
if not (run / "STATUS.txt").read_text().startswith("COMPLETE"):
    raise RuntimeError("Run incomplete")
if (run / "REPORT.md").exists():
    raise FileExistsError("Preserve existing report")
checks = []


def check(name, passed):
    checks.append({"check": name, "passed": bool(passed)})
    if not passed:
        raise AssertionError(name)


def table(df):
    def fmt(x):
        return f"{x:.6g}" if isinstance(x, (float, np.floating)) else str(x).replace("|", "/")
    return "\n".join(["| " + " | ".join(df.columns) + " |", "| " + " | ".join(["---"] * len(df.columns)) + " |"] +
                     ["| " + " | ".join(map(fmt, row)) + " |" for row in df.itertuples(index=False, name=None)])


for row in pd.read_csv(run / "input_manifest.csv").itertuples():
    check(f"unchanged input {row.path}", hashlib.md5(Path(row.path).read_bytes()).hexdigest() == row.md5)
status = pd.read_csv(run / "selection_status.csv")
effects = pd.read_csv(run / "marginal_effects.csv")
omnibus = pd.read_csv(run / "omnibus.csv")
attempts = pd.read_csv(run / "fit_attempts.csv") if (run / "fit_attempts.csv").exists() else pd.DataFrame()
cov = pd.read_csv(run / "covariance_diagnostics.csv") if (run / "covariance_diagnostics.csv").exists() else pd.DataFrame()
check("both experiments resolved", set(status.experiment) == {"E1", "E2"} and len(status) == 2)
check("R arithmetic and structure checks", pd.read_csv(run / "verification_checks.csv").passed.all())
refits = pd.read_csv(run / "independent_refit_comparison.csv") if (run / "independent_refit_comparison.csv").exists() else pd.DataFrame()
expected_refits = int(status.actor_valid.sum() + status.matched_pair_valid.sum())
check("independent refit coverage", len(refits) == expected_refits * 4)
if len(refits):
    check("independent refits exact and valid", refits.valid.all() and (refits.maximum_absolute_difference == 0).all())
if (run / "correlation_selection.csv").exists():
    correlation = pd.read_csv(run / "correlation_selection_verified.csv")
    check("corrected correlation table arithmetic", np.allclose(correlation.p_value, chi2.sf(correlation.chisq, correlation.df_difference), atol=1e-12))
for key, e in effects.groupby(["experiment", "model", "effect_type"]):
    n = 10 if key[2] == "emotion_marginal" else 1
    check(f"{key} complete family", len(e) == n and not e.contrast.duplicated().any() and (e.family_n == n).all())
    p = 2 * norm.sf(np.abs(e.estimate / e.SE))
    order = np.argsort(p)
    adjusted = np.empty(n)
    adjusted[order] = np.minimum(1, np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1])
    check(f"{key} p and BH", np.allclose(e["p.value"], p, atol=1e-12) and np.allclose(e.p_BH, adjusted, atol=1e-12))
    check(f"{key} Wald and OR", np.allclose(e["asymp.LCL"], e.estimate - norm.ppf(.975) * e.SE) and
          np.allclose(e["asymp.UCL"], e.estimate + norm.ppf(.975) * e.SE) and np.allclose(e.OR, np.exp(e.estimate)) and
          np.allclose(e.OR_lower, np.exp(e["asymp.LCL"])) and np.allclose(e.OR_upper, np.exp(e["asymp.UCL"])))
check("omnibus chi-square arithmetic", np.allclose(omnibus["Pr(>Chisq)"], chi2.sf(omnibus.Chisq, omnibus.Df), atol=1e-12))
coverage, comparison, influence, selected_cov, omnibus_ranges, failure_details = [], [], [], [], [], []
report = ["# Actor-preserving structure sensitivity", "", "Retrospective exploratory analysis; human verification: false.",
          "The primary models and reporting gate are unchanged. This stage selects participant random structure with an actor intercept forced throughout, then compares that structure with and without actor. It does not isolate actor adjustment by comparing a newly selected model only with the original primary model.",
          "Tables below are rounded to six significant digits for display. The stat-audit CSV files retain literal source fields and exact row keys.", "", table(status[["experiment", "actor_valid", "matched_pair_valid", "reason"]])]
for s in status.itertuples():
    ex = s.experiment
    e = effects[effects.experiment == ex]
    a = attempts[attempts.experiment == ex] if len(attempts) else attempts
    expected = {"original_primary"}
    check(f"{ex} selected actor object", (run / f"{ex}_actor_selected.rds").exists() == s.actor_valid)
    available = a.groupby("model").valid.any() if len(a) else pd.Series(dtype=bool)
    trials = pd.read_csv(primary / f"{ex}_analysis_trials.csv")
    check(f"{ex} matched status", bool(available.get("matched_no_actor", False)) == s.matched_pair_valid)
    for label, valid in available.items():
        aa = a[a.model == label]
        check(f"{ex} {label} two optimizers", len(aa) == 2 and set(aa.optimizer) == {"BFGS", "nlminb"})
        check(f"{ex} {label} valid-only model", (run / f"{ex}_{label}.rds").exists() == valid)
        dd = trials[trials.participant_id_internal != label.removeprefix("actor_delete_")] if label.startswith("actor_delete_") else trials
        check(f"{ex} {label} unchanged target rows", (aa.n_trials == len(dd)).all() and (aa.n_participants == dd.participant_id_internal.nunique()).all())
        if not valid:
            for attempt in aa.itertuples():
                zz = cov[(cov.experiment == ex) & (cov.model == label) & (cov.optimizer == attempt.optimizer)]
                tiny = zz[zz.min_variance < 1e-6]
                failure_details.append({"experiment": ex, "model": label, "optimizer": attempt.optimizer,
                    "variance_below_gate": "; ".join(tiny.group), "diagnostics": attempt.diagnostics,
                    "error": attempt.error, "warnings": attempt.warnings})
    report += ["", f"## {ex}", "", f"Selection status: {s.reason}."]
    if "NaN" in (run / f"{ex}_selection.log").read_text(encoding="utf-8", errors="replace"):
        report += ["The selection log contains undefined (NaN) LRT entries. Do not translate software messages into a claim that every retained random term was supported by a valid significant LRT; selection includes numerical-feasibility constraints. See the log and the current-run interpretation note for the affected stage."]
    if s.actor_valid:
        report += ["", "Actor formula:", "```r", s.actor_formula, "```", "Matched no-actor formula:", "```r", s.no_actor_formula, "```"]
        opts = {"actor_selected": None}
        if s.matched_pair_valid:
            opts.update({label: a[(a.model == label) & a.valid].sort_values("AIC").iloc[0].optimizer for label in available[available].index})
        for label, opt in opts.items():
            z = cov[(cov.experiment == ex) & (cov.model == label)]
            if opt is not None:
                z = z[z.optimizer == opt]
            participant_cov = z[z.group.str.startswith("participant_id_internal")]
            selected_cov.append({"experiment": ex, "model": label, "optimizer": z.optimizer.iloc[0],
                                 "minimum_random_variance": z.min_variance.min(), "minimum_relative_eigenvalue": z.relative_min_eigenvalue.min(),
                                 "participant_full_block_diagonal_eigen_ratio": participant_cov.min_eigenvalue.min() / participant_cov.max_eigenvalue.max(),
                                 "maximum_absolute_correlation": z.max_abs_correlation.max()})
    if s.matched_pair_valid:
        expected |= {"actor_selected", "matched_no_actor"}
        ids = pd.read_csv(run / f"{ex}_deletion_roster.csv").participant_id_internal.tolist()
        deletions = available[available.index.str.startswith("actor_delete_")]
        check(f"{ex} all deletion targets attempted", set(deletions.index) == {f"actor_delete_{x}" for x in ids})
        expected |= set(deletions[deletions].index)
        failed = [x.removeprefix("actor_delete_") for x in deletions[~deletions].index]
        coverage.append({"experiment": ex, "valid_deletions": int(deletions.sum()), "attempted_deletions": len(ids), "failed_ids": "; ".join(failed)})
        oo = omnibus[(omnibus.experiment == ex) & omnibus.model.str.startswith("actor_delete_")]
        for term, rows in oo.groupby("term"):
            omnibus_ranges.append({"experiment": ex, "term": term, "valid_deletions": len(rows),
                                   "minimum_p": rows["Pr(>Chisq)"].min(), "maximum_p": rows["Pr(>Chisq)"].max()})
        full = e[e.model == "actor_selected"].set_index("contrast")
        no = e[e.model == "matched_no_actor"].set_index("contrast")
        prim = e[e.model == "original_primary"].set_index("contrast")
        for contrast, row in full.iterrows():
            nr, pr = no.loc[contrast], prim.loc[contrast]
            comparison.append({"experiment": ex, "contrast": contrast, "actor_beta": row.estimate, "actor_OR": row.OR,
                "actor_CI_lower": row.OR_lower, "actor_CI_upper": row.OR_upper, "actor_p_BH": row.p_BH,
                "matched_no_actor_beta": nr.estimate, "matched_no_actor_OR": nr.OR, "matched_no_actor_p_BH": nr.p_BH,
                "actor_minus_matched_beta": row.estimate - nr.estimate, "original_primary_OR": pr.OR,
                "original_primary_p_BH": pr.p_BH, "same_BH_class_primary": (row.p_BH < .05) == (pr.p_BH < .05),
                "same_direction_matched": np.sign(row.estimate) == np.sign(nr.estimate),
                "same_BH_class_matched": (row.p_BH < .05) == (nr.p_BH < .05)})
            dd = e[(e.contrast == contrast) & e.model.str.startswith("actor_delete_")]
            if len(dd):
                influence.append({"experiment": ex, "contrast": contrast, "valid_deletions": len(dd),
                    "min_beta": dd.estimate.min(), "max_beta": dd.estimate.max(), "min_p_BH": dd.p_BH.min(), "max_p_BH": dd.p_BH.max(),
                    "direction_changes": int((np.sign(dd.estimate) != np.sign(row.estimate)).sum()),
                    "BH_class_changes": int(((dd.p_BH < .05) != (row.p_BH < .05)).sum()),
                    "max_shift_in_full_model_SE": (abs(dd.estimate - row.estimate) / row.SE).max()})
        report += ["", f"Valid fixed-structure actor deletions: {int(deletions.sum())}/{len(ids)}. Failed IDs: {', '.join(failed) or 'none'}.",
                   "These are deletion results for the newly selected actor model, not repairs of the earlier primary-model deletion sequence."]
    else:
        report += ["No inferential actor/no-actor comparison or actor deletion sequence is available under the recorded plan."]
    check(f"{ex} only authorized valid inference", set(e.model) == expected)
    reference = pd.read_csv(primary / "emotion_main_contrasts.csv")
    reference = reference[reference.experiment == ex].set_index("contrast")
    actual = e[(e.model == "original_primary") & (e.effect_type == "emotion_marginal")].set_index("contrast").loc[reference.index]
    check(f"{ex} original primary unchanged", np.allclose(actual.estimate, reference.estimate, atol=1e-12) and np.allclose(actual.p_BH, reference.p_BH, atol=1e-12))
    allocation_ref = pd.read_csv(primary / f"{ex}_allocation_main.csv")
    allocation_now = e[(e.model == "original_primary") & (e.effect_type == "allocation_marginal")]
    check(f"{ex} original allocation unchanged", all(np.allclose(allocation_now[n], allocation_ref[n], atol=1e-12) for n in ["estimate", "SE", "p.value"]))
    omnibus_ref = pd.read_csv(primary / "new_omnibus.csv")
    omnibus_ref = omnibus_ref[omnibus_ref.experiment == ex].set_index("term")
    omnibus_now = omnibus[(omnibus.experiment == ex) & (omnibus.model == "original_primary")].set_index("term").loc[omnibus_ref.index]
    check(f"{ex} original omnibus unchanged", all(np.allclose(omnibus_now[n], omnibus_ref[n], atol=1e-12) for n in ["Chisq", "Df", "Pr(>Chisq)"]))
    om = omnibus[(omnibus.experiment == ex) & ~omnibus.model.str.startswith("actor_delete_")]
    report += ["", "Full-sample omnibus tests (Type III Wald):", "", table(om)]

for name, rows in [("deletion_coverage", coverage), ("matched_comparisons", comparison), ("deletion_influence", influence), ("deletion_omnibus_ranges", omnibus_ranges), ("selected_covariance_summary", selected_cov), ("failed_fit_details", failure_details)]:
    df = pd.DataFrame(rows)
    df.to_csv(run / f"{name}.csv", index=False)
    if len(df):
        report += ["", f"## {name.replace('_', ' ').title()}", "", table(df)]
report += ["", "## Interpretation boundaries", "",
    "- Actor intercepts account for baseline variation among these actor identities; no actor-specific emotion slopes were estimated.",
    "- All emotion pairs are marginal over 5:5 and 6:4 with equal log-odds weighting. They cannot establish a 5:5-specific expression effect.",
    "- The allocation contrast is 5:5 minus 6:4, averaged equally over the five emotions on the log-odds scale. Its OR is an odds ratio, not a probability ratio.",
    "- BH correction uses ten emotion pairs within each experiment/model; allocation has one comparison. Wald confidence intervals are unadjusted, so interval exclusion and adjusted significance can disagree.",
    "- Matched coefficient changes are descriptive. Different significance classifications are not themselves a statistical test of a difference.",
    "- Sensitivity tests do not override the original primary interaction gate: E1 decomposition remains supplementary; E2 remains audit-only. This is the project's convention, not a universal rule.",
    "- Non-significant interaction is not evidence of equivalence or the absence of an interaction.",
    "- Selection was data-driven, not preregistered. Standard Wald uncertainty does not incorporate random-structure selection uncertainty.",
    "- The raw legacy correlation-selection table mislabels candidate parameter count as df_diff. Use correlation_selection_verified.csv for the actual difference. Independently recomputed p-values agree; no selection decision changed.",
    "- Passing the original fit gate is not proof of scientific reliability or full covariance rank. Inspect covariance eigenvalues and deletion coverage.",
    "- Failed fits supply no evidence that effects are absent. Any deletion summaries apply only to valid fits.",
    "- These behavioral results do not distinguish status defense, intentions, or partner valuation mechanisms.",
    "", "## Reproducibility", "", "See ANALYSIS_PLAN.md, actor_engine.diff, engine_provenance.json, SESSION_INFO.txt, selection logs and stage RDS objects. Full-sample independent refits are documented in REFIT_VERIFICATION.txt. Numerical validation does not replace human scientific review."]
pd.DataFrame(checks).to_csv(run / "python_verification_checks.csv", index=False)
audit = pd.read_csv(run / "marginal_effects.csv", dtype=str, keep_default_na=False)
audit["source"] = str(run / "marginal_effects.csv")
audit["source_key"] = audit.experiment + "/" + audit.model + "/" + audit.contrast
audit["claim_scope"] = "Marginal contrast for the named valid model only; not a simple effect or mechanism test"
audit.to_csv(run / "stat_audit_marginal.csv", index=False)
audit_o = pd.read_csv(run / "omnibus.csv", dtype=str, keep_default_na=False)
audit_o["source"] = str(run / "omnibus.csv")
audit_o["source_key"] = audit_o.experiment + "/" + audit_o.model + "/" + audit_o.term
audit_o.to_csv(run / "stat_audit_omnibus.csv", index=False)
(run / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
print(f"Passed {len(checks)} Python checks; report written to {run / 'REPORT.md'}")
