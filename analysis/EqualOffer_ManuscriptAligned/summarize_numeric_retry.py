"""Verify and summarize bounded numerical retries without modifying prior runs."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm

previous, run = map(lambda s: Path(s).resolve(), sys.argv[1:3])
if not (run / "STATUS.txt").read_text().startswith("COMPLETE"):
    raise RuntimeError("Retry run is incomplete")
if (run / "REPORT.md").exists():
    raise FileExistsError("Preserve the existing report")
checks = []


def check(name, value):
    checks.append({"check": name, "passed": bool(value)})
    if not value:
        raise AssertionError(name)


def table(df):
    def fmt(v):
        return f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v).replace("|", "/")
    return "\n".join(["| " + " | ".join(df.columns) + " |",
                      "| " + " | ".join(["---"] * len(df.columns)) + " |"] +
                     ["| " + " | ".join(map(fmt, row)) + " |" for row in df.itertuples(index=False, name=None)])


for r in pd.read_csv(run / "input_manifest.csv").itertuples():
    check(f"unchanged source {r.path}", hashlib.md5(Path(r.path).read_bytes()).hexdigest() == r.md5)
attempts = pd.read_csv(run / "fit_attempts.csv")
old = pd.read_csv(previous / "fit_attempts.csv")
status = pd.read_csv(run / "retry_status.csv")
cov = pd.read_csv(run / "covariance_diagnostics.csv")
effects = pd.read_csv(run / "marginal_effects.csv")
omnibus = pd.read_csv(run / "omnibus.csv")
f = pd.read_csv(run / "formula_checks.csv")
check("formula text and row counts unchanged", f.expression_equal.all() and (f.expected_rows == f.actual_rows).all())
check("R arithmetic checks passed", pd.read_csv(run / "verification_checks.csv").passed.all())
expected = old.groupby(["experiment", "model"]).valid.any()
expected = set(expected[~expected].index)
check("all and only failed targets", set(zip(status.experiment, status.model)) == expected)
for key, a in attempts.groupby(["experiment", "model"]):
    check(f"{key} exactly two optimizers", len(a) == 2 and set(a.optimizer) == {"BFGS", "nlminb"})
    previous_rows = old[(old.experiment == key[0]) & (old.model == key[1])]
    check(f"{key} unchanged trial/participant counts", set(a.n_trials) == set(previous_rows.n_trials) and
          set(a.n_participants) == set(previous_rows.n_participants))
    s = status[(status.experiment == key[0]) & (status.model == key[1])].iloc[0]
    check(f"{key} correct recovery status", s.recovered_original_gate == a.valid.any())
    if s.recovered_original_gate:
        check(f"{key} selected minimum AIC valid fit", s.selected_optimizer == a[a.valid].sort_values("AIC").iloc[0].optimizer)
    check(f"{key} valid-only selected object", (run / f"{key[0]}_{key[1]}.rds").exists() == s.recovered_original_gate)
check("complete attempt count", len(attempts) == 2 * len(expected))
for key, e in effects.groupby(["experiment", "model", "effect_type"]):
    n = 10 if key[2] == "emotion_marginal" else 1
    check(f"{key} full comparison family", len(e) == n and not e.contrast.duplicated().any())
    p = 2 * norm.sf(np.abs(e.estimate / e.SE))
    order = np.argsort(p)
    adj = np.empty(n)
    adj[order] = np.minimum(1, np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1])
    check(f"{key} p and correction", np.allclose(e["p.value"], p, atol=1e-12) and np.allclose(e.p_BH, adj, atol=1e-12))
    check(f"{key} CI and OR", np.allclose(e["asymp.LCL"], e.estimate - norm.ppf(.975) * e.SE) and
          np.allclose(e["asymp.UCL"], e.estimate + norm.ppf(.975) * e.SE) and np.allclose(e.OR, np.exp(e.estimate)))
coverage = []
selected_cov = []
for ex in ["E1", "E2"]:
    prev = old[(old.experiment == ex) & old.model.str.startswith("delete_")].groupby("model").valid.any()
    s = status[(status.experiment == ex) & status.model.str.startswith("delete_")]
    recovery = int(s.recovered_original_gate.sum())
    actor = status[(status.experiment == ex) & (status.model == "actor_intercept")].iloc[0]
    coverage.append({"experiment": ex, "previous_valid_deletions": int(prev.sum()), "recovered_deletions": recovery,
                     "valid_deletions_after_retry": int(prev.sum()) + recovery, "total_deletions": len(prev),
                     "actor_recovered_original_gate": actor.recovered_original_gate})
    available = set(effects[effects.experiment == ex].model)
    recovered = set(status[(status.experiment == ex) & status.recovered_original_gate].model)
    check(f"{ex} valid-only inference rows", available == recovered | {"primary"})
    for row in status[(status.experiment == ex) & status.recovered_original_gate].itertuples():
        z = cov[(cov.experiment == ex) & (cov.model == row.model) & (cov.optimizer == row.selected_optimizer)]
        selected_cov.append({"experiment": ex, "model": row.model, "optimizer": row.selected_optimizer,
                             "min_relative_covariance_eigenvalue": z.relative_min_eigenvalue.min(),
                             "minimum_random_variance": z.min_variance.min(), "max_abs_correlation": z.max_abs_correlation.max()})
pd.DataFrame(coverage).to_csv(run / "coverage_after_retry.csv", index=False)
selected_cov = pd.DataFrame(selected_cov)
selected_cov.to_csv(run / "selected_covariance_summary.csv", index=False)
pd.DataFrame(checks).to_csv(run / "python_verification_checks.csv", index=False)
# Preserve literal source values in the result audit, with row keys.
audit = effects.copy()
audit["source"] = "marginal_effects.csv; experiment + model + contrast"
audit["scope"] = "Original-gate-valid numerical retry only; primary unchanged; unadjusted Wald CIs"
audit.to_csv(run / "stat_source_audit.csv", index=False)
parts = ["# Bounded numerical retry report", "", "## Material Passport", "",
         "- Origin: academic-research-suite / experiment-agent; run + validate",
         "- Date: 2026-09-29",
         "- Verification status: ANALYZED; separate exact-refit verification if available",
         "- Human verification: false",
         "- Version: numeric_retry_v3", "", "## Scope", "",
         "Only previously failed targets were retried. BFGS maxit=2000; nlminb iter.max=2000, eval.max=4000. "
         "Default starts, original formulas, exclusions, factors, numerical acceptance criteria and primary "
         "reporting routes were preserved. This is a computational sensitivity check, not a new preregistered "
         "analysis. No random structure was simplified and no primary model was replaced.", "",
         "## Coverage", "", table(pd.DataFrame(coverage)), "", "## Per-target result", "",
         table(status[["experiment", "model", "recovered_original_gate", "selected_optimizer"]]), "",
         "## Covariance diagnostics for recovered fits", "",
         table(selected_cov) if not selected_cov.empty else "No recovered fits.", "",
         "Passing the original gate is not proof of full rank or inferential stability. Its variance check "
         "uses individual diagonal variances; small covariance eigen-directions may remain. These are "
         "reported without inventing an additional universal pass/fail threshold. Actor intercepts address "
         "baseline identity differences only. Do not equate failure to recover an actor fit with zero actor variability.", "",
         "## Optimizer agreement", "", table(pd.read_csv(run / "optimizer_comparison.csv")), "",
         "Objective/parameter comparisons include diagnostic-only invalid attempts and are not inferential "
         "findings. Similar objectives do not prove that an ill-conditioned covariance is well estimated.", "",
         "## Omnibus output for interpretable candidates", "", table(omnibus), "",
         "All 11 marginal contrasts per interpretable model, with ORs, unadjusted Wald intervals and correction "
         "families, are in marginal_effects.csv. Sensitivity p-values cannot unlock the primary decomposition gate.", "",
         "## Validation and limitations", "",
         f"All {len(checks)} Python checks and saved R input/matrix checks passed. See formula_checks.csv, "
         "input_manifest.csv and verification files. Any exact recovered-model refit is documented separately "
         "in REFIT_VERIFICATION.txt; do not assume the full retry sequence was independently repeated.", "",
         "The prior 11-category interpretation audit continues to apply: no individual-level or causal "
         "mechanism claim, no null/equivalence inference from failed or nonsignificant fits, no primary "
         "participant removal, no significance-driven selection, no claim that successful-fit subsets "
         "cover failures. Recruitment/generalizability and serial dependence were not assessed.", "",
         "Development v1/v2 stopped during formula-object attribute checks. v2 preserves identical expected "
         "and actual formula texts/trial counts. v3 verifies serialized formulas directly. Both earlier "
         "runs are retained and are not inferential sources.", ""]
(run / "REPORT.md").write_text("\n".join(parts), encoding="utf-8")
(run / "QA.json").write_text(json.dumps({"checks_passed": len(checks), "source_runs_unchanged": True}, indent=2))
print(json.dumps({"checks_passed": len(checks), "coverage": coverage}, indent=2, default=lambda x: bool(x)))
