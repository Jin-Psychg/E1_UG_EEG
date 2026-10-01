"""Compare verified aligned fits with the preserved Holm analysis; render the paired figure."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import numpy as np
import pandas as pd


def holm(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    result = np.empty_like(values)
    result[order] = np.minimum(1, np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    run = args.run.resolve(strict=True)
    root = next(p for p in run.parents if (p / "UG_ERP_Project.Rproj").exists())
    if (run / "STATUS.txt").read_text().splitlines()[0] != "COMPLETE" or not (run / "VERIFICATION.txt").exists():
        raise ValueError("Independent verification must finish first")
    if (run / "REPORT.md").exists():
        raise FileExistsError("Report already exists; preserve the completed output")
    old_run = root / "results/EqualOffer_replication_20260928_224854"
    old = pd.read_csv(old_run / "new_contrasts.csv", float_precision="round_trip")
    old = old[(old.model == "primary") & (old.scale == "log_odds")]
    new = pd.read_csv(run / "new_contrasts.csv", float_precision="round_trip")
    rows = []
    for _, o in old.iterrows():
        pair, context = o.contrast.split(" | ")
        a, b = pair.split(" - ")
        match = new[(new.experiment == o.experiment) & (new.context == context) & new.pair.isin([pair, f"{b} - {a}"])]
        if len(match) != 1:
            raise ValueError(f"Nonunique contrast: {o.experiment} {o.contrast}")
        n = match.iloc[0]
        sign = 1 if n.pair == pair else -1
        estimate = sign * n.estimate
        lower, upper = sorted([sign * n["asymp.LCL"], sign * n["asymp.UCL"]])
        rows.append(dict(experiment=o.experiment, contrast=o.contrast, family=o.family,
            old_estimate=o.estimate, old_SE=o.SE, old_p_raw=o["p.value"], old_p_Holm=o.p_holm,
            new_estimate=estimate, new_SE=n.SE, new_lower=lower, new_upper=upper,
            new_OR_or_ratio=np.exp(estimate), new_OR_lower=np.exp(lower), new_OR_upper=np.exp(upper),
            new_p_raw=n["p.value"], new_p_BH=n.p_BH, original_omnibus_route=n.original_omnibus_route))
    comparison = pd.DataFrame(rows)
    for family, index in comparison.groupby("family").groups.items():
        if len(index) != (12 if family == "focal" else 6):
            raise ValueError("Legacy comparison family is incomplete")
        comparison.loc[index, "new_model_p_legacy_Holm"] = holm(comparison.loc[index, "new_p_raw"])
    comparison.to_csv(run / "version_comparison.csv", index=False)
    # Exact-value source audit: no rounded number becomes the source of a claim.
    audit = comparison.copy()
    audit["old_source"] = str(old_run / "new_contrasts.csv")
    audit["new_source"] = str(run / "new_contrasts.csv")
    audit["location"] = "experiment + contrast; reversed expression pairs sign-aligned explicitly"
    audit.to_csv(run / "stat_source_audit.csv", index=False)
    diag = pd.read_csv(run / "aligned_diagnostics.csv")
    omnibus = pd.read_csv(run / "new_omnibus.csv", dtype=str)
    covariance = pd.read_csv(run / "covariance_diagnostics.csv", float_precision="round_trip")
    if not (run / "paired_figures").exists():
        subprocess.run([sys.executable, str(Path(__file__).with_name("plot_paired.py")), "--run", str(run)], check=True)
    else:
        import hashlib
        receipt = json.loads((run / "paired_figures/figure_provenance.json").read_text())
        for path, expected in receipt["input_hashes_after"].items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError("Existing figure inputs have changed")
    # Matched plot coordinates must be identical to the approved earlier visualization.
    old_points = pd.read_csv(old_run / "paired_figures/plotted_participant_coordinates.csv")
    new_points = pd.read_csv(run / "paired_figures/plotted_participant_coordinates.csv")
    pd.testing.assert_frame_equal(old_points, new_points, check_exact=True)
    readme = ["# Manuscript-aligned equal-offer analysis", "", "## Material Passport", "",
        "- Origin Skill: academic-research-suite / experiment-agent", "- Origin Mode: run + validate",
        "- Verification Status: VERIFIED (selected-fit reproduction and arithmetic only; full selection not rerun)",
        "- Human verification: false", "- Version: manuscript-aligned v2", "",
        "## Scope and interpretation", "",
        "This retrospective exploratory revision aligns the model-selection and multiplicity policy with the manuscript. "
        "The previous Holm results were known before this revision; it is not preregistered or an independent replication. "
        "The original analyses, data and Holm outputs are preserved. No mechanism or theory ranking is identified.", "",
        "The exact original five-step selection functions are reused, with fixed emotion-by-allocation effects. "
        "BH correction covers all ten expression pairs within each experiment/allocation, and a separate ten-test "
        "interaction family per experiment. Main-expression comparisons use another ten-test BH family. "
        "The allocation main comparison uses Tukey (one comparison). All intervals are unadjusted 95% Wald intervals.", "",
        "Targeted-follow-up exception: within-allocation and interaction contrasts are exported regardless of the "
        "omnibus p value. The original main/supplement/not-triggered routing is recorded but does not make this analysis confirmatory.", "",
        "## Selected structures", ""]
    for _, d in diag.iterrows():
        readme.extend([f"- {d.experiment}: `{d.formula}`. Optimizer: {d.optimizer}. {d.reason}"])
    readme.extend(["", "## Omnibus tests", "", "Exact stored values; Type III Wald chi-square.", "",
        "| Experiment | Term | Chi-square | df | p |", "|---|---|---|---|---|"])
    for _, r in omnibus.iterrows():
        readme.append(f"| {r.experiment} | {r.term} | {r.Chisq} | {r.Df} | {r['Pr(>Chisq)']} |")
    readme.extend(["", "## Previously specified contrasts", "",
        "Display values below are rounded (OR/CI: 3 significant digits; p: 4). Exact values and sign alignment are in "
        "version_comparison.csv and stat_source_audit.csv. Within-allocation effects are ORs; interaction effects are ratios of ORs.", "",
        "| Experiment | Contrast | New OR or ratio [95% CI] | Old model / old Holm p | New model / old Holm p | New model / new BH p |",
        "|---|---|---|---|---|---|"])
    for _, r in comparison.iterrows():
        name = r.contrast.replace("|", "/")
        readme.append(f"| {r.experiment} | {name} | {r.new_OR_or_ratio:.3g} [{r.new_OR_lower:.3g}, {r.new_OR_upper:.3g}] | {r.old_p_Holm:.4g} | {r.new_model_p_legacy_Holm:.4g} | {r.new_p_BH:.4g} |")
    readme.extend(["", "The middle p column is a diagnostic bridge only: it applies the OLD 12/6-test Holm families "
        "to the NEW model estimates, isolating model changes from the combined correction/family change. It is not another primary analysis. "
        "All 60 new pairwise tests, including non-highlighted pairs, are available in new_contrasts.csv.", "",
        "## Model diagnostics and limitations", "",
        "Both final fits passed the manuscript's convergence, Hessian, variance and finite-SE/likelihood gate. "
        "Failed candidate fits are recorded in selection logs and excluded by the original procedure. "
        "Passing this gate is not proof of adequate covariance dimensionality or complete robustness.", ""])
    for ex in ("E1", "E2"):
        sel = json.loads((run / f"{ex}_aligned_selection.json").read_text())
        pca = sel["step5_pca"]
        eff = next(i + 1 for i, r in enumerate(pca) if r["cumulative_variance"] >= .999)
        z = covariance[covariance.experiment == ex]
        readme.append(f"- {ex}: {eff}/{len(pca)} dimensions explain 99.9% of random-effect variance; "
            f"minimum within-block eigenvalue {z.min_eigenvalue.min()!r}; maximum absolute within-block correlation {z.max_absolute_correlation.max()!r}. "
            "A weak direction warrants caution and is not concealed by a successful optimizer flag.")
    readme.extend(["", "The old actor-intercept, participant-deletion and simulation checks belong to the old models. "
        "They have not been rerun for these newly selected structures and cannot be presented as validation of them. "
        "Wald uncertainty conditions on the selected structure and omits selection uncertainty. "
        "Large ORs can coexist with low absolute rejection probabilities and wide uncertainty.", "",
        "## Figure", "", "![Paired figure](paired_figures/equal_offers_paired.png)", "",
        (run / "paired_figures/caption.txt").read_text().strip(), "",
        "## Verification", "", (run / "VERIFICATION.txt").read_text().split("R version")[0].strip(), "",
        "Participant plot coordinates match the previously approved figure exactly. Only model diamonds/intervals and the model reference line change. "
        "See covariance_diagnostics.csv, independent_refit_comparison.csv, verification_checks.csv and figure_provenance.json.", "",
        "## Fallacy scan (11/11 considered)", "",
        "| Check | Assessment |", "|---|---|"])
    for name, detail in [
        ("Simpson's paradox", "Experiments and allocations are reported separately; no pooled direction claim."),
        ("Ecological fallacy", "Trial-level models and participant points are distinguished; no individual trait inference."),
        ("Berkson selection", "Original RT/response exclusions remain; scope is retained preprocessing trials."),
        ("Collider bias", "No reaction-conditioned subset or added post-treatment covariate; RT filtering remains a limitation."),
        ("Base-rate neglect", "Observed rates and conditional probabilities accompany ORs."),
        ("Regression to mean", "No extreme-participant selection or pre-post improvement claim."),
        ("Survivorship", "All participants including zero rejecters retained; earlier preprocessing attrition not re-audited."),
        ("Look-elsewhere", "All pairs exported and explicit BH families; not whole-paper error control."),
        ("Forking paths", "Old and new versions retained; revision follows observed results and stays exploratory."),
        ("Causation", "Expression/offer contrasts do not identify measured or unmeasured psychological mediation."),
        ("Reverse causality", "No inferred WTR/status/intentional mechanism direction claimed.")]:
        readme.append(f"| {name} | {detail} |")
    (run / "REPORT.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    (run / "report_QA.json").write_text(json.dumps({"participant_coordinates_exact_match": True,
        "comparison_rows": len(comparison), "all_new_contrasts": len(new), "statistical_files_read_only": True}, indent=2), encoding="utf-8")
    print(run / "REPORT.md")


if __name__ == "__main__":
    main()
