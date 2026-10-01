"""Apply the manuscript's original omnibus reporting gate to verified saved results.

No models, p-values, CIs or figures are recomputed. Prior reports remain audit records.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def route(p):
    return "main" if p < .05 else "supplement" if p < .10 else "audit_only"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def table(rows, columns):
    return ["| " + " | ".join(label for _, label in columns) + " |",
            "| " + " | ".join("---" for _ in columns) + " |"] + [
        "| " + " | ".join(str(row[key]).replace("|", "/") for key, _ in columns) + " |" for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    run = parser.parse_args().run.resolve(strict=True)
    dest = run / "standard_reporting"
    if dest.exists():
        raise FileExistsError("Preserve existing report; choose a fresh run for another generation")
    if (run / "STATUS.txt").read_text().splitlines()[0] != "COMPLETE" or not (run / "VERIFICATION.txt").exists():
        raise ValueError("Verified completed fits required")
    inputs = [run / name for name in ["new_omnibus.csv", "emotion_main_contrasts.csv", "new_contrasts.csv",
              "E1_allocation_main.csv", "E2_allocation_main.csv", "aligned_diagnostics.csv", "covariance_diagnostics.csv"]]
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    omnibus = read(run / "new_omnibus.csv")
    emotions = read(run / "emotion_main_contrasts.csv")
    contrasts = read(run / "new_contrasts.csv")
    gates = {}
    for ex in ("E1", "E2"):
        values = [r for r in omnibus if r["experiment"] == ex and r["term"] == "emotion:allocation"]
        if len(values) != 1:
            raise ValueError("Unique interaction test required")
        gates[ex] = route(float(values[0]["Pr(>Chisq)"]))
    dest.mkdir()
    header = ["# Standard manuscript-rule report", "", "Human verification: false.", "",
        "This report supersedes the presentation hierarchy in ../REPORT.md. The author withdrew the "
        "targeted-follow-up exception after reviewing the methods. This is a retrospective reporting-policy "
        "amendment, not a preregistration. The fitted models and every statistical value are unchanged.", "",
        "The design remains emotion (five levels) by offer ratio (5:5, 6:4). Main emotion comparisons "
        "are equally weighted across the two offer-ratio levels on the log-odds scale; trials were not "
        "collapsed and the interaction was not removed from the model.", "",
        "Reporting gate copied from the manuscript: interaction p<.05 routes decomposition to the main "
        "results; .05<=p<.10 routes it to supplementary exploratory results; p>=.10 leaves it audit-only. "
        "The .10 rule is a project convention, not a universal statistical requirement. No p value "
        "converts this retrospective analysis into a confirmatory test.", "",
        "All numbers below are literal CSV values. Log-odds estimates and unadjusted 95% Wald CIs are "
        "reported with adjusted p values. CI exclusion of zero does not replace the multiplicity-adjusted test.", "",
        "## Omnibus tests", ""]
    header += table([r for r in omnibus if r["term"] != "(Intercept)"],
        [("experiment","Experiment"),("term","Term"),("Chisq","Wald chi-square"),("Df","df"),("Pr(>Chisq)","p")])
    header += ["", "## Emotion marginal comparisons", "",
        "All ten pairs per experiment are shown, including nonsignificant pairs. BH-FDR is applied within "
        "each experiment's ten-test family. Positive estimates mean higher rejection log odds for the first expression.", ""]
    cols = [("experiment","Experiment"),("contrast","Contrast"),("estimate","Log-odds difference"),
            ("asymp.LCL","95% CI lower"),("asymp.UCL","95% CI upper"),("p_BH","BH p")]
    header += table(emotions,cols)
    allocation = []
    for ex in ("E1","E2"):
        allocation += [dict(r, experiment=ex) for r in read(run/f"{ex}_allocation_main.csv")]
    header += ["", "## Offer-ratio marginal comparison", "",
        "Equally weighted across five expressions on the log-odds scale. Tukey for one pair per experiment.", ""]
    header += table(allocation,cols[:-1]+[("p.value","p")])
    routed = []
    supp = ["# Supplementary exploratory decomposition", "",
        "Included only where .05<=interaction p<.10 under the original project rule. This does not "
        "establish an interaction or a trend. Comparisons use the unchanged BH families; CIs are unadjusted.", ""]
    for ex, gate in gates.items():
        subset = [r for r in contrasts if r["experiment"] == ex]
        if len(subset) != 30:
            raise ValueError("Expected all 30 simple and interaction comparisons per experiment")
        routed += [dict(r, reporting_route=gate) for r in subset]
        if gate == "main":
            header += ["", f"## {ex} interaction decomposition", ""] + table(subset,cols)
        elif gate == "supplement":
            supp += [f"## {ex}", "", "Interaction rows are differences of log-odds contrasts (6:4 minus 5:5).", ""] + table(subset,cols)
    header += ["", "## Routing and conclusion boundary", ""]
    header += [f"- {ex}: decomposition route = {gate}." for ex,gate in gates.items()]
    header += ["", "Interpret the main emotion result as an effect averaged across the two offer ratios. "
        "It does not by itself demonstrate a 5:5-specific effect. Nonsignificant interactions do not "
        "establish equality/additivity. Audit-only comparisons must not be promoted to report conclusions.", "",
        "[Supplementary exploratory decomposition](SUPPLEMENT.md). All previously calculated comparisons "
        "remain intact in ../new_contrasts.csv; reporting_routes.csv records their current status.", "",
        "## Diagnostics and next validation", "",
        "Both fits passed the original numerical gate and selected-fit reproduction checks. Some random-effect "
        "covariance directions are weak; see ../covariance_diagnostics.csv. The old actor-intercept, deletion "
        "and simulation checks do not validate these newer fits. Before manuscript use, assess influence and "
        "model adequacy for the new structures under a documented sensitivity plan; do not change the primary "
        "model or correction policy in response to significance. No such refits were run for this reporting edit.", "",
        "## Existing figure", "",
        "The preserved paired figure shows observations and conditional predictions, not permission to infer "
        "audit-only comparisons. No significance annotations were added.", "",
        "![Descriptive paired figure](../paired_figures/equal_offers_paired.png)"]
    (dest/"REPORT.md").write_text("\n".join(header)+"\n",encoding="utf-8")
    (dest/"SUPPLEMENT.md").write_text("\n".join(supp)+"\n",encoding="utf-8")
    with (dest/"reporting_routes.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(routed[0]));writer.writeheader();writer.writerows(routed)
    after = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    assert before == after
    (dest/"QA.json").write_text(json.dumps(dict(inputs_unchanged=True,input_sha256=after,
        routes=gates,main_emotion_rows=len(emotions),routed_rows=len(routed),no_refitting=True),indent=2),encoding="utf-8")
    print(dest/"REPORT.md")


if __name__ == "__main__":
    main()
