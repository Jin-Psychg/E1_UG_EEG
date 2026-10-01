"""Build the Round 52 final mixed-model specification table (R1.2 / SC-07) from the canonical
REPRODUCIBILITY.txt files of the live analysis tree (read only)."""
import glob
import os
import re

ROOT = r"C:/Code/UG_ERP_Project/results"
OUT = r"E:/MyPaper_withARS/draft/phase6_revision/round52_model_specification_table.md"

files = sorted(
    glob.glob(os.path.join(ROOT, "Behavior/*_TwoStage_Bates/*/REPRODUCIBILITY.txt"))
    + glob.glob(os.path.join(ROOT, "EEG/E*_TwoStage_Bates_*/Stage1_TrialLevel/*/REPRODUCIBILITY.txt"))
    + glob.glob(os.path.join(ROOT, "EEG/FacePhase_E1/**/REPRODUCIBILITY.txt"), recursive=True)
)


def field(text, key):
    m = re.search(key + r"\s*:?\s*(.+)", text)
    return m.group(1).strip() if m else "NA"


lines = [
    "# Round 52 - Final mixed-model specification table (R1.2 / SC-07)",
    "",
    "Source: `## Final formula` and convergence block of each `REPRODUCIBILITY.txt` under "
    "`C:/Code/UG_ERP_Project/results/` (2026-08-22 run). Random-slope terms use sum-contrast columns "
    "`RE_<factor>_k` (k = 1..levels-1), one one-dimensional term per contrast (zero-correlation "
    "parameterization), selected by the preregistered buildmer backward-LRT protocol (alpha = .20) with "
    "correlations tested in a final step. Gate = optimizer code 0, non-singular, positive-definite Hessian, "
    "finite SEs/likelihood. `dims` = retained / nominal random-effect dimensions.",
    "",
    "| Analysis | Engine | Gate passed | Optimizer | Singular | PD Hessian | dims | Inference | n | Final formula |",
    "|---|---|---|---|---|---|---|---|---|---|",
]
for f in files:
    t = open(f, encoding="utf-8", errors="replace").read()
    m = re.search(r"## Final formula:\s*\n(.+?)\n\s*\n", t, re.S)
    form = " ".join(m.group(1).split()) if m else "NA"
    eng = "glmmTMB" if "glmmTMB" in field(t, "Engine") else "lme4"
    opt = field(t, "Final optimizer")
    if opt == "NA":
        opt = field(t, "Optimizer")
    label = os.path.relpath(f, ROOT).replace("\\", "/").replace("/REPRODUCIBILITY.txt", "")
    dims = field(t, r"Step 5 \(final dims\)")
    cells = [label, eng, field(t, "Acceptance gate passed"), opt, field(t, "Singular fit"),
             field(t, "Hessian positive-definite"), dims, field(t, "ANOVA test"), field(t, "Observations"), "`" + form + "`"]
    lines.append("| " + " | ".join(cells) + " |")
open(OUT, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"{len(files)} models -> {OUT}")
