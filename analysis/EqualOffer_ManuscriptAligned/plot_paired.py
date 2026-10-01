"""Render existing exploratory results using the manuscript's paired-facet design.

No fitting, p-value adjustment, filtering or participant exclusion occurs here.
python analysis/EqualOffer_ManuscriptAligned/plot_paired.py --run results/RUN_NAME
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, help="Completed, verified analysis directory")
    parser.add_argument("--out", type=Path, help="New figure directory; defaults to RUN/paired_figures")
    args = parser.parse_args()
    run = args.run.resolve(strict=True)
    out = args.out.resolve() if args.out else run / "paired_figures"
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    if (run / "STATUS.txt").read_text().splitlines()[0] != "COMPLETE":
        raise ValueError("Analysis is not complete")
    if not (run / "VERIFICATION.txt").is_file():
        raise ValueError("Independent analysis verification is required")
    style_path = Path(__file__).with_name("manuscript_style.json")
    style = json.loads(style_path.read_text(encoding="utf-8"))
    participant_path, model_path = run / "participant_cells.csv", run / "new_probabilities.csv"
    before = {str(p): sha256(p) for p in (participant_path, model_path, run / "new_contrasts.csv")}
    p = pd.read_csv(participant_path)
    m = pd.read_csv(model_path)
    required_p = {"experiment", "participant_id", "emotion", "allocation", "valid_trials", "rejection_rate"}
    required_m = {"experiment", "emotion", "allocation", "prob", "asymp.LCL", "asymp.UCL"}
    if not required_p.issubset(p) or not required_m.issubset(m):
        raise ValueError("Missing required figure columns")
    if p.duplicated(["experiment", "participant_id", "emotion", "allocation"]).any():
        raise ValueError("Participant cells must be unique")
    if m.duplicated(["experiment", "emotion", "allocation"]).any():
        raise ValueError("Model cells must be unique")
    if set(p.emotion) != set(style["order"]) or set(p.allocation) != {"5:5", "6:4"}:
        raise ValueError("Unexpected expression/allocation labels")
    if not m[["prob", "asymp.LCL", "asymp.UCL"]].map(lambda x: 0 <= x <= 1).all().all():
        raise ValueError("Model estimates outside probability bounds")
    if not ((m["asymp.LCL"] <= m.prob) & (m.prob <= m["asymp.UCL"])).all():
        raise ValueError("Model confidence intervals do not bracket estimates")
    out.mkdir(parents=True)
    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7})
    fig, axes = plt.subplots(2, 5, figsize=(180 / 25.4, 140 / 25.4), sharey=True,
                             gridspec_kw={"hspace": .48, "wspace": .10})
    coordinates, cell_audit, sample_sizes = [], [], {}
    for row, ex in enumerate(("E1", "E2")):
        ep = p[p.experiment == ex]
        em = m[m.experiment == ex]
        ids = sorted(ep.participant_id.unique())
        sample_sizes[ex] = len(ids)
        rng = np.random.default_rng(style["jitter_seed"])
        offsets = dict(zip(ids, rng.uniform(-style["jitter_half_width"], style["jitter_half_width"], len(ids))))
        ref = em[(em.emotion == "neu") & (em.allocation == "6:4")]
        if len(ref) != 1:
            raise ValueError("Neutral-6:4 model reference is missing")
        ref_pct = 100 * float(ref.prob.iloc[0])
        for col, emotion in enumerate(style["order"]):
            ax = axes[row, col]
            rgb = mcolors.to_rgb(style["colors"][emotion])
            subset = ep[ep.emotion == emotion]
            pairs = subset.pivot(index="participant_id", columns="allocation", values="rejection_rate")
            lines_drawn, points_drawn = 0, 0
            for participant, values in pairs.iterrows():
                x = np.array([0, 1]) + offsets[participant]
                y = 100 * values.reindex(["5:5", "6:4"]).to_numpy(dtype=float)
                if np.isfinite(y).all():
                    ax.plot(x, y, color=rgb, alpha=style["line_alpha"], lw=.6, zorder=2,
                            solid_capstyle="round")
                    lines_drawn += 1
                for j, allocation in enumerate(("5:5", "6:4")):
                    if np.isfinite(y[j]):
                        ax.scatter(x[j], y[j], s=8, facecolor=(*rgb, style["dot_alpha"]),
                                   edgecolor="none", zorder=3)
                        points_drawn += 1
                        coordinates.append({"experiment": ex, "emotion": emotion, "participant_id": participant,
                            "allocation": allocation, "x": x[j], "observed_percent": y[j]})
            for j, allocation in enumerate(("5:5", "6:4")):
                cell = em[(em.emotion == emotion) & (em.allocation == allocation)]
                if len(cell) != 1:
                    raise ValueError(f"Missing model cell: {ex}, {emotion}, {allocation}")
                r = cell.iloc[0]
                value, low, high = 100 * r.prob, 100 * r["asymp.LCL"], 100 * r["asymp.UCL"]
                ax.errorbar(j, value, yerr=[[value - low], [high - value]], fmt="D", markersize=4.5,
                    markerfacecolor=rgb, markeredgecolor="#111111", markeredgewidth=.6,
                    ecolor="#111111", elinewidth=1, capsize=2.5, capthick=1, zorder=6)
            ax.axhline(ref_pct, color="#555555", lw=.6, ls=(0, (4, 3)), alpha=.55, zorder=1)
            luminance = .299 * rgb[0] + .587 * rgb[1] + .114 * rgb[2]
            ax.set_title(style["labels"][emotion], fontsize=8.5, fontweight="bold", pad=5,
                color="#1a1a1a" if luminance > .6 else "white",
                bbox=dict(facecolor=rgb, edgecolor="none", boxstyle="square,pad=0.28"))
            ax.set_xlim(-.5, 1.5); ax.set_ylim(-4, 104)
            ax.set_xticks([0, 1], ["5:5", "6:4"]); ax.set_yticks(np.arange(0, 101, 20))
            ax.tick_params(labelsize=8, length=3, width=.7)
            if col == 0:
                ax.set_ylabel("Rejection rate (%)", fontsize=8.5)
                ax.text(-.40, 1.22, f"{'ab'[row]}  Experiment {row + 1} (N = {len(ids)})",
                        transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left")
            cell_audit.append({"experiment": ex, "emotion": emotion, "participant_lines": lines_drawn,
                "participant_points": points_drawn, "model_diamonds": 2,
                "missing_rate_cells": int(subset.rejection_rate.isna().sum()), "color": style["colors"][emotion]})
    fig.subplots_adjust(left=.09, right=.985, top=.89, bottom=.19)
    fig.supxlabel("Proposer : responder allocation", fontsize=8.5, y=.105)
    grey = (.35, .35, .35)
    handles = [
        mlines.Line2D([], [], color=(*grey, .6), lw=.8, marker="o", markersize=3,
            markeredgecolor="none", label="Participant (5:5–6:4 pair)"),
        mlines.Line2D([], [], marker="D", linestyle="None", markersize=4.5,
            markerfacecolor=grey, markeredgecolor="#111111", markeredgewidth=.6, label="Model estimate ± 95% CI"),
        mlines.Line2D([], [], color="#555555", lw=.8, ls=(0, (4, 3)), label="Neutral–6:4 model reference")]
    fig.legend(handles=handles, loc="lower center", ncol=1, frameon=False, fontsize=7.5,
               bbox_to_anchor=(.54, -.015), handlelength=2.1, labelspacing=.35)
    fig.canvas.draw()
    clipped = []
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    for text in fig.findobj(match=matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if box.x0 < -2 or box.y0 < -2 or box.x1 > width + 2 or box.y1 > height + 2:
                clipped.append(text.get_text())
    # bbox_inches='tight' includes external text; record preliminary bounds for visual QA.
    stem = out / "equal_offers_paired"
    for ext, dpi in (("png", 300), ("svg", None), ("pdf", None), ("tiff", 600)):
        fig.savefig(stem.with_suffix("." + ext), dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    pd.DataFrame(coordinates).to_csv(out / "plotted_participant_coordinates.csv", index=False)
    pd.DataFrame(cell_audit).to_csv(out / "figure_cell_audit.csv", index=False)
    shutil.copy2(model_path, out / "source_model_probabilities.csv")
    shutil.copy2(style_path, out / style_path.name)
    shutil.copy2(Path(__file__), out / Path(__file__).name)
    after = {p: sha256(Path(p)) for p in before}
    if before != after:
        raise RuntimeError("A source result changed during rendering")
    caption = (
        "Exploratory rejection of equal and mildly unequal offers. a, Experiment 1; b, Experiment 2. "
        "Colored facet headers follow the manuscript expression order and palette. Thin lines connect each "
        "participant's observed rejection percentages at 5:5 and 6:4; dots mark the individual condition rates. "
        f"E1 N = {sample_sizes['E1']} and E2 N = {sample_sizes['E2']}; zero-rejection participants are retained. "
        "The same horizontal jitter is used for a participant at both allocations and across expression panels "
        "within an experiment. Diamonds and black whiskers show the saved fair-only GLMM predictions with "
        "random effects set to zero and unadjusted 95% confidence intervals. They are not arithmetic participant "
        "means or population-integrated probabilities. The dashed line marks that experiment's Neutral-6:4 "
        "prediction as a visual reference. No significance stars are shown. Estimates are from the manuscript-aligned "
        "five-step models; the companion contrast tables use BH-FDR within each experiment and comparison family."
    )
    (out / "caption.txt").write_text(caption + "\n", encoding="utf-8")
    (out / "FIGURE_NOTES.md").write_text(
        "# Manuscript-style figure revision\n\n## Figure contract\n\n"
        "- Claim: equal-offer and mildly unequal-offer rejection varies across expressions and participants; "
        "the figure does not identify a psychological mechanism.\n"
        "- Evidence: two experiment rows, five expression columns; paired observed rates plus existing model intervals.\n"
        "- Archetype: quantitative grid. Backend: Python/matplotlib, matching the manuscript workflow.\n"
        "- Export: 180 x 140 mm nominal canvas; tight bounding box includes the shared element key; "
        "editable SVG/PDF, 300-dpi PNG and 600-dpi TIFF.\n"
        "- Change: update model estimates to the new manuscript-aligned fits. "
        "Keep participant pairing, expression headers, diamond estimates and neutral reference from the manuscript.\n"
        "- Adaptation: rows denote experiments rather than choice/RT because this analysis only refitted choice; "
        "x-axis is 5:5/6:4; the reference is Neutral-6:4. Font sizes increased for legibility.\n"
        "- Risk: heavy zero overlap remains visible as overlapping points, not suppressed; "
        "wide between-person variation and conditional model estimates must not be conflated.\n"
        "- Statistics: this renderer preserves the new manuscript-aligned model outputs; it performs no fitting.\n\n"
        "## Figure\n\n![Paired figure](equal_offers_paired.png)\n\n" + caption + "\n",
        encoding="utf-8")
    receipt = {"input_hashes_before": before, "input_hashes_after": after, "unchanged": before == after,
        "style": style, "sample_sizes": sample_sizes, "cells": cell_audit,
        "pre_tight_export_external_text": clipped, "numerical_analysis_rerun": False,
        "script_sha256": sha256(Path(__file__))}
    (out / "figure_provenance.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
