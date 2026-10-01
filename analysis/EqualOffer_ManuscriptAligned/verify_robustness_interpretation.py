"""Audit the bounded interpretation against all completed result rows."""
from pathlib import Path
import sys
import pandas as pd

run = Path(sys.argv[1]).resolve()
dest = run / "validation_claim_audit.csv"
if dest.exists():
    raise FileExistsError("Preserve the existing audit")
influence = pd.read_csv(run / "influence_summary.csv")
omnibus = pd.read_csv(run / "omnibus.csv")
status = pd.read_csv(run / "sensitivity_status.csv")
simulation = pd.read_csv(run / "simulation_summary.csv")
refits = pd.read_csv(run / "independent_refit_comparison.csv")
attempts = pd.read_csv(run / "fit_attempts.csv")
audit = []


def record(ex, claim, value, source, rows):
    audit.append(dict(experiment=ex, claim=claim, value=value, source=source, row_filter=rows))


for ex in ["E1", "E2"]:
    row = status[status.experiment == ex].iloc[0]
    record(ex, "Valid / total deletion fits", f"{row.valid_deletions}/{row.total_deletions}",
           "sensitivity_status.csv", f"experiment={ex}")
    e = influence[(influence.experiment == ex) & (influence.contrast != "5:5 - 6:4")]
    significant = e[e.primary_p_BH < .05]
    flips = int(significant.direction_changes.sum())
    assert flips == 0
    record(ex, "Direction reversals among primary BH-significant emotion contrasts, valid fits only", flips,
           "influence_summary.csv", f"experiment={ex}; emotion pairs; primary_p_BH<.05")
    switches = int((((e.primary_p_BH < .05) & (e.p_BH_max >= .05)) |
                    ((e.primary_p_BH >= .05) & (e.p_BH_min < .05))).sum())
    assert switches > 0 if ex == "E1" else switches == 0
    record(ex, "Emotion pairs with any BH-threshold classification change, valid fits only", switches,
           "influence_summary.csv", f"experiment={ex}; emotion pairs; p_BH_min/max vs primary_p_BH")
    p = omnibus[(omnibus.experiment == ex) & omnibus.model.str.startswith("delete_") &
                (omnibus.term == "emotion")]["Pr(>Chisq)"]
    assert len(p) == row.valid_deletions and p.max() < .05
    record(ex, "Emotion omnibus p minimum, valid deletions", p.min(), "omnibus.csv", f"{ex}; delete_*; emotion")
    record(ex, "Emotion omnibus p maximum, valid deletions", p.max(), "omnibus.csv", f"{ex}; delete_*; emotion")
    z = simulation[simulation.experiment == ex]
    assert len(z) == 13 and not z.outside_pointwise_envelope.any()
    record(ex, "Summaries outside pointwise simulation envelope / checked", f"0/{len(z)}",
           "simulation_summary.csv", f"experiment={ex}; all rows")
    for r in attempts[(attempts.experiment == ex) & (attempts.model == "actor_intercept")].itertuples():
        assert not r.valid
        record(ex, f"Actor fit {r.optimizer} diagnostics", r.diagnostics, "fit_attempts.csv",
               f"experiment={ex}; actor_intercept; optimizer={r.optimizer}")
    r = refits[refits.experiment == ex]
    assert len(r) == 3 and r.valid.all() and (r.maximum_absolute_difference == 0).all()
    record(ex, "Independent targeted refit maximum numerical difference", r.maximum_absolute_difference.max(),
           "independent_refit_comparison.csv", f"experiment={ex}; all metrics")
pd.DataFrame(audit).to_csv(dest, index=False)
print(f"Interpretation verified against {len(audit)} source-linked claims; human verification remains false.")
