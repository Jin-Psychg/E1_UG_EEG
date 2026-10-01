"""Experiment 2 epoch-level feature extraction for Round 52 ERP sensitivities (R1.5, R1.6).

One read-only streaming pass over the canonical single-trial epoch files of the live
analysis tree (C:/Code/UG_ERP_Project/data/02_Pipeline_Output_E2/Method_Regression/
Stimulus_Locked/epochs/VP*_epo.csv; hu-neuro-pipeline output, -0.5 to 1.5 s, 500 Hz,
common-average reference, 0.05-30 Hz, FastICA-corrected, NOT baseline-corrected).

Per epoch it computes
  * ptp_max_uV  : maximum peak-to-peak amplitude over all EEG channels across the full
                  epoch (the quantity hu-neuro-pipeline thresholds: reject={'eeg': x})
  * registered-window ROI means (nominal windows from the Experiment 2 preregistration
    as transcribed in notes/close-reading/Arguments/prereg-reconciliation.md §1/§5):
      N2_reg     200-300 ms   ROI = [VALUE NEEDED: not transcribed] -> provisional FRN ROI
      N400_reg   350-450 ms   Fz, Cz, CPz, Pz
      LPP_reg    500-800 ms   Cz, C1, C2, CP1, CP2, Pz
    plus the -200-0 ms baseline of each ROI (Alday covariate)
  * LPP_pre_recomputed (400-600 ms, Pz Cz C1 C2 CP1 CP2) as an internal validation
    against the canonical trials.csv column LPP_pre.

Writes only to draft/sensitivity/erp-registered-E2_2026-08-23/.
"""
import glob
import hashlib
import os
import sys
import time

import numpy as np
import pandas as pd

SRC = r"C:/Code/UG_ERP_Project/data/02_Pipeline_Output_E2/Method_Regression/Stimulus_Locked"
OUT = r"E:/MyPaper_withARS/draft/sensitivity/erp-registered-E2_2026-08-23"
os.makedirs(OUT, exist_ok=True)

EXPECTED_TRIALS_MD5 = "7ebd8bab4d455856a0226d1d0c73307a"  # REPRODUCIBILITY.txt 2026-08-22

EOG = {"IO1", "VEOG", "HEOG"}
KEY = ["participant_id", "subject", "setting", "block", "index", "stim", "Offers_Other",
       "Offers_You", "Offer_Trigger", "reaction", "RT", "emotion", "event_id", "epoch"]
WINDOWS = {
    "N2_reg": (0.200, 0.300, ["F3", "Fz", "F4", "FC1", "FC2", "Cz"]),   # ROI provisional
    "N400_reg": (0.350, 0.450, ["Fz", "Cz", "CPz", "Pz"]),
    "LPP_reg": (0.500, 0.800, ["Cz", "C1", "C2", "CP1", "CP2", "Pz"]),
    "LPP_pre_recomputed": (0.400, 0.600, ["Pz", "Cz", "C1", "C2", "CP1", "CP2"]),
}
BASELINE = (-0.200, 0.0)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def window_mean(df, tmin, tmax, roi):
    # hu-neuro-pipeline compute_single_trials: mean over ROI channels and samples
    # with tmin <= time <= tmax (inclusive on both ends).
    sel = (df["time"] >= tmin - 1e-9) & (df["time"] <= tmax + 1e-9)
    return df.loc[sel, roi].to_numpy(dtype=np.float64).mean()


def main():
    trials_path = os.path.join(SRC, "trials.csv")
    got = md5(trials_path)
    if got != EXPECTED_TRIALS_MD5:
        sys.exit(f"trials.csv MD5 {got} != expected {EXPECTED_TRIALS_MD5}")
    trials = pd.read_csv(trials_path)
    files = sorted(glob.glob(os.path.join(SRC, "epochs", "VP*_epo.csv")))
    if len(files) != 30:
        sys.exit(f"expected 30 epoch files, found {len(files)}")

    rows = []
    t0 = time.time()
    for fi, path in enumerate(files, 1):
        pid = os.path.basename(path).split("_")[0]
        df = pd.read_csv(path, dtype={"stim": str}, low_memory=False)
        eeg_cols = [c for c in df.columns[df.columns.get_loc("epoch") + 1:]
                    if c not in EOG and not c.startswith("FRN_") and not c.startswith("LPP_")
                    and not c.startswith("P3_") and not c.startswith("Baseline_")]
        for ep, g in df.groupby("epoch", sort=True):
            x = g[eeg_cols].to_numpy(dtype=np.float64)
            ptp = (x.max(axis=0) - x.min(axis=0)).max()
            rec = {k: g[k].iloc[0] for k in KEY}
            rec["n_samples"] = len(g)
            rec["ptp_max_uV"] = ptp
            rec["ptp_max_uV_excl_A2"] = ((x[:, [i for i, c in enumerate(eeg_cols) if c != "A2"]].max(axis=0)
                                         - x[:, [i for i, c in enumerate(eeg_cols) if c != "A2"]].min(axis=0)).max())
            rec["LPP_pre_canonical"] = g["LPP_pre"].iloc[0]
            for name, (a, b, roi) in WINDOWS.items():
                rec[name] = window_mean(g, a, b, roi)
                rec["Baseline_" + name] = window_mean(g, BASELINE[0], BASELINE[1], roi)
            rows.append(rec)
        del df
        print(f"[{fi}/30] {pid} done, {time.time() - t0:.0f} s, eeg_cols={len(eeg_cols)}", flush=True)

    feat = pd.DataFrame(rows)
    feat.to_csv(os.path.join(OUT, "e2_epoch_features.csv"), index=False)

    # Validation 1: rejection at the canonical 200 uV threshold must reproduce trials.csv NaNs.
    merged = trials[["participant_id", "index", "LPP_pre"]].merge(
        feat[["participant_id", "index", "ptp_max_uV", "ptp_max_uV_excl_A2", "LPP_pre_recomputed"]],
        on=["participant_id", "index"], how="left")
    canon_rej = merged["LPP_pre"].isna()
    flag200 = merged["ptp_max_uV"] > 200.0
    flag200_exA2 = merged["ptp_max_uV_excl_A2"] > 200.0
    # Validation 2: recomputed LPP_pre equals canonical value on retained epochs.
    ok = ~canon_rej & merged["LPP_pre_recomputed"].notna()
    diff = (merged.loc[ok, "LPP_pre_recomputed"] - merged.loc[ok, "LPP_pre"]).abs()
    with open(os.path.join(OUT, "validation.txt"), "w", encoding="utf-8") as f:
        f.write(f"trials.csv rows: {len(trials)}; epochs extracted: {len(feat)}; unmatched trial rows: {merged['ptp_max_uV'].isna().sum()}\n")
        f.write(f"canonical rejected (LPP_pre NaN): {int(canon_rej.sum())}\n")
        f.write(f"ptp>200 incl. A2: {int(flag200.sum())}; agreement with canonical NaN: {int((flag200 == canon_rej).sum())}/{len(merged)}\n")
        f.write(f"ptp>200 excl. A2: {int(flag200_exA2.sum())}; agreement with canonical NaN: {int((flag200_exA2 == canon_rej).sum())}/{len(merged)}\n")
        f.write(f"ptp>100 incl. A2: {int((merged['ptp_max_uV'] > 100).sum())}; excl. A2: {int((merged['ptp_max_uV_excl_A2'] > 100).sum())}\n")
        f.write(f"LPP_pre recompute vs canonical on retained epochs: n={int(ok.sum())}, max|diff|={diff.max():.6f}, mean|diff|={diff.mean():.6f}\n")
        f.write(f"N2_reg ROI = provisional FRN ROI (F3 Fz F4 FC1 FC2 Cz) [VALUE NEEDED from E2 preregistration]\n")
    print(open(os.path.join(OUT, "validation.txt")).read())


if __name__ == "__main__":
    main()
