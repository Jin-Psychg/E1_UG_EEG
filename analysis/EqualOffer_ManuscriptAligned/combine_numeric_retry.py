"""Merge prior valid fits and recovered fits; verify rather than overwrite existing tables."""
from pathlib import Path
import sys
import pandas as pd

previous, run = map(lambda s: Path(s).resolve(), sys.argv[1:3])
for source, dest, keys in [
    ("marginal_effects.csv", "combined_valid_marginal_effects.csv", ["experiment", "model", "contrast"]),
    ("omnibus.csv", "combined_valid_omnibus.csv", ["experiment", "model", "term"]),
]:
    old = pd.read_csv(previous / source)
    new = pd.read_csv(run / source)
    new = new[new.model != "primary"]
    result = pd.concat([old, new], ignore_index=True)
    assert not result.duplicated(keys).any()
    if source == "marginal_effects.csv":
        result["source_run"] = [str(previous)] * len(old) + [str(run)] * len(new)
    path = run / dest
    if path.exists():
        pd.testing.assert_frame_equal(pd.read_csv(path), result, check_exact=False, rtol=1e-12, atol=1e-12)
    else:
        result.to_csv(path, index=False)
    print(f"Verified merged table: {dest}; {len(result)} rows")
