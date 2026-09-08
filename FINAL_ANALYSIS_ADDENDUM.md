# Final-analysis code addendum

The manuscript's final-analysis code is available in
[Final_Analysis_Code_Addendum_Review.zip](Final_Analysis_Code_Addendum_Review.zip).
Use the file page's **Download raw file** button to download the archive.
It contains 69 files (315,699 bytes compressed).

SHA-256: `05c92ce159c6eb85f7fa1e48e5473b6fc7b644194b8da75420f235b27535e6c0`

## Contents and use

Extract the ZIP into a separate directory. Start with its `README.md`,
`SOURCE_MAP.md` and `RUNNING.md`. The package contains the selected E1/E2
HDDM notebooks and configurations, completed sensitivity and recovery
follow-up code, provenance checks and selected aggregate reports.
Its `manifest.json` records the delivered file hashes.

The older notebooks in `hDDMdetails/` remain historical reference copies;
use the addendum for the manuscript's final HDDM analysis implementations.
The addendum supplements this repository's behavioral and EEG pipelines.
Keep the extracted directory intact so its internal relative links resolve.

## Validation and remaining limitations

Archive integrity and manifest hashes were checked before upload. Python
scripts passed syntax checks; notebook cell sources were preserved while
outputs were removed. No scientific analysis was rerun for this upload.
Participant-level data, EEG, fitted models, posterior draws, facial stimuli
and local agent instructions are excluded. Reproduction requires the
controlled-access inputs identified in `RUNNING.md`.

The original Docker image's registry digest and full package inventory
remain to be verified. This is a code-distribution snapshot, not a claim
that a fresh-environment reproduction has passed. The archive retains its
original `Review` filename and pre-upload README as an immutable packaging
record; this repository now distributes that exact archive. Future
environment-documentation updates should use a versioned successor and
record its new checksum.

The root `.gitignore` permits only this specific ZIP; other archive and
data exclusions remain in effect.
