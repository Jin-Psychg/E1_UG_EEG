# Custom study materials

This directory contains the author-developed materials needed to inspect the task implementation in Experiments 1 and 2. The materials are provided in their original German where that was the language presented to participants.

## Contents

- `task_code/experiment_1/experiment.sce`: Neurobehavioral Systems Presentation scenario for the sequential task in Experiment 1.
- `task_code/experiment_1/arrays.pcl`: practice trials and the ten counterbalanced 540-trial orders for Experiment 1.
- `task_code/experiment_2/experiment.sce`: Presentation scenario for the simultaneous face-and-offer task in Experiment 2.
- `task_code/experiment_2/arrays.pcl`: practice trials and the ten counterbalanced 840-trial orders for Experiment 2.
- `task_code/experiment_1/scale.png` and `task_code/experiment_2/scale.png`: the gradient used under the post-task visual-analogue scale.
- `participant_instructions.md`: source-language participant instructions, including the experiment-specific initial task screens.
- `post_task_ratings.md`: the three author-developed rating items, response scale, and English glosses.

The source scenarios were encoded in Windows-1252. The shared copies were converted to UTF-8, renamed generically, and had only their identifying header comments replaced. Executable task logic is otherwise unchanged. The archived source SHA-256 values were:

| File | Archived source SHA-256 |
|---|---|
| Experiment 1 scenario | `85BAB77FFCA61A42951478BDA04DC3C9339347D99D205E9DE8422D8DCE898AB5` |
| Experiment 1 trial arrays | `CAEED70AF27A06A29BECB005E269D3CF1E40535CEC07C622896CA2BBE4D449D2` |
| Experiment 2 scenario | `567336D84C004C34A5F44890BF27378F1B4CB477120F15AD2C5BFE1F201719A4` |
| Experiment 2 trial arrays | `88AFBF79EF8A6D3A918DC06E6FAA0085F2DF72DA907E79323FB1A2D50A567A78` |

## Materials not redistributed here

The facial image files are not redistributed. The trial arrays retain the relative image identifiers required to inspect counterbalancing and trial order, but running the task requires access to the corresponding images. Access to the facial stimuli is subject to permission from the source-database owners.

The item text of the Social Value Orientation Slider Measure and the modified Personality Inventory for DSM-5–Brief Form Plus is not duplicated here; readers should consult the original instruments cited in the manuscript. Participant logs, behavioral data, and EEG data are also not included in this materials directory.

Machine-specific Presentation `.exp` configuration files are omitted because they contain workstation paths and device mappings. The study logic and scheduled trial orders are contained in the shared `.sce` and `.pcl` files.
