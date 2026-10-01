# Retrospective exploratory analysis. Fixed before inspecting new model contrasts.
# Run identity is selected from the manuscript's formal-source audit, not timestamps.
experiments <- c("E1", "E2")
emotion_levels <- c("neu", "aff", "dis", "dom", "enj")
rt_bounds_ms <- c(300, 3000)
expected_input_md5 <- c(E1 = "154ae95b7a090c5a79074e197d34f9ac",
                        E2 = "7ebd8bab4d455856a0226d1d0c73307a")
reference_runs <- c(
  E1 = "results/RandomStructure_Behavior_E1_all_fit_20260915_145148/Behavior/E1_TwoStage_Bates",
  E2 = "results/RandomStructure_Behavior_E2_all_fit_20260915_145442/Behavior/E2_TwoStage_Bates")
seed <- 20260928L
random_sd_boundary <- 0.001
max_evaluations <- 10000L
fit_timeout_seconds <- 600
focal_pairs <- list("dom - aff" = c("dom", "aff"),
                    "dom - enj" = c("dom", "enj"))
secondary_pairs <- list("dis - enj" = c("dis", "enj"))
# Each pair contributes within-5:5, within-6:4, and (6:4 minus 5:5).
# Across both experiments: 12 focal and 6 secondary log-odds tests; Holm separately.
# Probability-scale estimates are descriptive companions, with unadjusted 95% CIs.
