# Run from the UG_ERP_Project root with Rscript --vanilla.
# Usage: Rscript --vanilla Run_Statistical_Revision.R Behavior E1 GLMM_rejection preflight
# Replace preflight with fit; use all to run all models/components for that experiment.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Supply method (Behavior/Alday/Traditional), experiment, target (or all), preflight/fit.")
method <- match.arg(args[1], c("Behavior", "Alday", "Traditional"))
experiment <- match.arg(args[2], c("E1", "E2", "Integrative"))
mode <- match.arg(args[4], c("preflight", "fit"))
if (method != "Behavior" && experiment == "Integrative") stop("EEG has no Integrative target.")
if (!file.exists("UG_ERP_Project.Rproj")) stop("Run from the project root.")
project_library <- file.path("renv", "library", paste0("R-", R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]), R.version$platform)
if (!dir.exists(project_library)) stop("Project renv library is missing; restore the recorded environment first.")
.libPaths(c(normalizePath(project_library), .libPaths()))
# Fail before pacman can install anything in the Rmd setup.
required <- c("pacman", "rmarkdown", "here", "readr", "readxl", "dplyr", "tidyr", "stringr", "tibble", "purrr", "lme4", "lmerTest", "glmmTMB", "pbkrtest", "car", "parameters", "performance", "buildmer", "R.utils", "emmeans", "effectsize", "sjPlot", "sjmisc", "ggplot2", "ggpubr", "gridExtra", "knitr", "kableExtra", "openxlsx")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing recorded dependencies: ", paste(missing, collapse = ", "))
expected_versions <- c(lme4 = "1.1-38", lmerTest = "3.1-3", glmmTMB = "1.1.14",
                       buildmer = "2.12", emmeans = "2.0.0", Matrix = "1.6-5")
if (getRversion() != "4.3.3") stop("This reproduction requires the recorded R 4.3.3 environment.")
for (p in names(expected_versions)) {
  if (packageVersion(p) != package_version(expected_versions[[p]])) stop("Recorded package version differs: ", p)
}
script <- if (method == "Behavior") "Sta_Behaviour_E1_E2_Integrative.Rmd" else paste0("Sta_EEG_OfferPhase_Ref", method, ".Rmd")
run_id <- paste("RandomStructure", method, experiment, args[3], mode, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
run_root <- file.path("results", run_id)
if (dir.exists(run_root)) stop("Output already exists: ", run_root)
dir.create(file.path(run_root, "Reports"), recursive = TRUE)
options(ug.run_id = run_id, ug.stages = experiment, ug.experiment = experiment,
        ug.verify_reference = TRUE, ug.preflight_only = mode == "preflight",
        ug.run_trait_screen = FALSE, ug.run_rating_control = FALSE,
        ug.run_baseline_diagnostics = FALSE)
if (args[3] != "all") {
  if (method == "Behavior") options(ug.analyses = args[3]) else options(ug.components = args[3])
}
writeLines(c(paste("Command:", paste(args, collapse = " ")),
             paste("Source MD5:", unname(tools::md5sum(script))),
             paste("Config MD5:", unname(tools::md5sum("config/random_structure_revision.R")))),
           file.path(run_root, "RUN_RECEIPT.txt"))
cat("RUN_ROOT:", normalizePath(run_root, winslash = "/"), "\n")
rmarkdown::render(script, output_dir = normalizePath(file.path(run_root, "Reports")),
                  intermediates_dir = normalizePath(file.path(run_root, "Reports")),
                  knit_root_dir = getwd(), envir = new.env(), quiet = FALSE)
ref_env <- new.env()
sys.source("config/random_structure_revision.R", envir = ref_env)
ids <- names(ref_env$revision_models)
ids <- ids[startsWith(ids, paste(method, experiment, "", sep = "__"))]
if (args[3] != "all") ids <- ids[endsWith(ids, paste0("__", args[3]))]
if (!length(ids)) stop("No configured target matched.")
expected <- vapply(ids, function(id) {
  path <- file.path(run_root, ref_env$revision_models[[id]]$path)
  if (mode == "preflight") file.path(dirname(path), "INPUT_FIXED_VERIFICATION.txt") else path
}, character(1))
if (!all(file.exists(expected))) stop("Missing accepted outputs: ", paste(ids[!file.exists(expected)], collapse = ", "))
cat("FINISHED:", normalizePath(run_root, winslash = "/"), "\n")
