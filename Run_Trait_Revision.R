# Run the covariate screen with full-model random-structure selection.
local({
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1] %in% c("E1", "E2")) stop("Supply E1 or E2.")
stage <- args[1]
if (!file.exists("UG_ERP_Project.Rproj")) stop("Run from the project root.")
.libPaths(c(normalizePath("renv/library/R-4.3/x86_64-w64-mingw32"), .libPaths()))
run_id <- paste0("TraitSelection_", stage, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
run_root <- file.path("results", run_id)
if (dir.exists(run_root)) stop("Output already exists.")
dir.create(file.path(run_root, "Reports"), recursive = TRUE)
options(ug.run_id = run_id, ug.stages = character(), ug.trait_stages = stage,
        ug.run_trait_screen = TRUE)
cat("RUN_ROOT:", normalizePath(run_root, winslash = "/"), "\n")
rmarkdown::render("Sta_Behaviour_E1_E2_Integrative.Rmd",
 output_dir = normalizePath(file.path(run_root, "Reports")),
 intermediates_dir = normalizePath(file.path(run_root, "Reports")),
 knit_root_dir = getwd(), envir = new.env(), quiet = FALSE)
cat("FINISHED:", normalizePath(run_root, winslash = "/"), "\n")
})
