# Re-run only the existing E1 trait screen using the verified revised base fits.
local({
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the revised primary run directory.")
if (!file.exists("UG_ERP_Project.Rproj")) stop("Run from the project root.")
.libPaths(c(normalizePath("renv/library/R-4.3/x86_64-w64-mingw32"), .libPaths()))
reference_root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
run_id <- paste0("TraitRevision_E1_", format(Sys.time(), "%Y%m%d_%H%M%S"))
run_root <- file.path("results", run_id)
if (dir.exists(run_root)) stop("Output already exists.")
dir.create(file.path(run_root, "Reports"), recursive = TRUE)
options(ug.run_id = run_id, ug.stages = character(), ug.trait_stages = "E1",
        ug.trait_model_root = reference_root, ug.run_trait_screen = TRUE)
cat("RUN_ROOT:", normalizePath(run_root, winslash = "/"), "\n")
rmarkdown::render("Sta_Behaviour_E1_E2_Integrative.Rmd",
 output_dir = normalizePath(file.path(run_root, "Reports")),
 intermediates_dir = normalizePath(file.path(run_root, "Reports")),
 knit_root_dir = getwd(), envir = new.env(), quiet = FALSE)
cat("FINISHED:", normalizePath(run_root, winslash = "/"), "\n")
})
