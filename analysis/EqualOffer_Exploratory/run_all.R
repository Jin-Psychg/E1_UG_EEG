#!/usr/bin/env Rscript
# One-command entry point. Run from the project root or from this script directory.
# Rscript --vanilla analysis/EqualOffer_Exploratory/run_all.R
flags <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(flags)) dirname(normalizePath(sub("^--file=", "", flags[1]), winslash = "/")) else getwd()
while (!file.exists(file.path(root, "UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if (parent == root) stop("Open UG_ERP_Project.Rproj before sourcing this entry point")
  root <- parent
}
code <- file.path(root, "analysis", "EqualOffer_Exploratory")
output <- file.path(root, "results", paste0("EqualOffer_replication_", format(Sys.time(), "%Y%m%d_%H%M%S")))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
steps <- list(
  c(shQuote(file.path(code, "run_analysis.R")), "--stage=all", shQuote(paste0("--out=", output))),
  c(shQuote(file.path(code, "verify_results.R")), shQuote(output)),
  c(shQuote(file.path(code, "render_report.R")), shQuote(output)))
for (step in steps) {
  status <- system2(rscript, c("--vanilla", step))
  if (status != 0L) stop("A stage failed (exit ", status, "). Preserve outputs and inspect the error above: ", output)
}
message("Analysis, independent verification and report complete: ", file.path(output, "REPORT.md"))
