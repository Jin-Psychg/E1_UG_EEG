#!/usr/bin/env Rscript
# Run from project root. Each execution creates a new timestamped output.
f <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(if(length(f)) dirname(sub("^--file=", "", f[1])) else getwd(), winslash="/")
while (!file.exists(file.path(root,"UG_ERP_Project.Rproj"))) {
  parent <- dirname(root); if(parent == root) stop("Project root not found"); root <- parent
}
setwd(root)
code <- file.path(root,"analysis","EqualOffer_ManuscriptAligned")
out <- file.path(root,"results",format(Sys.time(),"EqualOffer_Aligned_%Y%m%d_%H%M%S"))
rscript <- file.path(R.home("bin"),"Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"),"Rscript")
python <- Sys.getenv("UG_PYTHON", unset=Sys.which("python"))
if(!nzchar(python)) stop("Set UG_PYTHON to Python with numpy, pandas and matplotlib")
run <- function(command, args) {
  status <- system2(command, vapply(args, shQuote, character(1)))
  if(status != 0) stop("Stage failed; preserve its log/output. Exit code: ", status)
}
run(rscript,c("--vanilla",file.path(code,"run_analysis.R"),paste0("--out=",out)))
run(rscript,c("--vanilla",file.path(code,"verify_results.R"),out))
run(python,c(file.path(code,"plot_paired.py"),"--run",out))
run(python,c(file.path(code,"render_standard_report.py"),out))
cat("Completed report:",file.path(out,"standard_reporting","REPORT.md"),"\n")
