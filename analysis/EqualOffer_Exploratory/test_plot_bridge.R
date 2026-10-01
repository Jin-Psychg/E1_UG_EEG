root <- "C:/Code/UG_ERP_Project"
run <- file.path(root, "results/EqualOffer_replication_20260928_224854")
python <- Sys.getenv("UG_PYTHON", unset = Sys.which("python"))
if (!nzchar(python)) stop("Python is required for manuscript-style figures; set UG_PYTHON. See README.md.")
plot_script <- file.path(root, "analysis", "EqualOffer_Exploratory", "plot_paired.py")
plot_status <- system2(python, c(shQuote(plot_script), "--run", shQuote(run)))
if (plot_status != 0L) stop("Manuscript-style plotting failed; numerical outputs are preserved.")
