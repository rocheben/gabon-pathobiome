# Setup used by run_all.R and the analysis scripts: installs missing packages
# from CRAN, then loads pbapply for progress bars (shown also under Rscript).
pkgs <- c("lme4", "blme", "boot", "bnlearn", "Hmsc", "igraph", "ggplot2", "pbapply")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
suppressMessages(library(pbapply))
invisible(pboptions(type = "timer", char = "="))
