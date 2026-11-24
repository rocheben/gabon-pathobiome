#######################################################################
# Reproduces all analyses, tables (resultats/) and figures (figures/) of
# Abbate et al., "Exposure risk, not within-host interactions, structures
# a human pathogen community".
# Usage: source("code/run_all.R") or Rscript code/run_all.R. The working
# directory is set automatically to the project root (the folder that
# contains code/, data/, resultats/ and figures/); all paths are relative to it.
# Approximate total run time: about 1 h on a laptop (11_hmsc_jsdm.R is the longest).
#######################################################################
# ---- Locate the project root, whatever the current working directory
is_root <- function(p) file.exists(file.path(p, "code", "run_all.R")) && dir.exists(file.path(p, "data"))
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)       # source()
if (is.null(this_file)) {                                                                # Rscript
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) this_file <- normalizePath(f)
}
candidates <- c(if (!is.null(this_file)) dirname(dirname(this_file)), getwd(), dirname(getwd()))
root <- Filter(is_root, candidates)
if (length(root) == 0) stop("Project root not found: set the working directory to the folder ",
                            "that contains code/ and data/, then source(\"code/run_all.R\").")
setwd(root[1]); message("Working directory: ", getwd())
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)

dir.create("resultats", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
analyses <- c("01_data_preparation.R", "02_pairwise.R", "03_association_screening.R",
              "04_targeted_models_rarefaction.R", "05_bayesian_networks.R", "06_village_scale.R",
              "07_ebov_optical_density.R", "08_robustness.R", "09_null_models_type1.R",
              "10_null_models_power.R", "11_hmsc_jsdm.R", "12_scales.R")
figures <- c("fig1_community_structure.R", "fig2_scales.R", "fig3_pairwise.R", "fig4_association_screening.R",
             "figS1_rarefaction.R", "figS2_null_benchmark.R", "figS3_robustness.R", "figS4_bayesian_networks.R",
             "figS5_hmsc_residual.R", "figS6_village_EBOV_Plasmodium.R", "figS7_EBOV_optical_density.R")
all_scripts <- c(analyses, figures); t_start <- Sys.time()
for (i in seq_along(all_scripts)) {
  s <- all_scripts[i]; t0 <- Sys.time()
  message(sprintf("\n[%d/%d] %s", i, length(all_scripts), s))
  source(file.path("code", s), local = new.env())
  message(sprintf("  done in %s (total %s)", format(round(Sys.time() - t0, 1)), format(round(Sys.time() - t_start, 1))))
}
