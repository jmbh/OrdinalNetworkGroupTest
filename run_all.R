# Reproduce the complete analysis in the order shown below.
# Run from the repository root: Rscript run_all.R [--refit]
# Without --refit, the archived model is used to reproduce the exact outputs.

# Parse the only supported command-line option and pass it to the fitting script.
arguments <- commandArgs(trailingOnly = TRUE)
unknown <- setdiff(arguments, "--refit")
if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
Sys.setenv(REFIT_MODEL = if ("--refit" %in% arguments) "true" else "false")

source(file.path("R", "00_config.R"))
check_packages(c("bgms", "ggraph", "ggplot2", "patchwork", "readxl", "tidygraph", "xtable"))

scripts <- c(
  "01_prepare_data.R",
  "02_fit_model.R",
  "03_summarize_results.R",
  "04_plot_descriptives.R",
  "05_plot_model.R",
  "06_plot_priors.R",
  "07_posterior_predictive_checks.R",
  "08_prior_sensitivity.R"
)

# Each script runs in its own environment, preventing temporary objects from
# one analysis stage from leaking into the next.
for (script in scripts) {
  message("\n==> R/", script)
  source(file.path("R", script), local = new.env(parent = globalenv()))
}

# Record the software environment after all packages have been loaded.
capture.output(
  sessionInfo(),
  file = file.path("results", "diagnostics", "session_info.txt")
)
message("\nReproduction complete. See figures/ and results/.")
