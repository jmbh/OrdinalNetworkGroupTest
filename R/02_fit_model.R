# Fit or load the ordinal-network comparison --------------------------------
# Inputs:  data/derived/ppu_female.rds and ppu_male.rds
# Output:  results/model/bgms_fit.rds
# By default the archived fit is retained; set REFIT_MODEL=true via --refit to
# rerun the sampler from the prepared data.

source(file.path("R", "00_config.R"))
check_packages("bgms")

force_refit <- identical(Sys.getenv("REFIT_MODEL"), "true")
if (!file.exists(paths$fit) || force_refit) {
  female <- readRDS(paths$female_data)
  male <- readRDS(paths$male_data)

  # Group 1 is female and Group 2 is male. Edge differences are selected with
  # inclusion indicators; threshold differences are estimated without selection.
  fit <- bgms::bgmCompare(
    x = female,
    y = male,
    difference_selection = TRUE,
    main_difference_selection = FALSE,
    # Make the bgms 0.2.0.0 defaults explicit so the archived analysis remains
    # reproducible if package defaults change in a later release.
    difference_scale = 1,
    difference_family = "Normal",
    difference_prior = bgms::bernoulli_prior(0.5),
    interaction_prior = bgms::normal_prior(scale = 1),
    threshold_prior = bgms::beta_prime_prior(alpha = 0.5, beta = 0.5),
    na_action = "impute",
    iter = 20000,
    warmup = 10000,
    chains = 4,
    cores = 4,
    seed = 1,
    display_progress = "total"
  )
  saveRDS(fit, paths$fit, version = 3)
  cat("Model fitted and saved to", paths$fit, "\n")
} else {
  # Loading the stored fit reproduces the exact Monte Carlo output in the paper.
  fit <- readRDS(paths$fit)
  cat("Using archived fit. Pass --refit to rerun the sampler.\n")
}

# Guard against accidentally using a fit from another variable set or sample.
stopifnot(
  inherits(fit, "bgmCompare"),
  identical(fit$arguments$data_columnnames, unname(item_labels)),
  fit$arguments$num_cases == 907L,
  identical(as.character(fit$arguments$version), "0.2.0.0")
)
