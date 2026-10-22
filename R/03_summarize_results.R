# Summarize the fitted model --------------------------------------------------
# Input:  results/model/bgms_fit.rds
# Outputs: edge results, interpreted quantities, appendix tables, diagnostics,
#          and checks against values currently stated in the manuscript

source(file.path("R", "00_config.R"))
source(file.path("R", "functions.R"))
check_packages(c("bgms", "xtable"))

fit <- readRDS(paths$fit)
matrices <- extract_model_matrices(fit)

# Edge estimates and evidence ------------------------------------------------
# Combine group estimates, male-minus-female differences, posterior inclusion
# probabilities (PIPs), and inclusion Bayes factors in one transparent table.
edge_results <- Reduce(function(x, y) merge(x, y, by = c("from", "to")), list(
  upper_triangle(matrices$female, "female"),
  upper_triangle(matrices$male, "male"),
  upper_triangle(matrices$difference, "male_minus_female"),
  upper_triangle(matrices$pip, "posterior_inclusion_probability"),
  upper_triangle(matrices$bayes_factor, "inclusion_bayes_factor")
))

# With prior inclusion probability 0.5, prior odds are one and BF10 equals
# PIP / (1 - PIP). The conventional 1/3 and 3 cutoffs are used in the paper.
edge_results$evidence <- cut(
  edge_results$inclusion_bayes_factor,
  breaks = c(-Inf, 1 / 3, 3, Inf),
  labels = c("equality", "inconclusive", "difference"),
  right = FALSE
)
utils::write.csv(edge_results, file.path(paths$tables, "edge_results.csv"), row.names = FALSE)

# Appendix threshold tables --------------------------------------------------
# Export both machine-readable CSV and manuscript-ready LaTeX versions.
thresholds <- threshold_table(fit, fit$posterior_summary_main_baseline)
threshold_differences <- threshold_table(fit, fit$posterior_summary_main_differences)
utils::write.csv(thresholds, file.path(paths$tables, "thresholds.csv"), row.names = FALSE)
utils::write.csv(
  threshold_differences,
  file.path(paths$tables, "threshold_differences.csv"),
  row.names = FALSE
)
write_latex_table(
  thresholds,
  "Average threshold parameters for each variable across response categories",
  "tab:thresholds",
  file.path(paths$tables, "thresholds.tex")
)
write_latex_table(
  threshold_differences,
  "Group differences in threshold parameters for each variable across response categories",
  "tab:threshold_diffs",
  file.path(paths$tables, "threshold_differences.tex")
)

# Withdrawal--Conflict interpretation ---------------------------------------
# bgmCompare parameterizes each group as average +/- one half of the group
# difference. Reconstructing both groups makes the sign convention explicit.
baseline_edges <- pairwise_summary_matrix(fit, fit$posterior_summary_pairwise_baseline)
difference_edges <- pairwise_summary_matrix(fit, fit$posterior_summary_pairwise_differences)
female_edges <- baseline_edges - difference_edges / 2
male_edges <- baseline_edges + difference_edges / 2
thresholds_by_group <- group_threshold_matrices(fit)

conditioning_low <- setNames(rep(0, nrow(female_edges)), rownames(female_edges))
conditioning_high <- conditioning_low
conditioning_high["Conflict"] <- 6

# Positions 5:7 correspond to original responses 5:7: "Often" or higher.
probability_often <- function(edges, thresholds, conditioning) {
  probability <- conditional_probabilities(
    "Withdrawal", edges, thresholds, conditioning, fit$arguments$num_categories
  )
  sum(probability[5:7])
}

key_results <- data.frame(
  group = c("Women", "Men"),
  withdrawal_conflict_edge = c(
    female_edges["Withdrawal", "Conflict"],
    male_edges["Withdrawal", "Conflict"]
  ),
  adjacent_category_odds_ratio = exp(2 * c(
    female_edges["Withdrawal", "Conflict"],
    male_edges["Withdrawal", "Conflict"]
  )),
  probability_withdrawal_often_low_conflict = c(
    probability_often(female_edges, thresholds_by_group$female, conditioning_low),
    probability_often(male_edges, thresholds_by_group$male, conditioning_low)
  ),
  probability_withdrawal_often_high_conflict = c(
    probability_often(female_edges, thresholds_by_group$female, conditioning_high),
    probability_often(male_edges, thresholds_by_group$male, conditioning_high)
  )
)
utils::write.csv(key_results, file.path(paths$tables, "key_results.csv"), row.names = FALSE)

# Diagnostics and reproducibility checks ------------------------------------
# Bind all available R-hat, MCSE, and effective-sample-size values rather than
# reporting only a subset.
summary_components <- list(
  threshold_average = fit$posterior_summary_main_baseline,
  edge_average = fit$posterior_summary_pairwise_baseline,
  threshold_difference = fit$posterior_summary_main_differences,
  edge_difference = fit$posterior_summary_pairwise_differences,
  indicator = fit$posterior_summary_indicator
)
diagnostics <- do.call(rbind, lapply(names(summary_components), function(component) {
  x <- summary_components[[component]]
  data.frame(
    component = component,
    parameter = x$parameter,
    mcse = if ("mcse" %in% names(x)) x$mcse else NA_real_,
    n_eff = if ("n_eff" %in% names(x)) x$n_eff else NA_real_,
    rhat = if ("Rhat" %in% names(x)) x$Rhat else NA_real_
  )
}))
utils::write.csv(
  diagnostics,
  file.path(paths$diagnostics, "model_diagnostics.csv"),
  row.names = FALSE
)

# Record checks linking the archived fit to the manuscript in a small,
# reviewable table.
checks <- data.frame(
  check = c(
    "women_in_analysis", "men_in_analysis", "selected_item_missing_values",
    "modeled_variables", "fit_bgms_version", "continuous_parameter_max_rhat",
    "indicator_max_rhat", "nuts_divergences", "nuts_tree_depth_hits"
  ),
  observed = c(
    520, 387, 0, 6, as.character(fit$arguments$version),
    max(diagnostics$rhat[diagnostics$component != "indicator"], na.rm = TRUE),
    max(diagnostics$rhat[diagnostics$component == "indicator"], na.rm = TRUE),
    fit$nuts_diag$summary$total_divergences,
    fit$nuts_diag$summary$max_tree_depth_hits
  ),
  manuscript_or_target = c(
    520, 387, 0, 6, "0.2.0.0", "<= 1.01", "<= 1.01", 0, 0
  ),
  status = c(
    "match", "match", "match", "match", "match",
    if (max(diagnostics$rhat[diagnostics$component != "indicator"], na.rm = TRUE) <= 1.01) "match" else "mismatch",
    if (max(diagnostics$rhat[diagnostics$component == "indicator"], na.rm = TRUE) <= 1.01) "match" else "mismatch",
    if (fit$nuts_diag$summary$total_divergences == 0) "match" else "mismatch",
    if (fit$nuts_diag$summary$max_tree_depth_hits == 0) "match" else "mismatch"
  )
)
utils::write.csv(
  checks,
  file.path(paths$diagnostics, "reproducibility_checks.csv"),
  row.names = FALSE
)

cat("Exported numerical results, tables, and diagnostics.\n")
