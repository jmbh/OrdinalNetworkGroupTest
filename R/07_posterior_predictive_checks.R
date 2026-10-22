# Posterior predictive checks and indicator diagnostics ----------------------
# Inputs: archived bgmCompare fit and group-specific analysis data.
# Outputs: numerical PPC summaries, indicator diagnostics, and Appendix Figure A6.

source(file.path("R", "00_config.R"))
check_packages(c("bgms", "ggplot2", "patchwork"))

fit <- readRDS(paths$fit)
observed_groups <- list(
  Female = as.matrix(readRDS(paths$female_data)),
  Male = as.matrix(readRDS(paths$male_data))
)
labels <- fit$arguments$data_columnnames
num_categories <- fit$arguments$num_categories
p <- length(labels)
edge_names <- fit$raw_samples$parameter_names$pairwise_baseline

# The original response 5 ("Often") is represented by category code 4 after
# bgms recodes the seven response categories to 0,...,6.
high_response <- 4L
nrep <- 500L
simulation_iterations <- 1000L

pairwise_statistics <- function(data) {
  correlations <- stats::cor(data, method = "spearman")[lower.tri(diag(p))]
  upper_tail <- vapply(edge_names, function(edge) {
    variables <- strsplit(edge, "-", fixed = TRUE)[[1]]
    mean(data[, variables[1]] >= high_response & data[, variables[2]] >= high_response)
  }, numeric(1))
  data.frame(
    edge = rep(edge_names, 2),
    statistic = rep(c("Spearman correlation", "Upper-tail co-occurrence"), each = length(edge_names)),
    value = c(correlations, upper_tail),
    stringsAsFactors = FALSE
  )
}

burden_statistics <- function(data) {
  burden <- rowSums(data >= high_response)
  data.frame(
    burden = 0:p,
    value = tabulate(burden + 1L, nbins = p + 1L) / nrow(data)
  )
}

# Combine chains, then draw the same posterior iteration for both groups in a
# replicate. Columns contain baseline parameters followed by group differences.
main_draws <- do.call(rbind, fit$raw_samples$main)
pairwise_draws <- do.call(rbind, fit$raw_samples$pairwise)
n_main <- ncol(main_draws) / 2L
n_pairwise <- ncol(pairwise_draws) / 2L
stopifnot(n_main == sum(num_categories), n_pairwise == choose(p, 2))

set.seed(20261020)
draw_indices <- sample.int(nrow(main_draws), nrep, replace = FALSE)

simulate_group <- function(draw_index, group_index, replicate_index) {
  direction <- if (group_index == 1L) -0.5 else 0.5
  main_vector <- main_draws[draw_index, seq_len(n_main)] +
    direction * main_draws[draw_index, n_main + seq_len(n_main)]
  pairwise_vector <- pairwise_draws[draw_index, seq_len(n_pairwise)] +
    direction * pairwise_draws[draw_index, n_pairwise + seq_len(n_pairwise)]

  main_matrix <- matrix(main_vector, nrow = p, byrow = TRUE)
  pairwise_matrix <- matrix(0, p, p, dimnames = list(labels, labels))
  pairwise_matrix[lower.tri(pairwise_matrix)] <- pairwise_vector
  pairwise_matrix <- pairwise_matrix + t(pairwise_matrix)

  simulated <- bgms::simulate_mrf(
    num_states = nrow(observed_groups[[group_index]]),
    num_variables = p,
    num_categories = num_categories,
    pairwise = pairwise_matrix,
    main = main_matrix,
    variable_type = "ordinal",
    iter = simulation_iterations,
    seed = 20261020L + 2L * replicate_index + group_index
  )
  colnames(simulated) <- labels
  simulated
}

replicated_pairwise <- vector("list", nrep * 2L)
replicated_burden <- vector("list", nrep * 2L)
cursor <- 1L
for (replicate_index in seq_len(nrep)) {
  for (group_index in 1:2) {
    simulated <- simulate_group(draw_indices[replicate_index], group_index, replicate_index)
    group_name <- names(observed_groups)[group_index]
    pairwise_result <- pairwise_statistics(simulated)
    pairwise_result$group <- group_name
    pairwise_result$replicate <- replicate_index
    burden_result <- burden_statistics(simulated)
    burden_result$group <- group_name
    burden_result$replicate <- replicate_index
    replicated_pairwise[[cursor]] <- pairwise_result
    replicated_burden[[cursor]] <- burden_result
    cursor <- cursor + 1L
  }
  if (replicate_index %% 50L == 0L) message("Completed ", replicate_index, " of ", nrep, " PPC replicates.")
}
replicated_pairwise <- do.call(rbind, replicated_pairwise)
replicated_burden <- do.call(rbind, replicated_burden)

observed_pairwise <- do.call(rbind, lapply(names(observed_groups), function(group_name) {
  result <- pairwise_statistics(observed_groups[[group_name]])
  result$group <- group_name
  result
}))
observed_burden <- do.call(rbind, lapply(names(observed_groups), function(group_name) {
  result <- burden_statistics(observed_groups[[group_name]])
  result$group <- group_name
  result
}))

summarize_replications <- function(replicated, keys) {
  groups <- split(replicated, interaction(replicated[keys], drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(groups, function(x) {
    out <- x[1, keys, drop = FALSE]
    out$predictive_median <- stats::median(x$value)
    out$predictive_lower <- stats::quantile(x$value, 0.025, names = FALSE)
    out$predictive_upper <- stats::quantile(x$value, 0.975, names = FALSE)
    out
  }))
}

pairwise_summary <- summarize_replications(replicated_pairwise, c("group", "statistic", "edge"))
pairwise_summary <- merge(pairwise_summary, observed_pairwise, by = c("group", "statistic", "edge"))
names(pairwise_summary)[names(pairwise_summary) == "value"] <- "observed"
pairwise_summary$outside_95_interval <- with(
  pairwise_summary, observed < predictive_lower | observed > predictive_upper
)

burden_summary <- summarize_replications(replicated_burden, c("group", "burden"))
burden_summary <- merge(burden_summary, observed_burden, by = c("group", "burden"))
names(burden_summary)[names(burden_summary) == "value"] <- "observed"
burden_summary$outside_95_interval <- with(
  burden_summary, observed < predictive_lower | observed > predictive_upper
)

utils::write.csv(
  pairwise_summary, file.path(paths$diagnostics, "ppc_pairwise_summary.csv"), row.names = FALSE
)
utils::write.csv(
  burden_summary, file.path(paths$diagnostics, "ppc_symptom_burden_summary.csv"), row.names = FALSE
)

# Indicator diagnostics, including chain-specific inclusion rates. -----------
indicator_names <- fit$raw_samples$parameter_names$indicators
pairwise_indicator_index <- grep("\\(pairwise\\)$", indicator_names)
indicator_summary <- fit$posterior_summary_indicator[pairwise_indicator_index, ]
chain_rates <- vapply(
  fit$raw_samples$indicator,
  function(chain) colMeans(chain[, pairwise_indicator_index, drop = FALSE]),
  numeric(length(pairwise_indicator_index))
)
indicator_diagnostics <- data.frame(
  edge = sub(" \\(pairwise\\)$", "", indicator_summary$parameter),
  pip = indicator_summary$mean,
  inclusion_bf = indicator_summary$mean / (1 - indicator_summary$mean),
  chain_1 = chain_rates[, 1],
  chain_2 = chain_rates[, 2],
  chain_3 = chain_rates[, 3],
  chain_4 = chain_rates[, 4],
  mcse = indicator_summary$mcse,
  effective_sample_size = indicator_summary$n_eff,
  rhat = indicator_summary$Rhat
)
indicator_diagnostics <- indicator_diagnostics[order(-indicator_diagnostics$inclusion_bf), ]
utils::write.csv(
  indicator_diagnostics,
  file.path(paths$diagnostics, "indicator_diagnostics.csv"),
  row.names = FALSE
)

# Figure ---------------------------------------------------------------------
plot_pairwise <- function(statistic_name, x_label, panel_label) {
  data <- pairwise_summary[pairwise_summary$statistic == statistic_name, ]
  data$edge <- factor(data$edge, levels = rev(edge_names))
  ggplot2::ggplot(data, ggplot2::aes(y = edge)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = predictive_lower, xmax = predictive_upper),
      orientation = "y", width = 0, color = "#5B8DB8", linewidth = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(x = predictive_median), color = "#5B8DB8", shape = 1, size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(x = observed), color = "black", size = 1.8) +
    ggplot2::facet_wrap(~group, nrow = 1) +
    ggplot2::labs(x = x_label, y = NULL, title = panel_label) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 10)
    )
}

correlation_plot <- plot_pairwise(
  "Spearman correlation", "Spearman correlation", "A. Pairwise ordinal associations"
)
tail_plot <- plot_pairwise(
  "Upper-tail co-occurrence", "Proportion both rated Often or higher",
  "B. Upper-tail co-occurrence"
)

burden_plot_data <- burden_summary
burden_plot_data$group <- factor(burden_plot_data$group, levels = c("Female", "Male"))
burden_plot <- ggplot2::ggplot(
  burden_plot_data, ggplot2::aes(x = burden, group = group)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = predictive_lower, ymax = predictive_upper, fill = group),
    alpha = 0.22, color = NA
  ) +
  ggplot2::geom_line(ggplot2::aes(y = predictive_median, color = group), linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(y = observed, color = group), size = 2) +
  ggplot2::facet_wrap(~group, nrow = 1) +
  ggplot2::scale_x_continuous(breaks = 0:p) +
  ggplot2::scale_color_manual(values = c(Female = "#386CB0", Male = "#D95F02"), guide = "none") +
  ggplot2::scale_fill_manual(values = c(Female = "#386CB0", Male = "#D95F02"), guide = "none") +
  ggplot2::labs(
    x = "Number of symptoms rated Often or higher", y = "Proportion",
    title = "C. Symptom-burden distribution"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 10)
  )

combined_plot <- correlation_plot / tail_plot / burden_plot +
  patchwork::plot_layout(heights = c(1.45, 1.45, 0.8))
ggplot2::ggsave(
  file.path(paths$figures, "FigureA6_PPU_Posterior_Predictive_Checks.pdf"),
  combined_plot, width = 9, height = 11, units = "in"
)

cat("Created posterior predictive checks and indicator diagnostics.\n")
