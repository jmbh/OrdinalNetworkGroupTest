# Plot prior distributions discussed in the manuscript ----------------------
# Outputs: four appendix PDFs covering the current defaults and alternatives.

source(file.path("R", "00_config.R"))
check_packages("ggplot2")

prior_theme <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.25, color = "grey85"),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "grey85"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.spacing.x = grid::unit(0.8, "lines"),
    plot.margin = ggplot2::margin(6, 6, 6, 6)
  )

save_prior_plot <- function(plot, filename) {
  ggplot2::ggsave(
    filename = file.path(paths$figures, filename), plot = plot,
    device = grDevices::cairo_pdf, width = 9, height = 3, units = "in"
  )
}

shade_middle_half <- function(data, quantile_function) {
  intervals <- unique(data[c("panel", "value1", "value2")])
  intervals$lower <- mapply(quantile_function, 0.25, intervals$value1, intervals$value2)
  intervals$upper <- mapply(quantile_function, 0.75, intervals$value1, intervals$value2)
  merge(data, intervals[c("panel", "lower", "upper")], by = "panel")
}

make_settings <- function(value1, value2, labels) {
  data.frame(value1 = value1, value2 = value2, panel = factor(labels, levels = labels))
}

# Normal priors: current default for average interactions and differences ----
normal_labels <- c(
  "More concentrated\nscale = 0.5", "Default\nscale = 1",
  "More dispersed\nscale = 2"
)
normal_settings <- make_settings(c(0.5, 1, 2), 0, normal_labels)
normal_data <- merge(data.frame(parameter = seq(-6, 6, length.out = 3000)), normal_settings)
normal_data$density <- stats::dnorm(normal_data$parameter, sd = normal_data$value1)
normal_data <- shade_middle_half(
  normal_data, function(probability, scale, unused) stats::qnorm(probability, sd = scale)
)
normal_plot <- ggplot2::ggplot(normal_data, ggplot2::aes(parameter, density)) +
  ggplot2::geom_area(fill = "#DDE8F2", color = NA) +
  ggplot2::geom_area(
    data = subset(normal_data, parameter >= lower & parameter <= upper),
    fill = "#83A9CC", color = NA
  ) +
  ggplot2::geom_line(color = "#28597F", linewidth = 0.7) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::coord_cartesian(xlim = c(-6, 6), ylim = c(0, 0.82), expand = FALSE) +
  ggplot2::labs(x = "Parameter value", y = "Prior density") + prior_theme
save_prior_plot(normal_plot, "FigureA1_Normal_Priors.pdf")

# Cauchy priors: optional heavier-tailed alternative -------------------------
cauchy_labels <- c(
  "More concentrated\nscale = 0.5", "Constructor default\nscale = 1",
  "More dispersed\nscale = 2"
)
cauchy_settings <- make_settings(c(0.5, 1, 2), 0, cauchy_labels)
cauchy_data <- merge(data.frame(parameter = seq(-6, 6, length.out = 3000)), cauchy_settings)
cauchy_data$density <- stats::dcauchy(
  cauchy_data$parameter, location = 0, scale = cauchy_data$value1
)
cauchy_data <- shade_middle_half(
  cauchy_data,
  function(probability, scale, unused) stats::qcauchy(probability, scale = scale)
)
cauchy_plot <- ggplot2::ggplot(cauchy_data, ggplot2::aes(parameter, density)) +
  ggplot2::geom_area(fill = "#D8E6E7", color = NA) +
  ggplot2::geom_area(
    data = subset(cauchy_data, parameter >= lower & parameter <= upper),
    fill = "#82AEB0", color = NA
  ) +
  ggplot2::geom_line(color = "#1D4D4F", linewidth = 0.7) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::coord_cartesian(xlim = c(-6, 6), ylim = c(0, 0.66), expand = FALSE) +
  ggplot2::labs(x = "Parameter value", y = "Prior density") + prior_theme
save_prior_plot(cauchy_plot, "FigureA2_Cauchy_Priors.pdf")

# Beta-prime priors for average thresholds ----------------------------------
threshold_labels <- c(
  "More concentrated\nalpha = beta = 2", "Default\nalpha = beta = 0.5",
  "More dispersed\nalpha = beta = 0.25"
)
threshold_settings <- make_settings(c(2, 0.5, 0.25), c(2, 0.5, 0.25), threshold_labels)
threshold_data <- merge(
  data.frame(parameter = seq(-14, 14, length.out = 4000)), threshold_settings
)
threshold_data$probability <- stats::plogis(threshold_data$parameter)
threshold_data$density <- stats::dbeta(
  threshold_data$probability, threshold_data$value1, threshold_data$value2
) * threshold_data$probability * (1 - threshold_data$probability)
threshold_data <- shade_middle_half(
  threshold_data,
  function(probability, alpha, beta) stats::qlogis(stats::qbeta(probability, alpha, beta))
)
threshold_plot <- ggplot2::ggplot(threshold_data, ggplot2::aes(parameter, density)) +
  ggplot2::geom_area(fill = "#E9DFCF", color = NA) +
  ggplot2::geom_area(
    data = subset(threshold_data, parameter >= lower & parameter <= upper),
    fill = "#C8A978", color = NA
  ) +
  ggplot2::geom_line(color = "#6D4C1F", linewidth = 0.7) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::coord_cartesian(xlim = c(-14, 14), ylim = c(0, 0.4), expand = FALSE) +
  ggplot2::labs(x = expression("Threshold parameter " * mu[ih]), y = "Prior density") +
  prior_theme
save_prior_plot(threshold_plot, "FigureA3_BetaPrime_Threshold_Priors.pdf")

# Beta priors for the optional Beta--Bernoulli inclusion model ---------------
beta_labels <- c(
  "More concentrated\nalpha = beta = 4", "Default\nalpha = beta = 1",
  "More dispersed\nalpha = beta = 0.5"
)
beta_settings <- make_settings(c(4, 1, 0.5), c(4, 1, 0.5), beta_labels)
beta_data <- merge(data.frame(parameter = seq(0.005, 0.995, length.out = 3000)), beta_settings)
beta_data$density <- stats::dbeta(beta_data$parameter, beta_data$value1, beta_data$value2)
beta_data <- shade_middle_half(
  beta_data, function(probability, alpha, beta) stats::qbeta(probability, alpha, beta)
)
beta_plot <- ggplot2::ggplot(beta_data, ggplot2::aes(parameter, density)) +
  ggplot2::geom_area(fill = "#E7DAE9", color = NA) +
  ggplot2::geom_area(
    data = subset(beta_data, parameter >= lower & parameter <= upper),
    fill = "#B58BBC", color = NA
  ) +
  ggplot2::geom_line(color = "#613868", linewidth = 0.7) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 5), expand = FALSE) +
  ggplot2::labs(
    x = expression("Shared inclusion probability " * pi), y = "Prior density"
  ) + prior_theme
save_prior_plot(beta_plot, "FigureA4_BetaBernoulli_Priors.pdf")

cat("Created prior figures.\n")
