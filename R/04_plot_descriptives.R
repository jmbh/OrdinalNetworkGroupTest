# Plot observed marginal response distributions -----------------------------
# Inputs:  the two prepared group data files
# Output:  figures/Figure1_PPU_Marginal_Distributions.pdf

source(file.path("R", "00_config.R"))
check_packages("ggplot2")

female <- readRDS(paths$female_data)
male <- readRDS(paths$male_data)
combined <- rbind(female, male)

# Convert the six wide item columns to plotting data with within-item
# proportions. Factor levels retain the questionnaire's 1--7 category order.
plot_data <- do.call(rbind, lapply(names(combined), function(variable) {
  frequency <- prop.table(table(combined[[variable]]))
  data.frame(
    variable = factor(variable, levels = names(combined)),
    category = factor(names(frequency), levels = as.character(1:7)),
    relative_frequency = as.numeric(frequency)
  )
}))

# Facets share one scale so category frequencies remain comparable by eye.
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = category, y = relative_frequency)
) +
  ggplot2::geom_col(fill = "grey50", colour = "grey20", width = 0.8) +
  ggplot2::facet_wrap(~variable, nrow = 2, ncol = 3) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::labs(x = "Category", y = "Relative frequency") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 6, b = 1)),
    panel.spacing.y = grid::unit(1.5, "lines"),
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(paths$figures, "Figure1_PPU_Marginal_Distributions.pdf"),
  plot, width = 10, height = 6.5, units = "in"
)
cat("Created marginal-distribution figure.\n")
