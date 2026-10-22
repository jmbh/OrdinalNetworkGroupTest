# Plot model estimates and evidence -----------------------------------------
# Input:   results/model/bgms_fit.rds
# Outputs: the four-panel network figure and numerical heatmap appendix figure

source(file.path("R", "00_config.R"))
source(file.path("R", "functions.R"))
check_packages(c("bgms", "ggraph", "ggplot2", "patchwork", "tidygraph"))

fit <- readRDS(paths$fit)
matrices <- extract_model_matrices(fit)
labels <- fit$arguments$data_columnnames

edge_palette <- c(negative = "#C64747", positive = "#3B6FB6")
evidence_palette <- c(equality = "#eeb004", inconclusive = "#bfbfbf", difference = "#36648b")

# Rotate the circular layout so the same node appears at the top of each panel.
circle_layout <- function(network) {
  layout <- ggraph::create_layout(network, layout = "circle")
  top <- which.max(layout$y)
  angle <- pi / 2 - atan2(layout$y[top], layout$x[top])
  old_x <- layout$x
  old_y <- layout$y
  layout$x <- old_x * cos(angle) - old_y * sin(angle)
  layout$y <- old_x * sin(angle) + old_y * cos(angle)
  layout
}

# Draw one parameter network. Edge width encodes absolute magnitude, edge
# colour encodes sign, and labels are limited to the largest estimates.
network_plot <- function(matrix, title, label_n = 10) {
  edges <- upper_triangle(matrix)
  edges$sign <- ifelse(edges$value < 0, "negative", "positive")
  rank <- rank(-abs(edges$value), ties.method = "first")
  edges$label <- ifelse(rank <= label_n, format_number(edges$value, 3), "")
  graph <- tidygraph::tbl_graph(
    nodes = data.frame(name = labels), edges = edges, directed = FALSE
  )
  layout <- circle_layout(graph)

  # Place selected value labels at edge midpoints without labeling every edge.
  label_edges <- edges[nzchar(edges$label), ]
  from_index <- match(label_edges$from, layout$name)
  to_index <- match(label_edges$to, layout$name)
  label_edges$x <- (layout$x[from_index] + layout$x[to_index]) / 2
  label_edges$y <- (layout$y[from_index] + layout$y[to_index]) / 2

  ggraph::ggraph(layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = abs(value), colour = sign),
      alpha = 1, lineend = "round"
    ) +
    ggplot2::geom_label(
      data = label_edges,
      ggplot2::aes(x = x, y = y, label = label, colour = sign),
      inherit.aes = FALSE, fill = "white", linewidth = 0,
      label.padding = grid::unit(0.02, "lines"),
      label.r = grid::unit(0, "lines"), size = 2.6, show.legend = FALSE
    ) +
    ggraph::geom_node_point(size = 18.5, colour = "grey25") +
    ggraph::geom_node_text(ggplot2::aes(label = name), colour = "white", size = 2.3) +
    ggraph::scale_edge_width(range = c(0.5, 4), limits = c(0, 1), guide = "none") +
    ggraph::scale_edge_colour_manual(values = edge_palette, guide = "none") +
    ggplot2::scale_colour_manual(values = edge_palette, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.14)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.14)) +
    ggplot2::ggtitle(title) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "plain"),
      plot.margin = ggplot2::margin(5, 10, 5, 10)
    )
}

# Evidence panel -------------------------------------------------------------
# Unlike the parameter panels, these edge colours encode evidence categories.
evidence_edges <- upper_triangle(matrices$bayes_factor, "bayes_factor")
evidence_edges$evidence <- cut(
  evidence_edges$bayes_factor,
  breaks = c(-Inf, 1 / 3, 3, Inf),
  labels = c("equality", "inconclusive", "difference"),
  right = FALSE
)
evidence_graph <- tidygraph::tbl_graph(
  nodes = data.frame(name = labels), edges = evidence_edges, directed = FALSE
)
evidence_plot <- ggraph::ggraph(circle_layout(evidence_graph)) +
  ggraph::geom_edge_link(
    ggplot2::aes(colour = evidence), width = 1.25, alpha = 1, lineend = "round"
  ) +
  ggraph::geom_node_point(size = 18.5, colour = "grey25") +
  ggraph::geom_node_text(ggplot2::aes(label = name), colour = "white", size = 2.3) +
  ggraph::scale_edge_colour_manual(values = evidence_palette, guide = "none") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.14)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.14)) +
  ggplot2::ggtitle("Evidence Network") +
  ggraph::theme_graph(base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "plain"),
    plot.margin = ggplot2::margin(5, 10, 5, 10)
  )

# Assemble the manuscript's 2 x 2 network figure.
network_figure <- patchwork::wrap_plots(
  network_plot(matrices$female, "Female"),
  network_plot(matrices$male, "Male"),
  network_plot(matrices$difference, "Difference: Male - Female"),
  evidence_plot,
  ncol = 2
)
ggplot2::ggsave(
  file.path(paths$figures, "Figure2_PPU_Network_Group_Comparison.pdf"),
  network_figure, width = 8, height = 8, units = "in"
)

# Heatmap appendix figure ----------------------------------------------------
# Only the nonredundant upper triangle is shown because all matrices are symmetric.
heatmap_data <- function(matrix) {
  out <- upper_triangle(matrix)
  out$from <- factor(out$from, levels = rev(labels))
  out$to <- factor(out$to, levels = labels)
  out$label <- format_number(out$value, 2)
  out
}

# Female and male panels use a common colour range; the smaller differences
# receive their own symmetric range so their signs remain visible.
parameter_heatmap <- function(matrix, title, limit) {
  data <- heatmap_data(matrix)
  ggplot2::ggplot(data, ggplot2::aes(x = to, y = from, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2) +
    ggplot2::scale_fill_gradient2(
      low = edge_palette[["negative"]], mid = "white", high = edge_palette[["positive"]],
      midpoint = 0, limits = c(-limit, limit), guide = "none"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 13)
    )
}

# The evidence heatmap reports the exact BF10 in addition to its category.
evidence_heatmap_data <- evidence_edges
evidence_heatmap_data$from <- factor(evidence_heatmap_data$from, levels = rev(labels))
evidence_heatmap_data$to <- factor(evidence_heatmap_data$to, levels = labels)
evidence_heatmap_data$label <- format_number(evidence_heatmap_data$bayes_factor, 2)
evidence_heatmap <- ggplot2::ggplot(
  evidence_heatmap_data,
  ggplot2::aes(x = to, y = from, fill = evidence)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2) +
  ggplot2::scale_fill_manual(values = evidence_palette, drop = FALSE, name = expression(BF[10])) +
  ggplot2::coord_fixed() +
  ggplot2::labs(title = "Evidence: inclusion Bayes factors", x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_family = "sans") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(hjust = 0.5, size = 13)
  )

group_limit <- max(abs(c(matrices$female, matrices$male)))
difference_limit <- max(abs(matrices$difference))
heatmap_figure <- patchwork::wrap_plots(
  parameter_heatmap(matrices$female, "Female", group_limit),
  parameter_heatmap(matrices$male, "Male", group_limit),
  parameter_heatmap(matrices$difference, "Difference: Male - Female", difference_limit),
  evidence_heatmap,
  ncol = 2
)
ggplot2::ggsave(
  file.path(paths$figures, "FigureA3_PPU_Network_Heatmaps.pdf"),
  heatmap_figure, width = 10, height = 9, units = "in"
)
cat("Created network and heatmap figures.\n")
