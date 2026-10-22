# Check sensitivity to the Normal difference-slab scale ----------------------
# Inputs:  archived bgmCompare fit (or the refit produced by --refit)
# Outputs: a reusable sensitivity object, complete curve and anchor summaries,
#          a table of classification changes, and the appendix figure

source(file.path("R", "00_config.R"))
check_packages(c("bgms", "ggplot2"))

fit <- readRDS(paths$fit)
refit_requested <- identical(Sys.getenv("REFIT_MODEL", unset = "false"), "true")

# The endpoints match the concentrated and dispersed Normal scales shown in
# the prior appendix. Intermediate geometric anchors keep the importance-
# reweighted curve connected while changing only the difference-slab scale.
anchors <- c(0.5, sqrt(0.5), 1, sqrt(2), 2)

if (refit_requested || !file.exists(paths$prior_sensitivity)) {
  sensitivity <- bgms::prior_sensitivity_check(
    fit,
    anchors = anchors,
    evidence_threshold = 3,
    seed = 2401,
    cores = 4,
    keep_fits = FALSE,
    verbose = TRUE
  )
  saveRDS(sensitivity, paths$prior_sensitivity)
} else {
  sensitivity <- readRDS(paths$prior_sensitivity)
  message("Using archived prior-sensitivity analysis. Pass --refit to rerun it.")
}

stopifnot(
  isTRUE(all.equal(sensitivity$anchors, anchors)),
  identical(sensitivity$evidence_threshold, 3),
  isTRUE(all.equal(sensitivity$chosen_scale, 1))
)

# Export the complete pairwise curve in a transparent long format.
pairwise_index <- grep("\\(pairwise\\)$", sensitivity$edge_names)
curve <- do.call(rbind, lapply(pairwise_index, function(edge_index) {
  data.frame(
    edge = sub(" \\(pairwise\\)$", "", sensitivity$edge_names[edge_index]),
    difference_scale = sensitivity$scales,
    inclusion_bayes_factor = exp(sensitivity$log_bf[, edge_index]),
    log_bf_mcse = sensitivity$log_bf_mcse[, edge_index],
    verdict = sensitivity$verdict[, edge_index],
    importance_ess = sensitivity$curve$ess,
    anchor_used = sensitivity$curve$anchor_used,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  curve,
  file.path(paths$diagnostics, "prior_sensitivity_curve.csv"),
  row.names = FALSE
)

# Summarize the three scales displayed in the prior appendix.
reported_scales <- c(0.5, 1, 2)
reported_rows <- match(reported_scales, sensitivity$scales)
anchor_summary <- do.call(rbind, lapply(pairwise_index, function(edge_index) {
  data.frame(
    edge = sub(" \\(pairwise\\)$", "", sensitivity$edge_names[edge_index]),
    difference_scale = reported_scales,
    inclusion_bayes_factor = exp(sensitivity$log_bf[reported_rows, edge_index]),
    verdict = sensitivity$verdict[reported_rows, edge_index],
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  anchor_summary,
  file.path(paths$tables, "prior_sensitivity_anchors.csv"),
  row.names = FALSE
)

mover_names <- sensitivity$edges$edge[
  grepl("\\(pairwise\\)$", sensitivity$edges$edge) &
    sensitivity$edges$mover == "moved-beyond-wobble"
]
mover_names <- sub(" \\(pairwise\\)$", "", mover_names)
movers <- anchor_summary[anchor_summary$edge %in% mover_names, ]
utils::write.csv(
  movers,
  file.path(paths$tables, "prior_sensitivity_movers.csv"),
  row.names = FALSE
)

counts <- as.data.frame.matrix(table(
  factor(anchor_summary$difference_scale, levels = reported_scales),
  factor(anchor_summary$verdict, levels = c("presence", "undecided", "absence"))
))
counts$difference_scale <- reported_scales
counts <- counts[, c("difference_scale", "presence", "undecided", "absence")]
utils::write.csv(
  counts,
  file.path(paths$tables, "prior_sensitivity_counts.csv"),
  row.names = FALSE
)

# Fail visibly if regenerated outputs no longer match the values transcribed in
# the manuscript.
target_anchors <- anchor_summary[anchor_summary$edge == "Withdrawal-Conflict", ]
stopifnot(
  identical(round(target_anchors$inclusion_bayes_factor, 2), c(12.65, 2.64, 1.62)),
  identical(as.integer(counts$presence), c(3L, 0L, 0L)),
  identical(as.integer(counts$undecided), c(0L, 2L, 1L)),
  identical(as.integer(counts$absence), c(12L, 13L, 14L)),
  length(unique(movers$edge)) == 3L
)

# Plot the edge interpreted in the main text. The narrow ribbon displays
# Monte Carlo uncertainty in log(BF10), not posterior uncertainty in the edge.
target <- curve[curve$edge == "Withdrawal-Conflict", ]
target$lower <- exp(log(target$inclusion_bayes_factor) - 1.96 * target$log_bf_mcse)
target$upper <- exp(log(target$inclusion_bayes_factor) + 1.96 * target$log_bf_mcse)
target$is_anchor <- target$difference_scale %in% sensitivity$anchors

plot_sensitivity <- ggplot2::ggplot(
  target,
  ggplot2::aes(x = difference_scale, y = inclusion_bayes_factor)
) +
  ggplot2::annotate(
    "rect", xmin = 0.5, xmax = 2, ymin = 3, ymax = Inf,
    fill = "#3B82C4", alpha = 0.08
  ) +
  ggplot2::annotate(
    "rect", xmin = 0.5, xmax = 2, ymin = 1e-3, ymax = 1 / 3,
    fill = "#E3B341", alpha = 0.12
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower, ymax = upper),
    fill = "#315A7D", alpha = 0.18
  ) +
  ggplot2::geom_line(color = "#315A7D", linewidth = 0.9) +
  ggplot2::geom_point(
    data = target[target$is_anchor, ],
    color = "#315A7D", fill = "white", shape = 21, size = 2.5, stroke = 0.8
  ) +
  ggplot2::geom_hline(
    yintercept = c(1 / 3, 3), linetype = "dashed", color = "grey40"
  ) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dotted", color = "grey40") +
  ggplot2::scale_x_continuous(
    breaks = c(0.5, 1, 1.5, 2), limits = c(0.5, 2)
  ) +
  ggplot2::scale_y_log10(
    breaks = c(1 / 3, 1, 3, 10),
    labels = c("1/3", "1", "3", "10")
  ) +
  ggplot2::coord_cartesian(ylim = c(0.2, 20)) +
  ggplot2::labs(
    x = "Normal difference-slab scale",
    y = expression("Inclusion Bayes factor " * BF[10])
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(color = "black"),
    axis.text = ggplot2::element_text(color = "black")
  )

ggplot2::ggsave(
  file.path(paths$figures, "FigureA7_PPU_Prior_Sensitivity.pdf"),
  plot_sensitivity,
  width = 7,
  height = 4.5,
  device = grDevices::cairo_pdf
)

cat("Exported prior-sensitivity summaries and appendix figure.\n")
