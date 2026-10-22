# Helpers shared by the results and figure scripts --------------------------
# These functions contain the bookkeeping needed to turn bgmCompare's compact
# posterior summaries into symmetric matrices, tables, and probabilities.

# Reconstruct female, male, and male-minus-female edge matrices. Also convert
# posterior inclusion probabilities to BF10; prior inclusion odds equal one.
extract_model_matrices <- function(fit) {
  labels <- fit$arguments$data_columnnames
  p <- length(labels)
  estimates <- stats::coef(fit)

  female <- male <- matrix(0, p, p, dimnames = list(labels, labels))
  female[lower.tri(female)] <- estimates$pairwise_effects_groups[, 1]
  male[lower.tri(male)] <- estimates$pairwise_effects_groups[, 2]
  female <- female + t(female)
  male <- male + t(male)
  difference <- male - female

  # Verify the reconstruction against bgmCompare's stored difference summary.
  stopifnot(isTRUE(all.equal(
    difference[lower.tri(difference)],
    fit$posterior_summary_pairwise_differences$mean,
    tolerance = 1e-10,
    check.attributes = FALSE
  )))

  pip <- estimates$indicators
  bayes_factor <- pip / (1 - pip)
  diag(bayes_factor) <- NA_real_

  list(
    female = female,
    male = male,
    difference = difference,
    pip = pip,
    bayes_factor = bayes_factor
  )
}

# Convert a symmetric matrix to one row per unique, undirected edge.
upper_triangle <- function(x, value_name = "value") {
  index <- which(upper.tri(x), arr.ind = TRUE)
  out <- data.frame(
    from = rownames(x)[index[, 1]],
    to = colnames(x)[index[, 2]],
    value = x[index],
    stringsAsFactors = FALSE
  )
  names(out)[3] <- value_name
  out
}

# Format rounded labels while avoiding the visually confusing value "-0.00".
format_number <- function(x, digits = 2) {
  x[abs(x) < 0.5 * 10^-digits] <- 0
  sprintf(paste0("%.", digits, "f"), x)
}

# Reshape bgmCompare's stacked threshold summaries to variable-by-category form.
threshold_table <- function(fit, summary_table) {
  labels <- fit$arguments$data_columnnames
  counts <- fit$arguments$num_categories
  max_count <- max(counts)
  out <- matrix(
    NA_real_, length(labels), max_count,
    dimnames = list(labels, paste0("Cat. ", seq_len(max_count)))
  )
  starts <- cumsum(c(1L, head(counts, -1L)))
  for (i in seq_along(labels)) {
    rows <- seq.int(starts[i], length.out = counts[i])
    out[i, seq_along(rows)] <- summary_table$mean[rows]
  }
  out[abs(out) < 0.005] <- 0
  data.frame(Variable = labels, round(out, 2), check.names = FALSE)
}

# Write the same threshold table in the format included by the LaTeX manuscript.
write_latex_table <- function(x, caption, label, path) {
  tab <- xtable::xtable(
    x,
    caption = caption,
    label = label,
    align = c("l", "l", rep("r", ncol(x) - 1)),
    digits = c(0, 0, rep(2, ncol(x) - 1))
  )
  print(
    tab,
    file = path,
    table.placement = "H",
    include.rownames = FALSE,
    hline.after = c(-1, 0, nrow(x)),
    NA.string = "--",
    caption.placement = "top",
    sanitize.text.function = identity,
    sanitize.colnames.function = function(x) paste0("\\textbf{", x, "}"),
    comment = FALSE
  )
}

# Decode parameter names such as "Withdrawal-Conflict (diff1)" and place the
# corresponding posterior means in a symmetric matrix.
pairwise_summary_matrix <- function(fit, summary_table) {
  labels <- fit$arguments$data_columnnames
  out <- matrix(0, length(labels), length(labels), dimnames = list(labels, labels))
  for (i in seq_len(nrow(summary_table))) {
    parameter <- sub(" \\(diff1\\)$", "", as.character(summary_table$parameter[i]))
    variables <- strsplit(parameter, "-", fixed = TRUE)[[1]]
    out[variables[1], variables[2]] <- summary_table$mean[i]
    out[variables[2], variables[1]] <- summary_table$mean[i]
  }
  out
}

# Reconstruct group thresholds from average +/- one half of the difference.
group_threshold_matrices <- function(fit) {
  labels <- fit$arguments$data_columnnames
  counts <- fit$arguments$num_categories
  baseline <- fit$posterior_mean_main_baseline
  difference <- matrix(
    NA_real_, nrow(baseline), ncol(baseline), dimnames = dimnames(baseline)
  )
  cursor <- 1L
  for (i in seq_along(labels)) {
    rows <- cursor:(cursor + counts[i] - 1L)
    difference[i, seq_len(counts[i])] <-
      fit$posterior_summary_main_differences$mean[rows]
    cursor <- cursor + counts[i]
  }
  list(
    female = baseline - difference / 2,
    male = baseline + difference / 2
  )
}

# Evaluate the ordinal MRF conditional distribution for one target item while
# holding all other items at the supplied values. Category 0 is the reference.
conditional_probabilities <- function(
    target, pairwise, thresholds, conditioning_values, counts) {
  target_index <- match(target, rownames(pairwise))
  conditional_sum <- sum(
    pairwise[target_index, -target_index] * conditioning_values[-target_index]
  )
  eta <- c(
    0,
    vapply(seq_len(counts[target_index]), function(category) {
      thresholds[target_index, category] + 2 * category * conditional_sum
    }, numeric(1))
  )

  # Subtracting max(eta) prevents overflow and cancels during normalization.
  probability <- exp(eta - max(eta))
  probability / sum(probability)
}
