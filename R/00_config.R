# Shared configuration -------------------------------------------------------
# Sourced by every analysis script. It defines repository-relative paths,
# manuscript labels, dependency checks, and the required output folders.

# Centralizing paths keeps the remaining scripts free of user-specific paths.
paths <- list(
  raw_data = file.path("data", "raw", "PPU Data.xlsx"),
  female_data = file.path("data", "derived", "ppu_female.rds"),
  male_data = file.path("data", "derived", "ppu_male.rds"),
  fit = file.path("results", "model", "bgms_fit.rds"),
  prior_sensitivity = file.path("results", "model", "prior_sensitivity.rds"),
  tables = file.path("results", "tables"),
  diagnostics = file.path("results", "diagnostics"),
  figures = "figures"
)

# Raw column names on the left; concise manuscript labels on the right.
item_labels <- c(
  PPU_1 = "Salience",
  PPU_4 = "Tolerance",
  PPU_6 = "Withdrawal",
  PPU_11 = "Relapse",
  PPU_14 = "Mood Mod",
  PPU_15 = "Conflict"
)

# Fail early with a useful message when a script-specific dependency is absent.
check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the missing package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

# Output folders are safe to create repeatedly and may not exist in a fresh clone.
invisible(lapply(paths[c("tables", "diagnostics", "figures")], dir.create,
  recursive = TRUE, showWarnings = FALSE
))
dir.create(dirname(paths$female_data), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$fit), recursive = TRUE, showWarnings = FALSE)
