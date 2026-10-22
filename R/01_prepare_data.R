# Prepare the analysis data --------------------------------------------------
# Input:  data/raw/PPU Data.xlsx
# Output: one validated six-item RDS file for each comparison group

source(file.path("R", "00_config.R"))
check_packages("readxl")

# Read only once; all subsequent analysis uses the smaller derived files.
raw <- readxl::read_excel(paths$raw_data)
required_columns <- c("gender", names(item_labels))
stopifnot(all(required_columns %in% names(raw)))

# The paper compares respondents recorded as Female or Male. Other gender
# responses remain in the source workbook but are outside this comparison.
analysis_data <- as.data.frame(raw[raw$gender %in% c("Female", "Male"), required_columns])

# These checks make changes to the source data visible before model fitting.
stopifnot(
  nrow(analysis_data) == 907L,
  !anyNA(analysis_data),
  all(vapply(analysis_data[names(item_labels)], function(x) {
    all(x %in% 1:7)
  }, logical(1)))
)

# Split groups and replace questionnaire column codes with figure labels.
female <- analysis_data[analysis_data$gender == "Female", names(item_labels)]
male <- analysis_data[analysis_data$gender == "Male", names(item_labels)]
names(female) <- names(male) <- unname(item_labels)
row.names(female) <- row.names(male) <- NULL

# Store deterministic preprocessing outputs for inspection and reuse.
stopifnot(nrow(female) == 520L, nrow(male) == 387L)
saveRDS(female, paths$female_data, version = 3)
saveRDS(male, paths$male_data, version = 3)

cat("Prepared", nrow(female), "female and", nrow(male), "male observations.\n")
