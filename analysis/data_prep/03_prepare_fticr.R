#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "helpers.R"))
source(file.path(script_dir, "config.R"))
options <- parse_named_args(
  args,
  list(project_root = normalizePath(file.path(script_dir, "..", "..")), output_root = NA_character_),
  c("--project-root" = "project_root", "--output-root" = "output_root")
)
project_root <- normalizePath(options$project_root)
output_root <- if (is.na(options$output_root)) {
  under_project_root(prep_config$paths$output_root, project_root)
} else {
  under_project_root(options$output_root, project_root)
}
audit_dir <- file.path(output_root, "audit")
fticr_dir <- file.path(output_root, "fticr")
dir.create(fticr_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the R package 'readxl' to read the FT-ICR-MS workbook.", call. = FALSE)
}

ft_path <- under_project_root(prep_config$paths$fticr_workbook, project_root)
column_types <- c(
  rep("numeric", 9L), "text", "text", rep("numeric", 3L),
  rep("numeric", 62L)
)
raw <- suppressMessages(readxl::read_xlsx(
  ft_path,
  col_types = column_types,
  .name_repair = "minimal"
))
if (ncol(raw) != 76L) stop("Expected 76 FT-ICR-MS columns.", call. = FALSE)
if (!identical(
  names(raw)[prep_config$fticr$metadata_columns],
  prep_config$fticr$expected_metadata_headers
)) stop("FT-ICR-MS metadata columns do not match config.R.", call. = FALSE)

peak_id <- sprintf("FTICR_%05d", seq_len(nrow(raw)))
metadata <- as.data.frame(
  raw[prep_config$fticr$metadata_columns],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(metadata) <- prep_config$fticr$metadata_analysis_names
metadata <- data.frame(
  peak_id = peak_id,
  source_workbook_row = seq_len(nrow(raw)) + 1L,
  metadata,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

site_by_peak <- t(as.matrix(raw[prep_config$fticr$site_intensity_columns]))
storage.mode(site_by_peak) <- "numeric"
rownames(site_by_peak) <- names(raw)[prep_config$fticr$site_intensity_columns]
colnames(site_by_peak) <- peak_id
qc_by_peak <- t(as.matrix(raw[prep_config$fticr$lab_standard_columns]))
storage.mode(qc_by_peak) <- "numeric"
rownames(qc_by_peak) <- c("LAB_QC_1", "LAB_QC_2")
colnames(qc_by_peak) <- peak_id

filter_values <- prep_config$fticr$primary_filter
detected_site_count <- colSums(site_by_peak > 0, na.rm = TRUE)
within_mass_range <- metadata$measured_mz >= filter_values$minimum_mass &
  metadata$measured_mz <= filter_values$maximum_mass
c13_allowed <- !is.na(metadata$c13_indicator) &
  metadata$c13_indicator <= filter_values$maximum_c13_indicator
detected_in_minimum_sites <-
  detected_site_count >= filter_values$minimum_detected_site_count
formula_assigned <- !is.na(metadata$carbon_count) & metadata$carbon_count > 0
included_primary <- within_mass_range & c13_allowed &
  detected_in_minimum_sites & formula_assigned
exclusion_reason <- vapply(seq_along(peak_id), function(i) {
  reasons <- c(
    if (!within_mass_range[i]) "outside_mass_200_900" else NULL,
    if (!c13_allowed[i]) "c13_indicator_above_zero_or_missing" else NULL,
    if (!detected_in_minimum_sites[i]) "detected_in_fewer_than_2_sites" else NULL,
    if (!formula_assigned[i]) "no_formula_assignment_carbon_count_zero_or_missing" else NULL
  )
  if (!length(reasons)) "" else paste(reasons, collapse = ";")
}, FUN.VALUE = character(1))
filter_status <- data.frame(
  peak_id = peak_id,
  source_workbook_row = metadata$source_workbook_row,
  within_mass_range = within_mass_range,
  c13_allowed = c13_allowed,
  detected_site_count = detected_site_count,
  detected_in_minimum_sites = detected_in_minimum_sites,
  formula_assigned = formula_assigned,
  included_primary = included_primary,
  exclusion_reason = exclusion_reason,
  stringsAsFactors = FALSE
)

write_matrix_csv_gz(
  site_by_peak[, included_primary, drop = FALSE],
  "site_code",
  file.path(fticr_dir, "fticr_site_by_peak_primary.csv.gz")
)
write_matrix_csv_gz(
  qc_by_peak,
  "qc_sample_id",
  file.path(fticr_dir, "fticr_qc_by_peak.csv.gz")
)
metadata_connection <- gzfile(file.path(fticr_dir, "fticr_peak_metadata.csv.gz"), "wt")
write.csv(metadata, metadata_connection, row.names = FALSE, na = "")
close(metadata_connection)
filter_connection <- gzfile(file.path(fticr_dir, "fticr_peak_filter_status.csv.gz"), "wt")
write.csv(filter_status, filter_connection, row.names = FALSE, na = "")
close(filter_connection)
saveRDS(site_by_peak, file.path(fticr_dir, "fticr_site_by_peak.rds"))
saveRDS(
  site_by_peak[, included_primary, drop = FALSE],
  file.path(fticr_dir, "fticr_site_by_peak_primary.rds")
)
saveRDS(qc_by_peak, file.path(fticr_dir, "fticr_qc_by_peak.rds"))

filter_log <- data.frame(
  step = c(
    "raw_features", "mass_200_to_900", "then_c13_indicator_zero",
    "then_detected_in_at_least_2_sites", "then_formula_assigned_carbon_count_gt_zero"
  ),
  retained_features = c(
    nrow(metadata),
    sum(within_mass_range),
    sum(within_mass_range & c13_allowed),
    sum(within_mass_range & c13_allowed & detected_in_minimum_sites),
    sum(included_primary)
  ),
  rule = c(
    "none",
    "measured_mz >= 200 and <= 900",
    "c13_indicator <= 0",
    "site intensity > 0 in at least 2 of 60 site columns",
    "carbon_count > 0"
  ),
  stringsAsFactors = FALSE
)
site_summary <- data.frame(
  site_code = rownames(site_by_peak),
  detected_raw_features = rowSums(site_by_peak > 0, na.rm = TRUE),
  total_raw_intensity = rowSums(site_by_peak, na.rm = TRUE),
  detected_primary_features = rowSums(
    site_by_peak[, included_primary, drop = FALSE] > 0, na.rm = TRUE
  ),
  total_primary_intensity = rowSums(
    site_by_peak[, included_primary, drop = FALSE], na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)
write_audit_csv(filter_log, file.path(audit_dir, "fticr_filter_log.csv"))
write_audit_csv(site_summary, file.path(audit_dir, "fticr_site_summary.csv"))
stopifnot(
  nrow(metadata) == 33741L,
  nrow(site_by_peak) == 60L,
  ncol(site_by_peak) == 33741L,
  sum(!is.na(metadata$elemental_composition_source)) == 9677L,
  sum(included_primary) == 4760L
)
capture_session(file.path(audit_dir, "session_info_fticr.txt"))
message("FT-ICR preparation complete: 33,741 raw features; 4,760 in the primary table.")
