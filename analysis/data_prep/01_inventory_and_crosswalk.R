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
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the R package 'readxl' to read the FT-ICR-MS workbook.", call. = FALSE)
}

source_paths <- unlist(prep_config$paths[c(
  "master_environment", "fticr_workbook", "sample_manifest", "soil_metadata"
)])
inventory <- file_inventory(source_paths, project_root)
write_audit_csv(inventory, file.path(audit_dir, "source_file_inventory.csv"))

ft_path <- under_project_root(prep_config$paths$fticr_workbook, project_root)
ft_header <- suppressMessages(readxl::read_xlsx(
  ft_path, n_max = 0, .name_repair = "minimal"
))
ft_names <- names(ft_header)
if (length(ft_names) != 76L) {
  stop("Expected 76 FT-ICR-MS workbook columns, found ", length(ft_names), ".", call. = FALSE)
}
if (!identical(
  ft_names[prep_config$fticr$metadata_columns],
  prep_config$fticr$expected_metadata_headers
)) {
  stop("FT-ICR-MS metadata headers or positions differ from config.R.", call. = FALSE)
}

column_role <- rep(NA_character_, length(ft_names))
column_role[prep_config$fticr$metadata_columns] <- "molecular_metadata"
column_role[prep_config$fticr$site_intensity_columns] <- "site_intensity"
column_role[prep_config$fticr$lab_standard_columns] <- "lab_standard"
derived_identifier <- ft_names
derived_identifier[prep_config$fticr$metadata_columns] <-
  prep_config$fticr$metadata_analysis_names
derived_identifier[prep_config$fticr$lab_standard_columns] <- c("LAB_QC_1", "LAB_QC_2")
column_dictionary <- data.frame(
  source_workbook_column = seq_along(ft_names),
  source_header = ft_names,
  column_role = column_role,
  derived_identifier = derived_identifier,
  included_in_ecological_matrix = column_role == "site_intensity",
  mapping_rule = ifelse(
    column_role == "site_intensity",
    "exact source header to site_code",
    "declared position in config.R"
  ),
  stringsAsFactors = FALSE
)
write_audit_csv(
  column_dictionary,
  file.path(audit_dir, "fticr_column_dictionary.csv")
)
write_audit_csv(
  environment_variable_map,
  file.path(audit_dir, "environment_variable_dictionary.csv")
)

master <- read_source_csv(under_project_root(
  prep_config$paths$master_environment, project_root
))
assert_columns(master, c("Site Code", "Sample Type"), "master environmental CSV")
master$.source_workbook_row <- seq_len(nrow(master)) + 1L
master_records <- master[
  !is.na(master[["Site Code"]]) & master[["Site Code"]] != "" &
    master[["Sample Type"]] %in% c("Stream", "Hyporheic", "Sediment"),
  ,
  drop = FALSE
]
master_key <- paste(master_records[["Site Code"]], master_records[["Sample Type"]], sep = "::")
assert_unique_key(master_key, "master environmental site/sample-type rows")

manifest <- read_source_csv(under_project_root(
  prep_config$paths$sample_manifest, project_root
))
assert_columns(
  manifest,
  c("sample_id", "site_code", "sample_type", "include", "is_control"),
  "2016 sample manifest"
)
manifest <- manifest[manifest$include & !manifest$is_control, , drop = FALSE]
manifest_key <- paste(manifest$site_code, manifest$sample_type, sep = "::")
assert_unique_key(manifest_key, "included non-control sample manifest")

sample_id_for <- function(site_code, sample_type) {
  one_value_or_na(manifest$sample_id[
    manifest$site_code == site_code & manifest$sample_type == sample_type
  ])
}
master_row_for <- function(site_code, sample_type) {
  value <- master_records$.source_workbook_row[
    master_records[["Site Code"]] == site_code &
      master_records[["Sample Type"]] == sample_type
  ]
  if (!length(value)) return(NA_integer_)
  if (length(value) > 1L) stop("Duplicate master match for ", site_code, " / ", sample_type)
  value
}

ft_site_codes <- ft_names[prep_config$fticr$site_intensity_columns]
assert_unique_key(ft_site_codes, "FT-ICR-MS site headers")
crosswalk <- data.frame(
  fticr_source_column = prep_config$fticr$site_intensity_columns,
  site_code = ft_site_codes,
  sediment_sample_id = vapply(
    ft_site_codes, sample_id_for, sample_type = "stream sediment",
    FUN.VALUE = character(1)
  ),
  planktonic_sample_id = vapply(
    ft_site_codes, sample_id_for, sample_type = "planktonic streamwater",
    FUN.VALUE = character(1)
  ),
  hyporheic_sample_id = vapply(
    ft_site_codes, sample_id_for, sample_type = "hyporheic water",
    FUN.VALUE = character(1)
  ),
  stream_master_row = vapply(
    ft_site_codes, master_row_for, sample_type = "Stream",
    FUN.VALUE = integer(1)
  ),
  hyporheic_master_row = vapply(
    ft_site_codes, master_row_for, sample_type = "Hyporheic",
    FUN.VALUE = integer(1)
  ),
  sediment_master_row = vapply(
    ft_site_codes, master_row_for, sample_type = "Sediment",
    FUN.VALUE = integer(1)
  ),
  stringsAsFactors = FALSE
)
crosswalk$has_sediment_sequence <- !is.na(crosswalk$sediment_sample_id)
crosswalk$has_planktonic_sequence <- !is.na(crosswalk$planktonic_sample_id)
crosswalk$has_hyporheic_sequence <- !is.na(crosswalk$hyporheic_sample_id)
crosswalk$has_all_master_rows <- with(
  crosswalk,
  !is.na(stream_master_row) & !is.na(hyporheic_master_row) & !is.na(sediment_master_row)
)
crosswalk$include_sediment_multiblock <-
  crosswalk$has_sediment_sequence & crosswalk$has_all_master_rows
crosswalk$mapping_status <- ifelse(
  crosswalk$include_sediment_multiblock,
  "MATCHED_SEDIMENT_FTICR_ENVIRONMENT",
  ifelse(
    !crosswalk$has_sediment_sequence,
    "FTICR_SITE_WITHOUT_SEDIMENT_SEQUENCE",
    "MISSING_MASTER_ENVIRONMENT_ROW"
  )
)
write_audit_csv(crosswalk, file.path(audit_dir, "sediment_site_crosswalk.csv"))

sediment_sites <- manifest$site_code[manifest$sample_type == "stream sediment"]
issues <- rbind(
  data.frame(
    issue_type = "FTICR_SITE_WITHOUT_SEDIMENT_SEQUENCE",
    site_code = crosswalk$site_code[!crosswalk$has_sediment_sequence],
    detail = "Retained in the 60-site FT-ICR table; excluded from the paired sediment table.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    issue_type = "SEDIMENT_SEQUENCE_WITHOUT_FTICR_SITE",
    site_code = setdiff(sediment_sites, ft_site_codes),
    detail = "Retained in microbial-only analyses; excluded from the paired sediment table.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    issue_type = "SOIL_SAMPLE_WITH_UNRESOLVED_LOCATION",
    site_code = manifest$sample_id[
      manifest$sample_type == "terrestrial soil" &
        (is.na(manifest$site_code) | manifest$site_code == "")
    ],
    detail = "Retained as a regional soil source-pool sample with unresolved location.",
    stringsAsFactors = FALSE
  )
)
write_audit_csv(issues, file.path(audit_dir, "crosswalk_issues.csv"))

stopifnot(
  length(ft_site_codes) == 60L,
  sum(crosswalk$include_sediment_multiblock) == 44L,
  sum(crosswalk$include_sediment_multiblock & crosswalk$has_planktonic_sequence) == 15L,
  sum(crosswalk$include_sediment_multiblock & crosswalk$has_hyporheic_sequence) == 26L,
  sum(crosswalk$include_sediment_multiblock &
        crosswalk$has_planktonic_sequence &
        crosswalk$has_hyporheic_sequence) == 10L
)
capture_session(file.path(audit_dir, "session_info_inventory.txt"))
message("Inventory complete: 60 FT-ICR sites; 44 paired sediment sites.")
