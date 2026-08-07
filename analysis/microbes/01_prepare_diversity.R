#!/usr/bin/env Rscript

# Align the DADA2 ASV table with sample, environmental, soil, and FT-ICR-MS
# metadata. Raw inputs are never modified.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "helpers.R"))
source(file.path(script_dir, "config.R"))
args <- parse_microbe_args(commandArgs(trailingOnly = TRUE), script_dir)
project_root <- normalizePath(args$project_root, mustWork = TRUE)
output_root <- if (is.na(args$output_root)) {
  under_project_root(microbe_config$output, project_root)
} else {
  under_project_root(args$output_root, project_root)
}
derived_dir <- file.path(output_root, "derived")
audit_dir <- file.path(output_root, "audit")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

dada2_dir <- under_project_root(microbe_config$input$dada2, project_root)
input_dir <- under_project_root(microbe_config$input$analysis_inputs, project_root)
paths <- c(
  counts = file.path(dada2_dir, "sequence_table_asv.rds"),
  sequences = file.path(dada2_dir, "asv_sequences.csv"),
  taxonomy = file.path(dada2_dir, "asv_taxonomy.csv"),
  sample_metadata = file.path(dada2_dir, "sample_metadata_and_read_tracking.csv"),
  qc = file.path(dada2_dir, "qc", "sample_qc_flags.csv"),
  aquatic_environment = file.path(
    input_dir, "environment", "aquatic_59_site_environment.csv"
  ),
  soil_metadata = file.path(
    input_dir, "source_pools", "regional_soil_sample_metadata.csv"
  ),
  fticr_summary = file.path(input_dir, "audit", "fticr_site_summary.csv")
)
assert_files(paths, "Diversity input")

counts <- readRDS(paths[["counts"]])
sequences <- read_analysis_csv(paths[["sequences"]])
taxonomy <- read_analysis_csv(paths[["taxonomy"]])
metadata <- read_analysis_csv(paths[["sample_metadata"]])
qc <- read_analysis_csv(paths[["qc"]])
aquatic_environment <- read_analysis_csv(paths[["aquatic_environment"]])
soil_metadata <- read_analysis_csv(paths[["soil_metadata"]])
fticr_summary <- read_analysis_csv(paths[["fticr_summary"]])

assert_unique(rownames(counts), "DADA2 sample IDs")
assert_unique(sequences$asv_id, "ASV sequence IDs")
assert_unique(taxonomy$asv_id, "ASV taxonomy IDs")
assert_unique(metadata$sample_id, "Sample metadata IDs")
assert_unique(aquatic_environment$site_code, "Aquatic site IDs")
assert_unique(soil_metadata$sample_id, "Soil sample IDs")
assert_unique(fticr_summary$site_code, "FT-ICR-MS site IDs")

metadata <- metadata[match(rownames(counts), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id)) stop("Some ASV-table samples lack metadata.")

sequence_rows <- match(colnames(counts), sequences$sequence)
if (anyNA(sequence_rows)) stop("Some ASV sequences lack ASV identifiers.")
sequences <- sequences[sequence_rows, , drop = FALSE]
colnames(counts) <- sequences$asv_id
taxonomy <- taxonomy[match(colnames(counts), taxonomy$asv_id), , drop = FALSE]
if (anyNA(taxonomy$asv_id)) stop("Some ASVs lack taxonomy rows.")

qc <- qc[match(metadata$sample_id, qc$sample_id), , drop = FALSE]
if (anyNA(qc$sample_id)) stop("Some samples lack DADA2 QC rows.")
metadata$filter_retention <- qc$filter_retention
metadata$merge_retention <- qc$merge_retention
metadata$final_retention <- qc$final_retention
metadata$qc_flag <- qc$qc_flag

habitat_lookup <- c(
  "planktonic streamwater" = "planktonic",
  "hyporheic water" = "hyporheic",
  "stream sediment" = "sediment",
  "terrestrial soil" = "soil"
)
metadata$habitat <- unname(habitat_lookup[metadata$sample_type])
if (anyNA(metadata$habitat)) stop("Unexpected sample type in DADA2 metadata.")
metadata$habitat <- factor(metadata$habitat, levels = microbe_config$habitat_order)
metadata$habitat_label <- unname(
  microbe_config$habitat_labels[as.character(metadata$habitat)]
)

aquatic <- metadata$habitat != "soil"
environment_rows <- match(metadata$site_code, aquatic_environment$site_code)
environment_columns <- setdiff(
  names(aquatic_environment),
  c("site_code", "planktonic_sample_id", "hyporheic_sample_id", "sediment_sample_id")
)
for (column in environment_columns) {
  metadata[[column]] <- aquatic_environment[[column]][environment_rows]
  metadata[[column]][!aquatic] <- NA
}

soil_rows <- match(metadata$sample_id, soil_metadata$sample_id)
for (column in c(
  "soil_latitude", "soil_longitude", "soil_elevation_ft",
  "source_pool_role", "comparison_warning"
)) {
  metadata[[column]] <- soil_metadata[[column]][soil_rows]
}

fticr_rows <- match(metadata$site_code, fticr_summary$site_code)
for (column in setdiff(names(fticr_summary), "site_code")) {
  metadata[[column]] <- fticr_summary[[column]][fticr_rows]
  metadata[[column]][!aquatic] <- NA
}

metadata$environment_join_status <- ifelse(
  aquatic,
  ifelse(is.na(environment_rows), "UNMATCHED_AQUATIC_MASTER", "MATCHED_AQUATIC_MASTER"),
  ifelse(is.na(soil_rows), "SOIL_LOCATION_UNRESOLVED", "MATCHED_SOIL_METADATA")
)
if (any(metadata$environment_join_status == "UNMATCHED_AQUATIC_MASTER")) {
  stop("At least one aquatic sample failed the environmental-data join.")
}
metadata$fticr_profile_available <- aquatic &
  !is.na(metadata$detected_primary_features)
metadata$meets_rarefaction_depth <- metadata$nonchim >=
  microbe_config$rarefaction_depth

metadata_csv <- metadata
metadata_csv$habitat <- as.character(metadata_csv$habitat)
write_analysis_csv(
  metadata_csv,
  file.path(derived_dir, "analysis_sample_metadata.csv")
)
write_analysis_csv(sequences, file.path(derived_dir, "asv_sequence_metadata.csv"))
write_analysis_csv(taxonomy, file.path(derived_dir, "asv_taxonomy.csv"))

saveRDS(
  list(
    asv_counts = counts,
    sample_metadata = metadata,
    asv_sequences = sequences,
    taxonomy = taxonomy,
    parameters = microbe_config
  ),
  file.path(derived_dir, "diversity_input.rds"),
  compress = "xz"
)

join_audit <- metadata_csv[c(
  "sample_id", "site_code_raw", "site_code", "sample_type", "habitat",
  "mapping_status", "normalization_note", "environment_join_status",
  "fticr_profile_available", "qc_flag", "nonchim", "meets_rarefaction_depth"
)]
write_analysis_csv(join_audit, file.path(audit_dir, "sample_join_audit.csv"))

analysis_set_summary <- aggregate(
  cbind(
    libraries = rep(1L, nrow(metadata)),
    retained_at_rarefaction_depth = metadata$meets_rarefaction_depth
  ),
  list(habitat = as.character(metadata$habitat)),
  sum
)
analysis_set_summary <- rbind(
  analysis_set_summary,
  data.frame(
    habitat = "TOTAL",
    libraries = nrow(metadata),
    retained_at_rarefaction_depth = sum(metadata$meets_rarefaction_depth)
  )
)
write_analysis_csv(
  analysis_set_summary,
  file.path(audit_dir, "analysis_set_summary.csv")
)

stopifnot(
  nrow(counts) == 120L,
  ncol(counts) == 76020L,
  all(rowSums(counts) == metadata$nonchim),
  sum(metadata$habitat == "planktonic") == 24L,
  sum(metadata$habitat == "hyporheic") == 36L,
  sum(metadata$habitat == "sediment") == 45L,
  sum(metadata$habitat == "soil") == 15L
)
capture_analysis_session(file.path(audit_dir, "session_info_prepare.txt"))
message("Prepared microbial analysis inputs: ", normalizePath(output_root))
