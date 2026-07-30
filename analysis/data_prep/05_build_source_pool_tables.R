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
source_dir <- file.path(output_root, "source_pools")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

dada2_dir <- under_project_root(prep_config$paths$dada2_output, project_root)
completion_marker <- file.path(dada2_dir, "session_info.txt")
if (!file.exists(completion_marker)) {
  stop(
    "DADA2 completion marker is absent: ", completion_marker,
    ". This script will not read partial DADA2 output.",
    call. = FALSE
  )
}
crosswalk <- read_source_csv(file.path(audit_dir, "sediment_site_crosswalk.csv"))
crosswalk <- crosswalk[crosswalk$include_sediment_multiblock, , drop = FALSE]
manifest <- read_source_csv(under_project_root(
  prep_config$paths$sample_manifest, project_root
))
soil_raw <- read_source_csv(under_project_root(
  prep_config$paths$soil_metadata, project_root
))
asv_data <- read_source_csv(file.path(dada2_dir, "asv_count_table.csv"))
assert_columns(asv_data, "sample_id", "DADA2 ASV table")
assert_unique_key(asv_data$sample_id, "DADA2 ASV table")

paired_source_map <- crosswalk[c(
  "site_code", "sediment_sample_id", "planktonic_sample_id",
  "hyporheic_sample_id", "has_planktonic_sequence", "has_hyporheic_sequence"
)]
paired_source_map$source_comparison_scope <- ifelse(
  paired_source_map$has_planktonic_sequence & paired_source_map$has_hyporheic_sequence,
  "sediment + same-site planktonic + same-site hyporheic",
  ifelse(
    paired_source_map$has_planktonic_sequence,
    "sediment + same-site planktonic",
    ifelse(
      paired_source_map$has_hyporheic_sequence,
      "sediment + same-site hyporheic",
      "sediment only"
    )
  )
)
write_audit_csv(
  paired_source_map,
  file.path(source_dir, "sediment_44_same_site_source_map.csv")
)

aquatic_sample_ids <- unique(c(
  paired_source_map$sediment_sample_id,
  paired_source_map$planktonic_sample_id[!is.na(paired_source_map$planktonic_sample_id)],
  paired_source_map$hyporheic_sample_id[!is.na(paired_source_map$hyporheic_sample_id)]
))
aquatic_manifest <- manifest[match(aquatic_sample_ids, manifest$sample_id), , drop = FALSE]
aquatic_manifest <- aquatic_manifest[order(
  match(aquatic_manifest$site_code, paired_source_map$site_code),
  match(
    aquatic_manifest$sample_type,
    c("stream sediment", "planktonic streamwater", "hyporheic water")
  )
), , drop = FALSE]
aquatic_ids <- aquatic_manifest$sample_id
aquatic_asv_rows <- match(aquatic_ids, asv_data$sample_id)
aquatic_asv <- as.matrix(
  asv_data[aquatic_asv_rows, setdiff(names(asv_data), "sample_id"), drop = FALSE]
)
storage.mode(aquatic_asv) <- "numeric"
rownames(aquatic_asv) <- aquatic_ids

soil_manifest <- manifest[
  manifest$include & !manifest$is_control & manifest$sample_type == "terrestrial soil",
  ,
  drop = FALSE
]
soil_asv_rows <- match(soil_manifest$sample_id, asv_data$sample_id)
soil_asv <- as.matrix(
  asv_data[soil_asv_rows, setdiff(names(asv_data), "sample_id"), drop = FALSE]
)
storage.mode(soil_asv) <- "numeric"
rownames(soil_asv) <- soil_manifest$sample_id

soil_metadata <- soil_manifest
soil_match <- match(soil_metadata$site_code, soil_raw[["Site.Code"]])
for (column in names(soil_raw)) {
  soil_metadata[[paste0("soil_", column)]] <- soil_raw[[column]][soil_match]
}
soil_metadata$source_pool_role <- ifelse(
  soil_metadata$mapping_status == "MATCHED",
  "regional terrestrial source pool",
  "regional terrestrial source pool; location unresolved"
)
soil_metadata$comparison_warning <-
  "Regional soil context only; not spatially paired to aquatic site_code."

write_audit_csv(
  aquatic_manifest,
  file.path(source_dir, "same_site_aquatic_sample_metadata.csv")
)
write_audit_csv(
  soil_metadata,
  file.path(source_dir, "regional_soil_sample_metadata.csv")
)
saveRDS(aquatic_asv, file.path(source_dir, "same_site_aquatic_asv_counts.rds"))
saveRDS(soil_asv, file.path(source_dir, "regional_soil_asv_counts.rds"))
write_matrix_csv_gz(
  aquatic_asv, "sample_id",
  file.path(source_dir, "same_site_aquatic_asv_counts.csv.gz")
)
write_matrix_csv_gz(
  soil_asv, "sample_id",
  file.path(source_dir, "regional_soil_asv_counts.csv.gz")
)

summary <- data.frame(
  set = c(
    "paired sediment sites", "same-site planktonic sources",
    "same-site hyporheic sources", "sites with both aquatic sources",
    "regional soil samples", "regional soil samples with resolved locations"
  ),
  count = c(
    nrow(paired_source_map),
    sum(paired_source_map$has_planktonic_sequence),
    sum(paired_source_map$has_hyporheic_sequence),
    sum(paired_source_map$has_planktonic_sequence &
          paired_source_map$has_hyporheic_sequence),
    nrow(soil_manifest),
    sum(soil_manifest$mapping_status == "MATCHED")
  ),
  stringsAsFactors = FALSE
)
write_audit_csv(summary, file.path(audit_dir, "source_pool_summary.csv"))
stopifnot(
  !anyNA(aquatic_asv_rows),
  !anyNA(soil_asv_rows),
  nrow(paired_source_map) == 44L,
  sum(paired_source_map$has_planktonic_sequence) == 15L,
  sum(paired_source_map$has_hyporheic_sequence) == 26L,
  sum(paired_source_map$has_planktonic_sequence &
        paired_source_map$has_hyporheic_sequence) == 10L,
  nrow(soil_manifest) == 15L,
  sum(soil_manifest$mapping_status == "MATCHED") == 14L
)
capture_session(file.path(audit_dir, "session_info_source_pools.txt"))
message("Source-pool tables complete; no colonization direction was inferred.")
