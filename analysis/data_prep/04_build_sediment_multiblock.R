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
multiblock_dir <- file.path(output_root, "sediment_multiblock")
dir.create(multiblock_dir, recursive = TRUE, showWarnings = FALSE)

dada2_dir <- under_project_root(prep_config$paths$dada2_output, project_root)
completion_marker <- file.path(dada2_dir, "session_info.txt")
if (!file.exists(completion_marker)) {
  stop(
    "DADA2 completion marker is absent: ", completion_marker,
    ". This script will not read partial DADA2 output.",
    call. = FALSE
  )
}
required <- c(
  crosswalk = file.path(audit_dir, "sediment_site_crosswalk.csv"),
  environment = file.path(output_root, "environment", "sediment_44_environment.csv"),
  fticr = file.path(output_root, "fticr", "fticr_site_by_peak_primary.rds"),
  asv = file.path(dada2_dir, "asv_count_table.csv"),
  manifest = under_project_root(prep_config$paths$sample_manifest, project_root)
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing prerequisite file(s): ", paste(missing, collapse = ", "))

crosswalk <- read_source_csv(required[["crosswalk"]])
crosswalk <- crosswalk[crosswalk$include_sediment_multiblock, , drop = FALSE]
environment <- read_source_csv(required[["environment"]])
fticr_all <- readRDS(required[["fticr"]])
asv_data <- read_source_csv(required[["asv"]])
manifest <- read_source_csv(required[["manifest"]])
assert_columns(asv_data, "sample_id", "DADA2 ASV table")
assert_unique_key(asv_data$sample_id, "DADA2 ASV table")

site_order <- crosswalk$site_code
sample_order <- crosswalk$sediment_sample_id
environment <- environment[match(site_order, environment$site_code), , drop = FALSE]
fticr <- fticr_all[match(site_order, rownames(fticr_all)), , drop = FALSE]
asv_rows <- match(sample_order, asv_data$sample_id)
asv <- as.matrix(asv_data[asv_rows, setdiff(names(asv_data), "sample_id"), drop = FALSE])
storage.mode(asv) <- "numeric"
rownames(asv) <- site_order
rownames(fticr) <- site_order
fticr_presence_absence <- 1L * (fticr > 0)

manifest_rows <- match(sample_order, manifest$sample_id)
sample_metadata <- manifest[manifest_rows, c(
  "sample_id", "sample_year", "run_id", "plate", "lane_sample",
  "site_code", "sample_type", "metadata_source", "metadata_row",
  "mapping_status", "normalization_note"
), drop = FALSE]
sample_metadata$fticr_source_column <- crosswalk$fticr_source_column
sample_metadata$planktonic_sample_id <- crosswalk$planktonic_sample_id
sample_metadata$hyporheic_sample_id <- crosswalk$hyporheic_sample_id

enzyme_columns <- c(
  "site_code", "sediment_sample_id",
  "sediment_organic_content_percent",
  "sediment_nag_umol_g_hr", "sediment_lap_umol_g_hr",
  "sediment_glu_umol_g_hr", "sediment_ap_umol_g_hr"
)
assert_columns(environment, enzyme_columns, "44-site environmental table")
enzyme_table <- environment[enzyme_columns]

write_audit_csv(
  sample_metadata,
  file.path(multiblock_dir, "sediment_44_sample_metadata.csv")
)
write_audit_csv(
  environment,
  file.path(multiblock_dir, "sediment_44_environment.csv")
)
write_audit_csv(
  enzyme_table,
  file.path(multiblock_dir, "sediment_44_enzyme_activity.csv")
)
saveRDS(asv, file.path(multiblock_dir, "sediment_44_asv_counts.rds"))
saveRDS(fticr, file.path(multiblock_dir, "sediment_44_fticr_intensity.rds"))
saveRDS(
  fticr_presence_absence,
  file.path(multiblock_dir, "sediment_44_fticr_presence_absence.rds")
)
write_matrix_csv_gz(
  asv, "site_code",
  file.path(multiblock_dir, "sediment_44_asv_counts.csv.gz")
)
write_matrix_csv_gz(
  fticr, "site_code",
  file.path(multiblock_dir, "sediment_44_fticr_intensity.csv.gz")
)
write_matrix_csv_gz(
  fticr_presence_absence, "site_code",
  file.path(multiblock_dir, "sediment_44_fticr_presence_absence.csv.gz")
)

alignment <- data.frame(
  row_number = seq_along(site_order),
  site_code = site_order,
  sediment_sample_id = sample_order,
  environment_site_code = environment$site_code,
  asv_source_sample_id = asv_data$sample_id[asv_rows],
  asv_matrix_rowname = rownames(asv),
  fticr_source_site_code = rownames(fticr_all)[match(site_order, rownames(fticr_all))],
  fticr_matrix_rowname = rownames(fticr),
  all_identifiers_aligned = site_order == environment$site_code &
    site_order == rownames(asv) &
    site_order == rownames(fticr) &
    sample_order == asv_data$sample_id[asv_rows],
  stringsAsFactors = FALSE
)
write_audit_csv(alignment, file.path(audit_dir, "sediment_multiblock_alignment.csv"))
stopifnot(
  nrow(asv) == 44L,
  nrow(fticr) == 44L,
  nrow(environment) == 44L,
  all(alignment$all_identifiers_aligned),
  !anyNA(asv_rows),
  !anyNA(match(site_order, rownames(fticr_all)))
)
capture_session(file.path(audit_dir, "session_info_sediment_multiblock.txt"))
message(
  "Sediment multiblock complete: 44 aligned rows x ",
  ncol(asv), " ASVs x ", ncol(fticr), " FT-ICR features."
)
