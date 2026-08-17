# Prepare the 2016 microbial data used by downstream analyses.
# Run this file from the hja-synoptic repository root.

library(vegan)

rarefaction_depth <- 10000L
random_seed <- 2016L
output_dir <- "data/derived/microbial_diversity_2016"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the DADA2 counts and sample information.
count_table <- read.csv(
  "results/dada2_2016/asv_count_table.csv",
  check.names = FALSE
)
rownames(count_table) <- count_table$sample_id
count_table$sample_id <- NULL
raw_counts <- as.matrix(count_table)
storage.mode(raw_counts) <- "integer"

metadata <- read.csv(
  "results/dada2_2016/sample_metadata_and_read_tracking.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
if (!setequal(rownames(raw_counts), metadata$sample_id)) {
  stop("The ASV table and sample metadata contain different sample IDs.")
}
metadata <- metadata[match(rownames(raw_counts), metadata$sample_id), ]

# Add the DADA2 quality summaries.
qc <- read.csv(
  "results/dada2_2016/qc/sample_qc_flags.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
qc <- qc[match(metadata$sample_id, qc$sample_id), ]
metadata <- cbind(
  metadata,
  qc[c("filter_retention", "merge_retention", "final_retention", "qc_flag")]
)

# Give the four sample types short habitat names.
habitat_names <- c(
  "planktonic streamwater" = "planktonic",
  "hyporheic water" = "hyporheic",
  "stream sediment" = "sediment",
  "terrestrial soil" = "soil"
)
metadata$habitat <- unname(habitat_names[metadata$sample_type])
if (anyNA(metadata$habitat)) stop("An unexpected sample type was found.")

# Add site-level environmental data to aquatic samples.
environment <- read.csv(
  "data/derived/analysis_inputs/environment/aquatic_59_site_environment.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
aquatic <- metadata$habitat != "soil"
environment_rows <- match(metadata$site_code, environment$site_code)
if (anyNA(environment_rows[aquatic])) {
  stop("At least one aquatic sample does not match the environmental table.")
}
environment_columns <- setdiff(
  names(environment),
  c("site_code", "planktonic_sample_id", "hyporheic_sample_id", "sediment_sample_id")
)
environment_values <- environment[environment_rows, environment_columns]
environment_values[!aquatic, ] <- NA
metadata <- cbind(metadata, environment_values)

# Add the available soil coordinates and source-pool labels.
soil <- read.csv(
  "data/derived/analysis_inputs/source_pools/regional_soil_sample_metadata.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
soil_columns <- c(
  "soil_latitude", "soil_longitude", "soil_elevation_ft",
  "source_pool_role", "comparison_warning"
)
soil_values <- soil[match(metadata$sample_id, soil$sample_id), soil_columns]
metadata <- cbind(metadata, soil_values)

# Record whether an FT-ICR-MS profile is available at each aquatic site.
fticr <- read.csv(
  "data/derived/analysis_inputs/audit/fticr_site_summary.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
fticr_values <- fticr[
  match(metadata$site_code, fticr$site_code),
  setdiff(names(fticr), "site_code")
]
fticr_values[!aquatic, ] <- NA
metadata <- cbind(metadata, fticr_values)
metadata$fticr_profile_available <- aquatic &
  !is.na(metadata$detected_primary_features)

# Summarize the raw sequencing results and select the 10K analysis set.
metadata$sequence_reads <- rowSums(raw_counts)
metadata$observed_asvs <- specnumber(raw_counts)
metadata$singleton_asvs <- rowSums(raw_counts == 1)
metadata$goods_coverage <- 1 - metadata$singleton_asvs / metadata$sequence_reads
metadata$included_10k <- metadata$sequence_reads >= rarefaction_depth

if (!all(metadata$sequence_reads == metadata$nonchim)) {
  stop("ASV counts do not agree with the DADA2 read totals.")
}

unrarefied_counts <- raw_counts[metadata$included_10k, , drop = FALSE]
unrarefied_counts <- unrarefied_counts[
  , colSums(unrarefied_counts) > 0, drop = FALSE
]

set.seed(random_seed)
counts_10k <- rrarefy(unrarefied_counts, sample = rarefaction_depth)
counts_10k <- counts_10k[, colSums(counts_10k) > 0, drop = FALSE]

# ASVs are rows in the written count tables so the CSVs can be opened easily.
unrarefied_table <- data.frame(
  asv_id = colnames(unrarefied_counts),
  t(unrarefied_counts),
  check.names = FALSE
)
counts_10k_table <- data.frame(
  asv_id = colnames(counts_10k),
  t(counts_10k),
  check.names = FALSE
)

write.csv(
  metadata,
  file.path(output_dir, "sample_metadata.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  unrarefied_table,
  file.path(output_dir, "asv_counts_unrarefied.csv"),
  row.names = FALSE
)
write.csv(
  counts_10k_table,
  file.path(output_dir, "asv_counts_10k.csv"),
  row.names = FALSE
)

print(table(metadata$habitat, metadata$included_10k))
message("Prepared flat analysis tables in: ", output_dir)
