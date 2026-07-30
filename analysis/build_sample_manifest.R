#!/usr/bin/env Rscript

# Reconcile the 2016 HJA FASTQ files with the local sample list and
# environmental metadata. Source files are read only; all outputs are written
# to data/derived (or the directory supplied with --out-dir).

parse_args <- function(args) {
  values <- list(
    project_root = ".",
    sequence_dir = "sequences",
    sample_list = "data/hja-synoptic_sequence-sample-list.csv",
    aquatic_metadata = "data/hja-env_data_clean.csv",
    soil_metadata = "data/hja-synoptic_env-data-soils.csv",
    out_dir = "data/derived"
  )
  flags <- c(
    "--project-root" = "project_root",
    "--sequence-dir" = "sequence_dir",
    "--sample-list" = "sample_list",
    "--aquatic-metadata" = "aquatic_metadata",
    "--soil-metadata" = "soil_metadata",
    "--out-dir" = "out_dir"
  )
  i <- 1L
  while (i <= length(args)) {
    flag <- args[[i]]
    if (!flag %in% names(flags) || i == length(args)) {
      stop("Unknown or incomplete argument: ", flag, call. = FALSE)
    }
    values[[flags[[flag]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values
}

is_absolute <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

under_root <- function(path, root) {
  if (is_absolute(path)) path else file.path(root, path)
}

read_csv_preserve <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN"),
    strip.white = FALSE
  )
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(args$project_root, mustWork = TRUE)
sequence_dir <- under_root(args$sequence_dir, root)
sample_list_path <- under_root(args$sample_list, root)
aquatic_path <- under_root(args$aquatic_metadata, root)
soil_path <- under_root(args$soil_metadata, root)
out_dir <- under_root(args$out_dir, root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(sample_list_path, aquatic_path, soil_path)
missing_inputs <- required_files[!file.exists(required_files)]
if (length(missing_inputs)) {
  stop("Missing input file(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}
if (!dir.exists(sequence_dir)) {
  stop("Missing sequence directory: ", sequence_dir, call. = FALSE)
}

fastq_paths <- sort(list.files(
  sequence_dir,
  pattern = "_R[12]_001[.]fastq[.]gz$",
  full.names = TRUE
))
if (!length(fastq_paths)) {
  stop("No paired FASTQ files found in ", sequence_dir, call. = FALSE)
}

fastq_names <- basename(fastq_paths)
sample_id <- sub(
  ".*(hja2016_[0-9]{3}).*",
  "\\1",
  fastq_names
)
read_direction <- sub(
  ".*_(R[12])_001[.]fastq[.]gz$",
  "\\1",
  fastq_names
)
plate <- sub(
  ".*-(Plate[0-9]+)-.*",
  "\\1",
  fastq_names
)
lane_sample <- sub(
  ".*_(S[0-9]+)_L[0-9]+_R[12]_001[.]fastq[.]gz$",
  "\\1",
  fastq_names
)

if (any(sample_id == fastq_names) || any(!read_direction %in% c("R1", "R2"))) {
  stop("At least one FASTQ filename did not match the expected naming convention.", call. = FALSE)
}

file_index <- data.frame(
  sample_id = sample_id,
  read_direction = read_direction,
  plate = plate,
  lane_sample = lane_sample,
  filename = fastq_names,
  stringsAsFactors = FALSE
)
duplicate_file_keys <- duplicated(file_index[c("sample_id", "read_direction")]) |
  duplicated(file_index[c("sample_id", "read_direction")], fromLast = TRUE)
if (any(duplicate_file_keys)) {
  stop(
    "Duplicate sample/read FASTQ keys: ",
    paste(
      apply(file_index[duplicate_file_keys, c("sample_id", "read_direction")], 1, paste, collapse = "/"),
      collapse = ", "
    ),
    call. = FALSE
  )
}

forward <- file_index[file_index$read_direction == "R1", c(
  "sample_id", "plate", "lane_sample", "filename"
)]
names(forward)[names(forward) == "filename"] <- "forward_filename"
reverse <- file_index[file_index$read_direction == "R2", c(
  "sample_id", "plate", "lane_sample", "filename"
)]
names(reverse)[names(reverse) == "filename"] <- "reverse_filename"
pairs <- merge(
  forward,
  reverse,
  by = c("sample_id", "plate", "lane_sample"),
  all = TRUE,
  sort = TRUE
)
if (anyNA(pairs$forward_filename) || anyNA(pairs$reverse_filename)) {
  stop("At least one library is missing its R1 or R2 file.", call. = FALSE)
}
pairs$run_id <- paste0("HJA2016_", pairs$plate)
pairs$forward_path <- file.path(args$sequence_dir, pairs$forward_filename)
pairs$reverse_path <- file.path(args$sequence_dir, pairs$reverse_filename)

sample_list <- read_csv_preserve(sample_list_path)
required_sample_columns <- c("Sample Name", "Site Code", "Sample Type")
if (!all(required_sample_columns %in% names(sample_list))) {
  stop("Sample list is missing required columns.", call. = FALSE)
}
sample_list <- sample_list[required_sample_columns]
names(sample_list) <- c("sample_id", "site_code_raw", "sample_type_original")
sample_list$sample_id <- clean_text(sample_list$sample_id)
sample_list$site_code_raw <- clean_text(sample_list$site_code_raw)
sample_list$sample_type_original <- clean_text(sample_list$sample_type_original)
if (anyDuplicated(sample_list$sample_id)) {
  stop("Duplicate sample IDs in the sample list.", call. = FALSE)
}

crosswalk <- merge(pairs, sample_list, by = "sample_id", all.x = TRUE, sort = TRUE)
crosswalk$site_code <- crosswalk$site_code_raw
crosswalk$site_code[crosswalk$site_code_raw == "WS1_S8"] <- "WS1-8"
crosswalk$normalization_note <- NA_character_
crosswalk$normalization_note[crosswalk$site_code_raw == "WS1_S8"] <-
  "Normalized WS1_S8 to master-data site code WS1-8."

# hja2016_200 is missing from both the CSV and the original Numbers sheet.
# The public archive identifies it only as a soil metagenome (SAMN12129905).
unresolved_200 <- crosswalk$sample_id == "hja2016_200" &
  is.na(crosswalk$sample_type_original)
crosswalk$sample_type_original[unresolved_200] <- "Terrestrial soil"
crosswalk$normalization_note[unresolved_200] <-
  paste(
    "Missing from local sample list; NCBI BioSample SAMN12129905 identifies",
    "a soil metagenome, but the soil site is unresolved."
  )

sample_type_map <- c(
  "Stream" = "planktonic streamwater",
  "Hyporheic" = "hyporheic water",
  "Sediment" = "stream sediment",
  "Terrestrial soil" = "terrestrial soil"
)
crosswalk$sample_type <- unname(sample_type_map[crosswalk$sample_type_original])
crosswalk$sample_year <- 2016L
crosswalk$include <- TRUE
crosswalk$is_control <- FALSE
crosswalk$control_type <- NA_character_
crosswalk$biosample_accession <- NA_character_
crosswalk$sra_run_accession <- NA_character_
crosswalk$biosample_accession[crosswalk$sample_id == "hja2016_200"] <- "SAMN12129905"
crosswalk$sra_run_accession[crosswalk$sample_id == "hja2016_200"] <- "SRR9592656"

aquatic <- read_csv_preserve(aquatic_path)
if (!all(c("Site Code", "Sample Type") %in% names(aquatic))) {
  stop("Aquatic metadata is missing Site Code or Sample Type.", call. = FALSE)
}
aquatic$`Site Code` <- clean_text(aquatic$`Site Code`)
aquatic$`Sample Type` <- clean_text(aquatic$`Sample Type`)
aquatic$.metadata_row <- seq_len(nrow(aquatic)) + 1L
aquatic <- aquatic[
  !is.na(aquatic$`Site Code`) & !is.na(aquatic$`Sample Type`),
  ,
  drop = FALSE
]
aquatic_key <- paste(aquatic$`Site Code`, aquatic$`Sample Type`, sep = "\r")
if (anyDuplicated(aquatic_key)) {
  stop("Aquatic metadata has duplicate Site Code / Sample Type keys.", call. = FALSE)
}

soils <- read_csv_preserve(soil_path)
if (!"Site.Code" %in% names(soils)) {
  stop("Soil metadata is missing Site.Code.", call. = FALSE)
}
soils$Site.Code <- clean_text(soils$Site.Code)
soils$.metadata_row <- seq_len(nrow(soils)) + 1L
if (anyDuplicated(soils$Site.Code[!is.na(soils$Site.Code)])) {
  stop("Soil metadata has duplicate Site.Code values.", call. = FALSE)
}

crosswalk$metadata_source <- ifelse(
  crosswalk$sample_type_original == "Terrestrial soil",
  "data/hja-synoptic_env-data-soils.csv",
  "data/hja-env_data_clean.csv"
)
crosswalk$metadata_row <- NA_integer_
crosswalk$mapping_status <- "REVIEW_NO_METADATA_MATCH"

aquatic_rows <- crosswalk$sample_type_original %in% c("Stream", "Hyporheic", "Sediment")
aquatic_match <- match(
  paste(
    crosswalk$site_code[aquatic_rows],
    crosswalk$sample_type_original[aquatic_rows],
    sep = "\r"
  ),
  aquatic_key
)
crosswalk$metadata_row[aquatic_rows] <- aquatic$.metadata_row[aquatic_match]
matched_aquatic <- aquatic_rows
matched_aquatic[aquatic_rows] <- !is.na(aquatic_match)
crosswalk$mapping_status[matched_aquatic] <- "MATCHED"
normalized_match <- matched_aquatic &
  !is.na(crosswalk$site_code_raw) &
  crosswalk$site_code_raw != crosswalk$site_code
crosswalk$mapping_status[normalized_match] <- "MATCHED_NORMALIZED_SITE_CODE"

soil_rows <- crosswalk$sample_type_original == "Terrestrial soil"
soil_match <- match(crosswalk$site_code[soil_rows], soils$Site.Code)
crosswalk$metadata_row[soil_rows] <- soils$.metadata_row[soil_match]
matched_soils <- soil_rows
matched_soils[soil_rows] <- !is.na(soil_match)
crosswalk$mapping_status[matched_soils] <- "MATCHED"
crosswalk$mapping_status[unresolved_200] <- "REVIEW_UNRESOLVED_SITE"

manifest_columns <- c(
  "sample_id", "sample_year", "run_id", "plate", "lane_sample",
  "site_code_raw", "site_code", "sample_type_original", "sample_type",
  "forward_filename", "reverse_filename", "forward_path", "reverse_path",
  "include", "is_control", "control_type", "metadata_source", "metadata_row",
  "mapping_status", "normalization_note", "biosample_accession",
  "sra_run_accession"
)
crosswalk <- crosswalk[order(crosswalk$sample_id), manifest_columns]
aquatic_rows <- crosswalk$sample_type_original %in% c("Stream", "Hyporheic", "Sediment")
soil_rows <- crosswalk$sample_type_original == "Terrestrial soil"
normalized_match <- crosswalk$mapping_status == "MATCHED_NORMALIZED_SITE_CODE"
unresolved_200 <- crosswalk$sample_id == "hja2016_200"

aquatic_crosswalk <- crosswalk[aquatic_rows, , drop = FALSE]
aq_idx <- match(
  paste(
    aquatic_crosswalk$site_code,
    aquatic_crosswalk$sample_type_original,
    sep = "\r"
  ),
  aquatic_key
)
aquatic_prefixed <- aquatic[aq_idx, , drop = FALSE]
names(aquatic_prefixed) <- paste0("env_", names(aquatic_prefixed))
aquatic_joined <- cbind(aquatic_crosswalk, aquatic_prefixed)

soil_crosswalk <- crosswalk[soil_rows, , drop = FALSE]
soil_idx <- match(soil_crosswalk$site_code, soils$Site.Code)
soil_prefixed <- soils[soil_idx, , drop = FALSE]
names(soil_prefixed) <- paste0("soil_", names(soil_prefixed))
soil_joined <- cbind(soil_crosswalk, soil_prefixed)

aquatic_types <- c("Stream", "Hyporheic", "Sediment")
aquatic_sites <- sort(unique(aquatic_crosswalk$site_code))
site_completeness <- data.frame(
  site_code = aquatic_sites,
  planktonic_streamwater = FALSE,
  hyporheic_water = FALSE,
  stream_sediment = FALSE,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(site_completeness))) {
  types <- aquatic_crosswalk$sample_type_original[
    aquatic_crosswalk$site_code == site_completeness$site_code[[i]]
  ]
  site_completeness$planktonic_streamwater[[i]] <- "Stream" %in% types
  site_completeness$hyporheic_water[[i]] <- "Hyporheic" %in% types
  site_completeness$stream_sediment[[i]] <- "Sediment" %in% types
}
presence_columns <- c(
  "planktonic_streamwater", "hyporheic_water", "stream_sediment"
)
site_completeness$sample_type_count <- rowSums(site_completeness[presence_columns])
site_completeness$complete_triplet <- site_completeness$sample_type_count == 3L
site_completeness$missing_sample_types <- apply(
  site_completeness[presence_columns],
  1,
  function(present) {
    labels <- c("planktonic streamwater", "hyporheic water", "stream sediment")
    paste(labels[!present], collapse = "; ")
  }
)

issues <- data.frame(
  issue_type = character(),
  sample_id = character(),
  site_code = character(),
  detail = character(),
  action = character(),
  stringsAsFactors = FALSE
)
add_issue <- function(issue_type, sample_id = NA_character_, site_code = NA_character_,
                      detail, action) {
  data.frame(
    issue_type = issue_type,
    sample_id = sample_id,
    site_code = site_code,
    detail = detail,
    action = action,
    stringsAsFactors = FALSE
  )
}
if (any(normalized_match)) {
  issues <- rbind(
    issues,
    add_issue(
      "NORMALIZED_SITE_CODE",
      crosswalk$sample_id[normalized_match],
      crosswalk$site_code[normalized_match],
      "Sequence sample list used WS1_S8; master metadata uses WS1-8.",
      "Automatically normalized to WS1-8."
    )
  )
}
if (any(unresolved_200)) {
  issues <- rbind(
    issues,
    add_issue(
      "UNRESOLVED_SOIL_SITE",
      "hja2016_200",
      NA_character_,
      paste(
        "The local CSV and original Numbers sheet omit this library.",
        "NCBI SAMN12129905 confirms soil metagenome but does not identify the site."
      ),
      "Keep in DADA2; do not attach site-level environmental metadata until identity is confirmed."
    )
  )
}
soil_sites_with_sequences <- unique(soil_crosswalk$site_code[!is.na(soil_crosswalk$site_code)])
soil_without_sequences <- setdiff(soils$Site.Code, soil_sites_with_sequences)
if (length(soil_without_sequences)) {
  issues <- rbind(
    issues,
    add_issue(
      "SOIL_METADATA_WITHOUT_SEQUENCE",
      NA_character_,
      soil_without_sequences,
      "Soil metadata row has no mapped FASTQ library.",
      "Retain metadata; exclude from sequence analyses unless a library is identified."
    )
  )
}
incomplete_sites <- site_completeness[!site_completeness$complete_triplet, , drop = FALSE]
if (nrow(incomplete_sites)) {
  issues <- rbind(
    issues,
    data.frame(
      issue_type = "INCOMPLETE_AQUATIC_TRIPLET",
      sample_id = NA_character_,
      site_code = incomplete_sites$site_code,
      detail = paste0("Missing: ", incomplete_sites$missing_sample_types),
      action = "Treat the design as unbalanced; do not synthesize missing libraries.",
      stringsAsFactors = FALSE
    )
  )
}

summary <- data.frame(
  metric = c(
    "FASTQ files", "paired libraries", "sample-list rows",
    "mapped without normalization", "mapped after site-code normalization",
    "unresolved libraries", "planktonic streamwater libraries",
    "hyporheic-water libraries", "stream-sediment libraries",
    "terrestrial-soil libraries", "aquatic sites represented",
    "complete aquatic triplets", "incomplete aquatic sites",
    "soil metadata rows without sequences"
  ),
  value = c(
    nrow(file_index),
    nrow(crosswalk),
    nrow(sample_list),
    sum(crosswalk$mapping_status == "MATCHED"),
    sum(crosswalk$mapping_status == "MATCHED_NORMALIZED_SITE_CODE"),
    sum(grepl("^REVIEW", crosswalk$mapping_status)),
    sum(crosswalk$sample_type == "planktonic streamwater"),
    sum(crosswalk$sample_type == "hyporheic water"),
    sum(crosswalk$sample_type == "stream sediment"),
    sum(crosswalk$sample_type == "terrestrial soil"),
    nrow(site_completeness),
    sum(site_completeness$complete_triplet),
    sum(!site_completeness$complete_triplet),
    length(soil_without_sequences)
  ),
  stringsAsFactors = FALSE
)

write_csv(crosswalk, file.path(out_dir, "hja_2016_sample_manifest.csv"))
write_csv(
  aquatic_joined,
  file.path(out_dir, "hja_2016_aquatic_sequence_metadata.csv")
)
write_csv(
  soil_joined,
  file.path(out_dir, "hja_2016_soil_sequence_metadata.csv")
)
write_csv(
  site_completeness,
  file.path(out_dir, "hja_2016_site_completeness.csv")
)
write_csv(issues, file.path(out_dir, "hja_2016_mapping_issues.csv"))
write_csv(summary, file.path(out_dir, "hja_2016_mapping_summary.csv"))

message("Wrote 2016 mapping outputs to: ", normalizePath(out_dir))
print(summary, row.names = FALSE)

if (nrow(crosswalk) != 120L ||
    sum(crosswalk$mapping_status == "MATCHED") != 118L ||
    sum(crosswalk$mapping_status == "MATCHED_NORMALIZED_SITE_CODE") != 1L ||
    sum(crosswalk$mapping_status == "REVIEW_UNRESOLVED_SITE") != 1L) {
  stop("Mapping totals differ from the audited expected values.", call. = FALSE)
}
