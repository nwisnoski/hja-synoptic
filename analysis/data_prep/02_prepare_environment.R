#!/usr/bin/env Rscript

library(here)

# 1. Settings and shared source-column definitions.
source(here("analysis", "data_prep", "config.R"), local = TRUE)
source(here("analysis", "data_prep", "helpers.R"), local = TRUE)
output_root <- prep_config$paths$output_root

# 2. Read inputs and check their identities.
audit_dir <- file.path(output_root, "audit")
environment_dir <- file.path(output_root, "environment")
dir.create(environment_dir, recursive = TRUE, showWarnings = FALSE)

crosswalk_path <- file.path(audit_dir, "sediment_site_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Run 01_inventory_and_crosswalk.R first.", call. = FALSE)
}
crosswalk <- read_source_csv(crosswalk_path)
master <- read_source_csv(prep_config$paths$master_environment)
manifest <- read_source_csv(prep_config$paths$sample_manifest)
assert_columns(master, c("Site Code", "Sample Type"), "master environmental CSV")
assert_columns(
  manifest,
  c("sample_id", "site_code", "sample_type", "include", "is_control"),
  "2016 sample manifest"
)
assert_columns(
  master, unique(environment_variable_map$source_column),
  "master environmental CSV"
)
master$.source_workbook_row <- seq_len(nrow(master)) + 1L

# The same declared column mapping builds both site inventories.
master_key <- paste(master[["Site Code"]], master[["Sample Type"]], sep = "::")
build_environment_table <- function(site_codes) {
  result <- data.frame(site_code = site_codes, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(environment_variable_map))) {
    rule <- environment_variable_map[i, , drop = FALSE]
    target <- paste(site_codes, rule$source_sample_type, sep = "::")
    matched <- match(target, master_key)
    result[[rule$analysis_column]] <- master[[rule$source_column]][matched]
  }
  result
}

aquatic_sample_types <- c(
  "planktonic streamwater", "hyporheic water", "stream sediment"
)
aquatic_manifest <- manifest[
  manifest$include & !manifest$is_control &
    manifest$sample_type %in% aquatic_sample_types,
  ,
  drop = FALSE
]
aquatic_site_codes <- unique(aquatic_manifest$site_code)
# 3. Index the available habitats at each sequenced aquatic site.
# The survey is unbalanced, so absent same-site libraries remain missing.
aquatic_key <- paste(aquatic_manifest$site_code, aquatic_manifest$sample_type, sep = "::")
assert_unique_key(aquatic_key, "included aquatic site/sample-type rows")
aquatic_index <- data.frame(site_code = aquatic_site_codes, stringsAsFactors = FALSE)
sample_columns <- c("planktonic_sample_id", "hyporheic_sample_id", "sediment_sample_id")
for (i in seq_along(aquatic_sample_types)) {
  target_key <- paste(aquatic_site_codes, aquatic_sample_types[i], sep = "::")
  matched <- match(target_key, aquatic_key)
  aquatic_index[[sample_columns[i]]] <- aquatic_manifest$sample_id[matched]
}
aquatic_index$fticr_source_column <- crosswalk$fticr_source_column[
  match(aquatic_index$site_code, crosswalk$site_code)
]
aquatic_index$has_fticr_profile <- !is.na(aquatic_index$fticr_source_column)
# 4. Build environmental tables in the explicit site order.
environment_aquatic <- build_environment_table(aquatic_site_codes)
environment_aquatic <- cbind(
  aquatic_index,
  environment_aquatic[setdiff(names(environment_aquatic), "site_code")],
  stringsAsFactors = FALSE
)

environment_60 <- build_environment_table(crosswalk$site_code)
environment_60 <- cbind(
  crosswalk[c(
    "site_code", "fticr_source_column", "sediment_sample_id",
    "planktonic_sample_id", "hyporheic_sample_id",
    "stream_master_row", "hyporheic_master_row", "sediment_master_row",
    "include_sediment_multiblock"
  )],
  environment_60[setdiff(names(environment_60), "site_code")],
  stringsAsFactors = FALSE
)
environment_44 <- environment_60[
  environment_60$include_sediment_multiblock,
  ,
  drop = FALSE
]
write_audit_csv(
  environment_60,
  file.path(environment_dir, "fticr_60_site_environment.csv")
)
write_audit_csv(
  environment_aquatic,
  file.path(environment_dir, "aquatic_59_site_environment.csv")
)
write_audit_csv(
  environment_44,
  file.path(environment_dir, "sediment_44_environment.csv")
)

# 5. Report missingness for every declared variable in both analysis pools.
missingness_rows <- vector("list", nrow(environment_variable_map))
aquatic_missingness_rows <- vector("list", nrow(environment_variable_map))
for (i in seq_len(nrow(environment_variable_map))) {
  rule <- environment_variable_map[i, , drop = FALSE]
  values <- environment_44[[rule$analysis_column]]
  missingness_rows[[i]] <- data.frame(
    analysis_column = rule$analysis_column,
    source_sample_type = rule$source_sample_type,
    source_column = rule$source_column,
    predictor_block = rule$predictor_block,
    primary_44_site = rule$primary_44_site,
    n_sites = nrow(environment_44),
    n_observed = sum(!is.na(values)),
    n_missing = sum(is.na(values)),
    proportion_missing = mean(is.na(values)),
    stringsAsFactors = FALSE
  )
  values <- environment_aquatic[[rule$analysis_column]]
  aquatic_missingness_rows[[i]] <- data.frame(
    analysis_column = rule$analysis_column,
    source_sample_type = rule$source_sample_type,
    source_column = rule$source_column,
    predictor_block = rule$predictor_block,
    n_sites = nrow(environment_aquatic),
    n_observed = sum(!is.na(values)),
    n_missing = sum(is.na(values)),
    proportion_missing = mean(is.na(values)),
    stringsAsFactors = FALSE
  )
}
missingness <- do.call(rbind, missingness_rows)
aquatic_missingness <- do.call(rbind, aquatic_missingness_rows)
write_audit_csv(
  missingness,
  file.path(audit_dir, "environment_missingness_44_sites.csv")
)
write_audit_csv(
  aquatic_missingness,
  file.path(audit_dir, "environment_missingness_aquatic_59_sites.csv")
)

# 6. Verify site coverage and record package versions.
stopifnot(
  nrow(environment_aquatic) == 59L,
  nrow(environment_60) == 60L,
  nrow(environment_44) == 44L,
  !anyDuplicated(environment_aquatic$site_code),
  !anyDuplicated(environment_44$site_code),
  !anyNA(environment_44$sediment_sample_id),
  setequal(
    environment_aquatic$site_code,
    unique(aquatic_manifest$site_code)
  )
)
capture_session(file.path(audit_dir, "session_info_environment.txt"))
message(
  "Environment preparation complete: 59 sequenced aquatic sites, ",
  "60 FT-ICR sites, and 44 paired sediment sites."
)
