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
environment_dir <- file.path(output_root, "environment")
dir.create(environment_dir, recursive = TRUE, showWarnings = FALSE)

crosswalk_path <- file.path(audit_dir, "sediment_site_crosswalk.csv")
if (!file.exists(crosswalk_path)) {
  stop("Run 01_inventory_and_crosswalk.R first.", call. = FALSE)
}
crosswalk <- read_source_csv(crosswalk_path)
master <- read_source_csv(under_project_root(
  prep_config$paths$master_environment, project_root
))
manifest <- read_source_csv(under_project_root(
  prep_config$paths$sample_manifest, project_root
))
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

build_environment_table <- function(site_codes) {
  result <- data.frame(site_code = site_codes, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(environment_variable_map))) {
    rule <- environment_variable_map[i, , drop = FALSE]
    key <- paste(master[["Site Code"]], master[["Sample Type"]], sep = "::")
    target <- paste(site_codes, rule$source_sample_type, sep = "::")
    matched <- match(target, key)
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
sample_id_for <- function(site_code, sample_type) {
  matches <- aquatic_manifest$sample_id[
    aquatic_manifest$site_code == site_code &
      aquatic_manifest$sample_type == sample_type
  ]
  if (!length(matches)) return(NA_character_)
  if (length(matches) > 1L) {
    stop("Duplicate aquatic samples for ", site_code, " / ", sample_type)
  }
  matches
}
aquatic_index <- data.frame(
  site_code = aquatic_site_codes,
  planktonic_sample_id = vapply(
    aquatic_site_codes, sample_id_for,
    sample_type = "planktonic streamwater", FUN.VALUE = character(1)
  ),
  hyporheic_sample_id = vapply(
    aquatic_site_codes, sample_id_for,
    sample_type = "hyporheic water", FUN.VALUE = character(1)
  ),
  sediment_sample_id = vapply(
    aquatic_site_codes, sample_id_for,
    sample_type = "stream sediment", FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)
aquatic_index$fticr_source_column <- crosswalk$fticr_source_column[
  match(aquatic_index$site_code, crosswalk$site_code)
]
aquatic_index$has_fticr_profile <- !is.na(aquatic_index$fticr_source_column)
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

missingness <- do.call(rbind, lapply(seq_len(nrow(environment_variable_map)), function(i) {
  rule <- environment_variable_map[i, , drop = FALSE]
  values <- environment_44[[rule$analysis_column]]
  data.frame(
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
}))
write_audit_csv(
  missingness,
  file.path(audit_dir, "environment_missingness_44_sites.csv")
)
aquatic_missingness <- do.call(
  rbind,
  lapply(seq_len(nrow(environment_variable_map)), function(i) {
    rule <- environment_variable_map[i, , drop = FALSE]
    values <- environment_aquatic[[rule$analysis_column]]
    data.frame(
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
  })
)
write_audit_csv(
  aquatic_missingness,
  file.path(audit_dir, "environment_missingness_aquatic_59_sites.csv")
)
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
