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
assert_columns(master, c("Site Code", "Sample Type"), "master environmental CSV")
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
stopifnot(
  nrow(environment_60) == 60L,
  nrow(environment_44) == 44L,
  !anyDuplicated(environment_44$site_code),
  !anyNA(environment_44$sediment_sample_id)
)
capture_session(file.path(audit_dir, "session_info_environment.txt"))
message("Environment preparation complete: 60 FT-ICR sites and 44 paired sediment sites.")
