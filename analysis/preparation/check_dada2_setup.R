#!/usr/bin/env Rscript

# Fast preflight for the 2016 HJA DADA2 workflow. This script does not read
# sequence contents or alter FASTQs.

parse_args <- function(args) {
  values <- list(
    project_root = ".",
    manifest = "data/derived/hja_2016_sample_manifest.csv",
    run_config = "config/dada2_run_config.csv",
    output = "results/dada2_2016"
  )
  flags <- c(
    "--project-root" = "project_root",
    "--manifest" = "manifest",
    "--run-config" = "run_config",
    "--output" = "output"
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

is_absolute <- function(path) grepl("^(/|[A-Za-z]:[/\\\\])", path)
under_root <- function(path, root) {
  if (is_absolute(path)) path else file.path(root, path)
}
as_flag <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "YES", "Y", "1")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(args$project_root, mustWork = TRUE)
manifest_path <- under_root(args$manifest, root)
config_path <- under_root(args$run_config, root)
output_dir <- under_root(args$output, root)

required_inputs <- c(manifest_path, config_path)
if (any(!file.exists(required_inputs))) {
  stop(
    "Missing required input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", "),
    call. = FALSE
  )
}

manifest <- read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
config <- read.csv(
  config_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
required_manifest <- c(
  "sample_id", "run_id", "forward_path", "reverse_path", "include"
)
required_config <- c(
  "run_id", "trim_left_f", "trim_left_r", "trunc_len_f", "trunc_len_r",
  "max_ee_f", "max_ee_r", "trunc_q", "min_len", "max_n", "rm_phix",
  "min_overlap", "max_mismatch", "pool", "threads"
)
missing_manifest_columns <- setdiff(required_manifest, names(manifest))
missing_config_columns <- setdiff(required_config, names(config))
if (length(missing_manifest_columns) || length(missing_config_columns)) {
  stop(
    "Missing columns. Manifest: ",
    paste(missing_manifest_columns, collapse = ", "),
    "; config: ",
    paste(missing_config_columns, collapse = ", "),
    call. = FALSE
  )
}

manifest <- manifest[as_flag(manifest$include), , drop = FALSE]
forward <- vapply(
  manifest$forward_path,
  under_root,
  root = root,
  FUN.VALUE = character(1)
)
reverse <- vapply(
  manifest$reverse_path,
  under_root,
  root = root,
  FUN.VALUE = character(1)
)
all_fastqs <- c(forward, reverse)
fastq_info <- file.info(all_fastqs)
missing_fastqs <- all_fastqs[is.na(fastq_info$size)]
configured_runs <- unique(config$run_id)
manifest_runs <- unique(manifest$run_id)

package_status <- c(
  dada2 = requireNamespace("dada2", quietly = TRUE),
  ggplot2 = requireNamespace("ggplot2", quietly = TRUE)
)
package_versions <- vapply(
  names(package_status),
  function(package) {
    if (package_status[[package]]) {
      as.character(utils::packageVersion(package))
    } else {
      "NOT INSTALLED"
    }
  },
  character(1)
)

disk_free_gib <- NA_real_
if (.Platform$OS.type == "unix") {
  disk_lines <- tryCatch(
    system2("df", c("-Pk", root), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(disk_lines) >= 2L) {
    fields <- strsplit(trimws(tail(disk_lines, 1L)), "[[:space:]]+")[[1]]
    if (length(fields) >= 4L) disk_free_gib <- as.numeric(fields[[4]]) / 1024^2
  }
}

checks <- data.frame(
  check = c(
    paste0("package_", names(package_status)),
    "unique_sample_ids",
    "all_fastq_files_present",
    "two_fastqs_per_library",
    "all_manifest_runs_configured",
    "config_run_ids_unique",
    "expected_2016_library_count",
    "available_disk_gib",
    "raw_fastq_gib",
    "ready_for_quality_review"
  ),
  value = c(
    package_versions,
    as.character(!anyDuplicated(manifest$sample_id)),
    as.character(!length(missing_fastqs)),
    as.character(length(all_fastqs) == 2L * nrow(manifest)),
    as.character(all(manifest_runs %in% configured_runs)),
    as.character(!anyDuplicated(config$run_id)),
    as.character(nrow(manifest) == 120L),
    ifelse(is.na(disk_free_gib), "not measured", sprintf("%.1f", disk_free_gib)),
    sprintf("%.2f", sum(fastq_info$size, na.rm = TRUE) / 1024^3),
    as.character(
      all(package_status) &&
        !anyDuplicated(manifest$sample_id) &&
        !length(missing_fastqs) &&
        length(all_fastqs) == 2L * nrow(manifest) &&
        all(manifest_runs %in% configured_runs) &&
        !anyDuplicated(config$run_id) &&
        nrow(manifest) == 120L
    )
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  checks,
  file.path(output_dir, "setup_check.csv"),
  row.names = FALSE,
  na = ""
)
print(checks, row.names = FALSE)

if (tail(checks$value, 1L) != "TRUE") {
  stop("Setup check failed. Review the FALSE rows above.", call. = FALSE)
}
