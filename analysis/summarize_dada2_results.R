#!/usr/bin/env Rscript

# Build compact QC tables after the DADA2 run. Flags are prompts for review,
# not automatic sample-exclusion rules.

parse_args <- function(args) {
  values <- list(
    project_root = ".",
    input = "results/dada2_2016",
    min_final_retention = 0.20
  )
  flags <- c(
    "--project-root" = "project_root",
    "--input" = "input",
    "--min-final-retention" = "min_final_retention"
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
  values$min_final_retention <- as.numeric(values$min_final_retention)
  if (is.na(values$min_final_retention) ||
      values$min_final_retention < 0 ||
      values$min_final_retention > 1) {
    stop("--min-final-retention must be between 0 and 1.", call. = FALSE)
  }
  values
}

is_absolute <- function(path) grepl("^(/|[A-Za-z]:[/\\\\])", path)
under_root <- function(path, root) {
  if (is_absolute(path)) path else file.path(root, path)
}
write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(args$project_root, mustWork = TRUE)
input_dir <- under_root(args$input, root)
tracking_path <- file.path(input_dir, "read_tracking.csv")
sequences_path <- file.path(input_dir, "asv_sequences.csv")
required <- c(tracking_path, sequences_path)
if (any(!file.exists(required))) {
  stop(
    "DADA2 output is incomplete; missing: ",
    paste(required[!file.exists(required)], collapse = ", "),
    call. = FALSE
  )
}

tracking <- read.csv(
  tracking_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
asvs <- read.csv(
  sequences_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_tracking <- c(
  "sample_id", "run_id", "input", "filtered", "denoised_forward",
  "denoised_reverse", "merged", "nonchim"
)
if (!all(required_tracking %in% names(tracking))) {
  stop(
    "Read tracking is missing: ",
    paste(setdiff(required_tracking, names(tracking)), collapse = ", "),
    call. = FALSE
  )
}
if (!all(c("asv_id", "sequence_length", "total_abundance") %in% names(asvs))) {
  stop("ASV sequence table is missing required columns.", call. = FALSE)
}

tracking$filter_retention <- ifelse(
  tracking$input > 0,
  tracking$filtered / tracking$input,
  NA_real_
)
tracking$merge_retention <- ifelse(
  tracking$filtered > 0,
  tracking$merged / tracking$filtered,
  NA_real_
)
tracking$final_retention <- ifelse(
  tracking$input > 0,
  tracking$nonchim / tracking$input,
  NA_real_
)
tracking$qc_flag <- ifelse(
  tracking$nonchim == 0,
  "NO_NONCHIMERIC_READS",
  ifelse(
    tracking$final_retention < args$min_final_retention,
    "LOW_FINAL_RETENTION",
    "OK"
  )
)

summarize_run <- function(run_rows) {
  data.frame(
    run_id = unique(run_rows$run_id),
    samples = nrow(run_rows),
    input_reads = sum(run_rows$input),
    filtered_reads = sum(run_rows$filtered),
    merged_reads = sum(run_rows$merged),
    nonchim_reads = sum(run_rows$nonchim),
    median_filter_retention = stats::median(
      run_rows$filter_retention,
      na.rm = TRUE
    ),
    median_merge_retention = stats::median(
      run_rows$merge_retention,
      na.rm = TRUE
    ),
    median_final_retention = stats::median(
      run_rows$final_retention,
      na.rm = TRUE
    ),
    flagged_samples = sum(run_rows$qc_flag != "OK"),
    stringsAsFactors = FALSE
  )
}
run_summary <- do.call(
  rbind,
  lapply(split(tracking, tracking$run_id), summarize_run)
)
rownames(run_summary) <- NULL

length_counts <- as.data.frame(table(asvs$sequence_length))
names(length_counts) <- c("sequence_length", "asv_count")
length_counts$sequence_length <- as.integer(as.character(
  length_counts$sequence_length
))
length_abundance <- aggregate(
  asvs$total_abundance,
  by = list(sequence_length = asvs$sequence_length),
  FUN = sum
)
names(length_abundance)[[2]] <- "total_abundance"
length_summary <- merge(
  length_counts,
  length_abundance,
  by = "sequence_length",
  all = TRUE,
  sort = TRUE
)

qc_dir <- file.path(input_dir, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(run_summary, file.path(qc_dir, "run_qc_summary.csv"))
write_csv(
  tracking[order(tracking$qc_flag, tracking$sample_id), ],
  file.path(qc_dir, "sample_qc_flags.csv")
)
write_csv(length_summary, file.path(qc_dir, "asv_length_summary.csv"))

cat("\nRun-level read retention\n")
print(run_summary, row.names = FALSE)
cat(
  "\nASVs:", nrow(asvs),
  "\nFlagged samples:", sum(tracking$qc_flag != "OK"),
  "\nQC tables:", normalizePath(qc_dir), "\n"
)
