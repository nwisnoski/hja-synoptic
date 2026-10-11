#!/usr/bin/env Rscript

# Create representative input quality profiles before the long DADA2 run.
# Samples are spaced across each run so the review is not biased toward the
# first FASTQ files.

parse_args <- function(args) {
  values <- list(
    project_root = ".",
    manifest = "data/derived/hja_2016_sample_manifest.csv",
    output = "results/dada2_2016",
    samples_per_run = 12L
  )
  flags <- c(
    "--project-root" = "project_root",
    "--manifest" = "manifest",
    "--output" = "output",
    "--samples-per-run" = "samples_per_run"
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
  values$samples_per_run <- as.integer(values$samples_per_run)
  if (is.na(values$samples_per_run) || values$samples_per_run < 1L) {
    stop("--samples-per-run must be a positive integer.", call. = FALSE)
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
if (!requireNamespace("dada2", quietly = TRUE)) {
  stop("The Bioconductor package dada2 is required.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The CRAN package ggplot2 is required.", call. = FALSE)
}

root <- normalizePath(args$project_root, mustWork = TRUE)
manifest_path <- under_root(args$manifest, root)
output_dir <- under_root(args$output, root)
if (!file.exists(manifest_path)) {
  stop("Missing manifest: ", manifest_path, call. = FALSE)
}

manifest <- read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
required <- c("sample_id", "run_id", "forward_path", "reverse_path", "include")
if (!all(required %in% names(manifest))) {
  stop(
    "Manifest is missing: ",
    paste(setdiff(required, names(manifest)), collapse = ", "),
    call. = FALSE
  )
}
manifest <- manifest[as_flag(manifest$include), , drop = FALSE]
manifest$forward_file <- vapply(
  manifest$forward_path,
  under_root,
  root = root,
  FUN.VALUE = character(1)
)
manifest$reverse_file <- vapply(
  manifest$reverse_path,
  under_root,
  root = root,
  FUN.VALUE = character(1)
)
missing_fastqs <- unique(c(
  manifest$forward_file[!file.exists(manifest$forward_file)],
  manifest$reverse_file[!file.exists(manifest$reverse_file)]
))
if (length(missing_fastqs)) {
  stop("Missing FASTQ file(s): ", paste(missing_fastqs, collapse = ", "), call. = FALSE)
}

selected_rows <- list()
for (run_id in unique(manifest$run_id)) {
  run_samples <- manifest[manifest$run_id == run_id, , drop = FALSE]
  run_samples <- run_samples[order(run_samples$sample_id), , drop = FALSE]
  n_select <- min(args$samples_per_run, nrow(run_samples))
  selected_index <- unique(as.integer(round(seq(
    1,
    nrow(run_samples),
    length.out = n_select
  ))))
  selected <- run_samples[selected_index, , drop = FALSE]
  run_dir <- file.path(output_dir, "runs", run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(run_dir, "quality_profiles_input.pdf")

  grDevices::pdf(pdf_path, width = 9, height = 6, onefile = TRUE)
  for (i in seq_len(nrow(selected))) {
    print(
      dada2::plotQualityProfile(selected$forward_file[[i]]) +
        ggplot2::ggtitle(paste(selected$sample_id[[i]], "forward"))
    )
    print(
      dada2::plotQualityProfile(selected$reverse_file[[i]]) +
        ggplot2::ggtitle(paste(selected$sample_id[[i]], "reverse"))
    )
  }
  grDevices::dev.off()

  selected_rows[[run_id]] <- data.frame(
    run_id = run_id,
    sample_id = selected$sample_id,
    forward_path = selected$forward_path,
    reverse_path = selected$reverse_path,
    quality_pdf = file.path(
      args$output,
      "runs",
      run_id,
      "quality_profiles_input.pdf"
    ),
    stringsAsFactors = FALSE
  )
  message("Wrote ", normalizePath(pdf_path))
}

quality_index <- do.call(rbind, selected_rows)
rownames(quality_index) <- NULL
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  quality_index,
  file.path(output_dir, "quality_profile_samples.csv"),
  row.names = FALSE,
  na = ""
)
message(
  "\nQuality review is ready. Inspect every PDF before running dada2_pipeline.R."
)
