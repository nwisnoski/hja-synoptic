# Shared base-R helpers for the HJA analysis-input preparation scripts.

read_source_csv <- function(path) {
  read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA", "NaN"), strip.white = FALSE
  )
}
write_audit_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
}
write_matrix_csv_gz <- function(matrix, id_name, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  table <- data.frame(
    setNames(list(rownames(matrix)), id_name), matrix,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(table, connection, row.names = FALSE, na = "")
}
assert_columns <- function(data, required, source_label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      source_label, " is missing required column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
}
assert_unique_key <- function(key, source_label) {
  duplicates <- unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)])
  if (length(duplicates)) {
    stop(
      source_label, " has duplicate key(s): ",
      paste(head(duplicates, 20L), collapse = ", "), call. = FALSE
    )
  }
}
# Keep relative source names in the audit, while resolving files through here.
file_inventory <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing source file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  source_paths <- substring(paths, nchar(here::here()) + 2L)
  info <- file.info(paths)
  data.frame(
    source_path = unname(source_paths),
    absolute_path = normalizePath(paths),
    bytes = unname(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}
capture_session <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}
