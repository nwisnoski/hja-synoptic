# Shared base-R helpers for the HJA analysis-input preparation scripts.

parse_named_args <- function(args, defaults, flags) {
  values <- defaults
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

is_absolute_path <- function(path) grepl("^(/|[A-Za-z]:[/\\\\])", path)
under_project_root <- function(path, project_root) {
  if (is_absolute_path(path)) path else file.path(project_root, path)
}
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
one_value_or_na <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (!length(x)) return(NA_character_)
  if (length(x) > 1L) {
    stop("Expected at most one value, found: ", paste(x, collapse = ", "), call. = FALSE)
  }
  x
}
file_inventory <- function(paths, project_root) {
  resolved <- vapply(
    paths, under_project_root, project_root = project_root,
    FUN.VALUE = character(1)
  )
  missing <- resolved[!file.exists(resolved)]
  if (length(missing)) {
    stop("Missing source file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  info <- file.info(resolved)
  data.frame(
    source_path = unname(paths),
    absolute_path = normalizePath(resolved),
    bytes = unname(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    md5 = unname(tools::md5sum(resolved)),
    stringsAsFactors = FALSE
  )
}
capture_session <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}
