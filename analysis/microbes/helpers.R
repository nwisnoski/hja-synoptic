# Small utilities shared by the microbial-diversity scripts.

parse_microbe_args <- function(args, script_dir) {
  values <- list(
    project_root = normalizePath(file.path(script_dir, "..", "..")),
    output_root = NA_character_
  )
  flags <- c("--project-root" = "project_root", "--output-root" = "output_root")
  i <- 1L
  while (i <= length(args)) {
    if (!args[[i]] %in% names(flags) || i == length(args)) {
      stop("Unknown or incomplete argument: ", args[[i]], call. = FALSE)
    }
    values[[flags[[args[[i]]]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values
}

under_project_root <- function(path, project_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(project_root, path)
}

read_analysis_csv <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

write_analysis_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
}

assert_files <- function(paths, label = "Required input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(label, " missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

assert_unique <- function(x, label) {
  if (anyNA(x) || any(x == "") || anyDuplicated(x)) {
    stop(label, " must contain complete, unique IDs.", call. = FALSE)
  }
}

capture_analysis_session <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}

load_diversity_input <- function(output_root) {
  path <- file.path(output_root, "derived", "diversity_input.rds")
  if (!file.exists(path)) {
    stop("Run 01_prepare_diversity.R first: ", path, call. = FALSE)
  }
  readRDS(path)
}
