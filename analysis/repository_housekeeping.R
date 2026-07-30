#!/usr/bin/env Rscript

# Reversibly move explicitly approved, unused files into archive/.
# Default behavior is a dry run:
#   Rscript analysis/repository_housekeeping.R
# Apply the plan:
#   Rscript analysis/repository_housekeeping.R --apply
# Restore every archived file:
#   Rscript analysis/repository_housekeeping.R --restore

command <- commandArgs()
script_file <- sub("^--file=", "", grep("^--file=", command, value = TRUE)[1])
args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), ".."))
project_index <- match("--project-root", args)
if (!is.na(project_index)) {
  if (project_index == length(args)) {
    stop("--project-root requires a path.", call. = FALSE)
  }
  project_root <- normalizePath(args[project_index + 1L], mustWork = TRUE)
  args <- args[-c(project_index, project_index + 1L)]
}
if (length(setdiff(args, c("--apply", "--restore")))) {
  stop(
    "Allowed options are --apply, --restore, and --project-root PATH.",
    call. = FALSE
  )
}
if (all(c("--apply", "--restore") %in% args)) {
  stop("Choose either --apply or --restore, not both.", call. = FALSE)
}
mode <- if ("--apply" %in% args) {
  "apply"
} else if ("--restore" %in% args) {
  "restore"
} else {
  "dry_run"
}

plan_path <- file.path(project_root, "config", "repository_tidy_plan.csv")
plan <- read.csv(plan_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "source_path", "archive_path", "category", "reason", "expected_md5"
)
if (!all(required %in% names(plan))) {
  stop("The tidy plan is missing required columns.", call. = FALSE)
}
if (anyDuplicated(plan$source_path) || anyDuplicated(plan$archive_path)) {
  stop("The tidy plan contains duplicate paths.", call. = FALSE)
}

source_absolute <- file.path(project_root, plan$source_path)
archive_absolute <- file.path(project_root, plan$archive_path)
status <- character(nrow(plan))

for (i in seq_len(nrow(plan))) {
  source_exists <- file.exists(source_absolute[i])
  archive_exists <- file.exists(archive_absolute[i])
  if (source_exists && archive_exists) {
    stop(
      "Both source and archive destination exist for: ", plan$source_path[i],
      call. = FALSE
    )
  }
  existing <- if (source_exists) source_absolute[i] else archive_absolute[i]
  if (!file.exists(existing)) {
    stop(
      "Neither source nor archive destination exists for: ",
      plan$source_path[i], call. = FALSE
    )
  }
  observed_md5 <- unname(tools::md5sum(existing))
  if (!identical(observed_md5, plan$expected_md5[i])) {
    stop("Checksum mismatch for: ", existing, call. = FALSE)
  }

  if (mode == "dry_run") {
    status[i] <- if (source_exists) "would_archive" else "already_archived"
  } else if (mode == "apply") {
    if (archive_exists) {
      status[i] <- "already_archived"
    } else {
      dir.create(
        dirname(archive_absolute[i]), recursive = TRUE, showWarnings = FALSE
      )
      if (!file.rename(source_absolute[i], archive_absolute[i])) {
        stop("Could not archive: ", plan$source_path[i], call. = FALSE)
      }
      status[i] <- "archived"
    }
  } else {
    if (source_exists) {
      status[i] <- "already_restored"
    } else {
      dir.create(
        dirname(source_absolute[i]), recursive = TRUE, showWarnings = FALSE
      )
      if (!file.rename(archive_absolute[i], source_absolute[i])) {
        stop("Could not restore: ", plan$source_path[i], call. = FALSE)
      }
      status[i] <- "restored"
    }
  }
}

report <- data.frame(
  plan,
  status = status,
  action_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
print(report[c("source_path", "archive_path", "status")], row.names = FALSE)
if (mode != "dry_run") {
  manifest_path <- file.path(project_root, "archive", "repository_tidy_manifest.csv")
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(report, manifest_path, row.names = FALSE, na = "")
  message("Manifest written to: ", manifest_path)
} else {
  message("Dry run only; no files were moved.")
}
