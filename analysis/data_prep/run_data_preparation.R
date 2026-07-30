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
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
run_step <- function(script_name) {
  message("\n--- ", script_name, " ---")
  status <- system2(
    rscript,
    c(
      file.path(script_dir, script_name),
      "--project-root", project_root,
      "--output-root", output_root
    )
  )
  if (!identical(status, 0L)) {
    stop(script_name, " failed with exit status ", status, ".", call. = FALSE)
  }
}

base_steps <- c(
  "01_inventory_and_crosswalk.R",
  "02_prepare_environment.R",
  "03_prepare_fticr.R"
)
for (step in base_steps) run_step(step)

dada2_marker <- file.path(
  under_project_root(prep_config$paths$dada2_output, project_root),
  "session_info.txt"
)
status <- data.frame(
  stage = c(
    "inventory_and_crosswalk", "environment", "fticr",
    "sediment_multiblock", "source_pool_tables"
  ),
  status = c(
    rep("complete", 3L),
    if (file.exists(dada2_marker)) rep("pending", 2L) else rep("waiting_for_dada2", 2L)
  ),
  prerequisite = c(
    rep("", 3L),
    rep("results/dada2_2016/session_info.txt", 2L)
  ),
  stringsAsFactors = FALSE
)
if (file.exists(dada2_marker)) {
  run_step("04_build_sediment_multiblock.R")
  status$status[4] <- "complete"
  run_step("05_build_source_pool_tables.R")
  status$status[5] <- "complete"
} else {
  message(
    "\nDADA2 is still incomplete. Steps 04 and 05 were not run; rerun this ",
    "same command after the completion marker appears."
  )
}
write_audit_csv(status, file.path(output_root, "preparation_status.csv"))
capture_session(file.path(output_root, "audit", "session_info_runner.txt"))
message("\nPreparation status written to: ", file.path(output_root, "preparation_status.csv"))
