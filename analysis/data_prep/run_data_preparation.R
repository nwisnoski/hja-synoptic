#!/usr/bin/env Rscript

# Prepare the 2016 environmental, molecular, and sample-identity tables.
# Each step can also be sourced on its own in Positron.

library(here)

# 1. Settings.
source(here("analysis", "data_prep", "config.R"), local = TRUE)
output_root <- prep_config$paths$output_root
dada2_marker <- file.path(prep_config$paths$dada2_output, "session_info.txt")
steps <- c(
  "01_inventory_and_crosswalk.R",
  "02_prepare_environment.R",
  "03_prepare_fticr.R",
  "04_build_sediment_multiblock.R",
  "05_build_source_pool_tables.R"
)
status <- data.frame(
  stage = c(
    "inventory_and_crosswalk", "environment", "fticr",
    "sediment_multiblock", "source_pool_tables"
  ),
  status = "pending",
  prerequisite = c(rep("", 3L), rep("results/dada2_2016/session_info.txt", 2L)),
  stringsAsFactors = FALSE
)

# 2. Run each preparation step in order, stopping if a step fails.
# Separate environments keep intermediate tables from masking later objects.
for (i in seq_along(steps)) {
  if (i > 3L && !file.exists(dada2_marker)) {
    status$status[i] <- "waiting_for_dada2"
    next
  }
  message("Preparing ", steps[i])
  source(
    here("analysis", "data_prep", steps[i]),
    local = new.env(parent = globalenv())
  )
  status$status[i] <- "complete"
}

# 3. Record completion and the local package versions.
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
write.csv(
  status,
  file.path(output_root, "preparation_status.csv"),
  row.names = FALSE,
  na = ""
)
capture.output(
  sessionInfo(),
  file = file.path(output_root, "audit", "session_info_runner.txt")
)
if (any(status$status == "waiting_for_dada2")) {
  message("Steps 04 and 05 require completed DADA2 outputs; rerun after completion.")
}
message("Preparation status written to: ", file.path(output_root, "preparation_status.csv"))
