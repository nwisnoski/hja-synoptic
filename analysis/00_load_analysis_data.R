# Load and validate all HJA analysis inputs from their documented derived files.
#
# Interactive use from anywhere inside the repository:
#   source("analysis/00_load_analysis_data.R")
#   hja <- load_hja_analysis_data()
#   hja
#
# This script never writes files. Named transformed views are derived in memory
# from raw counts/intensities and are recorded in hja$transformation_log.

find_hja_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "hja-synoptic.Rproj"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not find hja-synoptic.Rproj above: ", start,
        ". Supply project_root explicitly.",
        call. = FALSE
      )
    }
    current <- parent
  }
}

read_analysis_csv <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

read_matrix_csv <- function(path, id_column = 1L) {
  table <- read_analysis_csv(path)
  ids <- table[[id_column]]
  values <- as.matrix(table[-id_column])
  storage.mode(values) <- "numeric"
  rownames(values) <- ids
  values
}

read_dada2_asv_rds <- function(sequence_table_path, sequence_map_path) {
  values <- readRDS(sequence_table_path)
  sequence_map <- read_analysis_csv(sequence_map_path)
  required <- c("asv_id", "sequence")
  if (!all(required %in% names(sequence_map))) {
    stop("DADA2 sequence map is missing asv_id or sequence.", call. = FALSE)
  }
  asv_rows <- match(colnames(values), sequence_map$sequence)
  if (anyNA(asv_rows)) {
    stop("Some DADA2 sequence-table columns lack ASV identifiers.", call. = FALSE)
  }
  if (anyDuplicated(sequence_map$asv_id[asv_rows])) {
    stop("DADA2 ASV identifiers are not unique.", call. = FALSE)
  }
  colnames(values) <- sequence_map$asv_id[asv_rows]
  values
}

assert_file_set <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      label, " is incomplete. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

assert_unique_ids <- function(ids, label) {
  if (anyNA(ids) || any(ids == "")) {
    stop(label, " contains missing identifiers.", call. = FALSE)
  }
  duplicates <- unique(ids[duplicated(ids)])
  if (length(duplicates)) {
    stop(
      label, " contains duplicate identifiers: ",
      paste(head(duplicates, 20L), collapse = ", "),
      call. = FALSE
    )
  }
}

row_relative <- function(x) {
  totals <- rowSums(x, na.rm = TRUE)
  denominators <- totals
  denominators[denominators <= 0] <- NA_real_
  sweep(x, 1L, denominators, "/")
}

hellinger_transform <- function(x) sqrt(row_relative(x))

presence_absence <- function(x, threshold = 0) {
  result <- 1L * (x > threshold)
  dimnames(result) <- dimnames(x)
  result
}

make_compositional_views <- function(x, data_type, threshold = 0) {
  totals <- rowSums(x, na.rm = TRUE)
  list(
    raw = x,
    row_relative = row_relative(x),
    hellinger = hellinger_transform(x),
    presence_absence = presence_absence(x, threshold = threshold),
    row_totals = totals,
    zero_total_rows = rownames(x)[totals <= 0],
    data_type = data_type
  )
}

load_hja_analysis_data <- function(
    project_root = find_hja_project_root(),
    require_dada2 = FALSE,
    prepare_views = TRUE,
    verify_source_checksums = TRUE) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  input_root <- file.path(project_root, "data", "derived", "analysis_inputs")
  audit_root <- file.path(input_root, "audit")
  environment_root <- file.path(input_root, "environment")
  fticr_root <- file.path(input_root, "fticr")
  dada2_root <- file.path(project_root, "results", "dada2_2016")

  base_paths <- c(
    status = file.path(input_root, "preparation_status.csv"),
    inventory = file.path(audit_root, "source_file_inventory.csv"),
    crosswalk = file.path(audit_root, "sediment_site_crosswalk.csv"),
    variable_dictionary = file.path(audit_root, "environment_variable_dictionary.csv"),
    missingness = file.path(audit_root, "environment_missingness_44_sites.csv"),
    aquatic_missingness = file.path(
      audit_root, "environment_missingness_aquatic_59_sites.csv"
    ),
    fticr_filter_log = file.path(audit_root, "fticr_filter_log.csv"),
    environment_aquatic = file.path(
      environment_root, "aquatic_59_site_environment.csv"
    ),
    environment_60 = file.path(environment_root, "fticr_60_site_environment.csv"),
    environment_44 = file.path(environment_root, "sediment_44_environment.csv"),
    fticr_raw = file.path(fticr_root, "fticr_site_by_peak.rds"),
    fticr_primary = file.path(fticr_root, "fticr_site_by_peak_primary.rds"),
    fticr_qc = file.path(fticr_root, "fticr_qc_by_peak.rds"),
    peak_metadata = file.path(fticr_root, "fticr_peak_metadata.csv.gz"),
    peak_filter_status = file.path(fticr_root, "fticr_peak_filter_status.csv.gz")
  )
  assert_file_set(base_paths, "Base environmental/FT-ICR analysis input set")

  status <- read_analysis_csv(base_paths[["status"]])
  inventory <- read_analysis_csv(base_paths[["inventory"]])
  if (verify_source_checksums) {
    assert_file_set(inventory$absolute_path, "Recorded raw source set")
    current_md5 <- unname(tools::md5sum(inventory$absolute_path))
    if (!identical(current_md5, inventory$md5)) {
      changed <- inventory$source_path[current_md5 != inventory$md5]
      stop(
        "Raw source checksum mismatch: ", paste(changed, collapse = ", "),
        call. = FALSE
      )
    }
  }

  crosswalk <- read_analysis_csv(base_paths[["crosswalk"]])
  variable_dictionary <- read_analysis_csv(base_paths[["variable_dictionary"]])
  environmental_missingness <- read_analysis_csv(base_paths[["missingness"]])
  aquatic_environmental_missingness <- read_analysis_csv(
    base_paths[["aquatic_missingness"]]
  )
  environment_aquatic <- read_analysis_csv(base_paths[["environment_aquatic"]])
  environment_60 <- read_analysis_csv(base_paths[["environment_60"]])
  environment_44 <- read_analysis_csv(base_paths[["environment_44"]])
  fticr_raw <- readRDS(base_paths[["fticr_raw"]])
  fticr_primary <- readRDS(base_paths[["fticr_primary"]])
  fticr_qc <- readRDS(base_paths[["fticr_qc"]])
  peak_metadata <- read_analysis_csv(base_paths[["peak_metadata"]])
  peak_filter_status <- read_analysis_csv(base_paths[["peak_filter_status"]])

  assert_unique_ids(crosswalk$site_code, "FT-ICR site crosswalk")
  assert_unique_ids(
    environment_aquatic$site_code, "59-site aquatic environmental table"
  )
  assert_unique_ids(environment_60$site_code, "60-site environmental table")
  assert_unique_ids(environment_44$site_code, "44-site environmental table")
  assert_unique_ids(peak_metadata$peak_id, "FT-ICR peak metadata")
  if (!identical(environment_60$site_code, rownames(fticr_primary))) {
    stop("The 60-site environment and FT-ICR rows are not aligned.", call. = FALSE)
  }
  primary_peak_ids <- peak_filter_status$peak_id[peak_filter_status$included_primary]
  if (!identical(primary_peak_ids, colnames(fticr_primary))) {
    stop("Primary FT-ICR feature IDs do not match the filter audit.", call. = FALSE)
  }
  if (!identical(peak_metadata$peak_id, colnames(fticr_raw))) {
    stop("Raw FT-ICR feature IDs do not match peak metadata.", call. = FALSE)
  }

  primary_variables <- variable_dictionary$analysis_column[
    variable_dictionary$primary_44_site
  ]
  if (!all(primary_variables %in% names(environment_44))) {
    stop("Declared primary environmental variables are missing.", call. = FALSE)
  }
  sensitivity_variables <- setdiff(
    variable_dictionary$analysis_column, primary_variables
  )
  complete_primary_variables <- environmental_missingness$analysis_column[
    environmental_missingness$primary_44_site &
      environmental_missingness$n_missing == 0
  ]
  incomplete_primary_variables <- setdiff(
    primary_variables, complete_primary_variables
  )

  transformation_log <- data.frame(
    object = c(
      "fticr_60_primary$row_relative",
      "fticr_60_primary$hellinger",
      "fticr_60_primary$presence_absence"
    ),
    input = "fticr_site_by_peak_primary.rds",
    operation = c(
      "x / rowSums(x)",
      "sqrt(x / rowSums(x))",
      "1 * (x > 0)"
    ),
    parameters = c("", "", "threshold=0"),
    stringsAsFactors = FALSE
  )

  if (!prepare_views) transformation_log <- transformation_log[0, , drop = FALSE]
  result <- list(
    project_root = project_root,
    preparation_status = status,
    source_checksums_verified = verify_source_checksums,
    audit = list(
      source_inventory = inventory,
      site_crosswalk = crosswalk,
      environmental_variables = variable_dictionary,
      environmental_missingness_aquatic_59 =
        aquatic_environmental_missingness,
      environmental_missingness_44 = environmental_missingness,
      fticr_filter_log = read_analysis_csv(base_paths[["fticr_filter_log"]])
    ),
    environment = list(
      aquatic_59_sites = environment_aquatic,
      fticr_60_sites = environment_60,
      sediment_44_sites = environment_44,
      sediment_44_declared_primary_predictors = environment_44[
        c("site_code", primary_variables)
      ],
      sediment_44_complete_primary_predictors = environment_44[
        c("site_code", complete_primary_variables)
      ],
      sediment_44_incomplete_primary_predictors = environment_44[
        c("site_code", incomplete_primary_variables)
      ],
      sediment_44_sensitivity_predictors = environment_44[
        c("site_code", sensitivity_variables)
      ],
      predictor_blocks = split(
        variable_dictionary$analysis_column,
        variable_dictionary$predictor_block
      )
    ),
    fticr = list(
      intensity_60_raw_features = fticr_raw,
      intensity_60_primary_features = fticr_primary,
      laboratory_standards = fticr_qc,
      peak_metadata = peak_metadata,
      peak_filter_status = peak_filter_status
    ),
    microbes = NULL,
    sediment_multiblock = NULL,
    source_pools = NULL,
    dada2_ready = FALSE,
    transformation_log = transformation_log
  )
  if (prepare_views) {
    result$views <- list(
      fticr_60_primary = make_compositional_views(
        fticr_primary, "FT-ICR-MS peak intensity"
      )
    )
  } else {
    result$views <- NULL
  }

  dada2_marker <- file.path(dada2_root, "session_info.txt")
  sediment_root <- file.path(input_root, "sediment_multiblock")
  source_root <- file.path(input_root, "source_pools")
  microbial_paths <- c(
    asv_all = file.path(dada2_root, "sequence_table_asv.rds"),
    asv_sequences = file.path(dada2_root, "asv_sequences.csv"),
    sample_metadata = file.path(dada2_root, "sample_metadata_and_read_tracking.csv"),
    sediment_metadata = file.path(sediment_root, "sediment_44_sample_metadata.csv"),
    sediment_environment = file.path(sediment_root, "sediment_44_environment.csv"),
    sediment_asv = file.path(sediment_root, "sediment_44_asv_counts.rds"),
    sediment_fticr = file.path(sediment_root, "sediment_44_fticr_intensity.rds"),
    sediment_fticr_pa = file.path(sediment_root, "sediment_44_fticr_presence_absence.rds"),
    source_map = file.path(source_root, "sediment_44_same_site_source_map.csv"),
    aquatic_metadata = file.path(source_root, "same_site_aquatic_sample_metadata.csv"),
    aquatic_asv = file.path(source_root, "same_site_aquatic_asv_counts.rds"),
    soil_metadata = file.path(source_root, "regional_soil_sample_metadata.csv"),
    soil_asv = file.path(source_root, "regional_soil_asv_counts.rds")
  )

  if (!file.exists(dada2_marker)) {
    if (require_dada2) {
      stop("DADA2 is not complete; completion marker is absent.", call. = FALSE)
    }
    class(result) <- "hja_analysis_data"
    return(result)
  }
  assert_file_set(
    microbial_paths,
    paste0(
      "DADA2 is complete, but integrated inputs are absent. Rerun ",
      "analysis/data_prep/run_data_preparation.R; the microbial input set"
    )
  )

  asv_all <- read_dada2_asv_rds(
    microbial_paths[["asv_all"]],
    microbial_paths[["asv_sequences"]]
  )
  sample_metadata <- read_analysis_csv(microbial_paths[["sample_metadata"]])
  assert_unique_ids(rownames(asv_all), "DADA2 ASV count-table rows")
  assert_unique_ids(sample_metadata$sample_id, "DADA2 sample metadata")
  sample_metadata_rows <- match(rownames(asv_all), sample_metadata$sample_id)
  if (anyNA(sample_metadata_rows)) {
    stop("Some ASV-table samples are absent from DADA2 sample metadata.", call. = FALSE)
  }
  sample_metadata <- sample_metadata[sample_metadata_rows, , drop = FALSE]
  sediment_asv <- readRDS(microbial_paths[["sediment_asv"]])
  sediment_fticr <- readRDS(microbial_paths[["sediment_fticr"]])
  sediment_fticr_pa <- readRDS(microbial_paths[["sediment_fticr_pa"]])
  sediment_environment <- read_analysis_csv(
    microbial_paths[["sediment_environment"]]
  )
  expected_sediment_sites <- sediment_environment$site_code
  if (!identical(rownames(sediment_asv), expected_sediment_sites) ||
      !identical(rownames(sediment_fticr), expected_sediment_sites) ||
      !identical(rownames(sediment_fticr_pa), expected_sediment_sites)) {
    stop("The 44-site sediment data blocks are not aligned.", call. = FALSE)
  }

  aquatic_asv <- readRDS(microbial_paths[["aquatic_asv"]])
  soil_asv <- readRDS(microbial_paths[["soil_asv"]])
  aquatic_metadata <- read_analysis_csv(microbial_paths[["aquatic_metadata"]])
  soil_metadata <- read_analysis_csv(microbial_paths[["soil_metadata"]])
  aquatic_metadata <- aquatic_metadata[
    match(rownames(aquatic_asv), aquatic_metadata$sample_id), , drop = FALSE
  ]
  soil_metadata <- soil_metadata[
    match(rownames(soil_asv), soil_metadata$sample_id), , drop = FALSE
  ]
  if (anyNA(aquatic_metadata$sample_id) || anyNA(soil_metadata$sample_id)) {
    stop("Source-pool ASV rows are absent from their metadata.", call. = FALSE)
  }
  result$microbes <- list(
    asv_counts_all = asv_all,
    sample_metadata = sample_metadata,
    asv_sequences = if (file.exists(file.path(dada2_root, "asv_sequences.csv"))) {
      read_analysis_csv(file.path(dada2_root, "asv_sequences.csv"))
    } else NULL,
    taxonomy = if (file.exists(file.path(dada2_root, "asv_taxonomy.csv"))) {
      read_analysis_csv(file.path(dada2_root, "asv_taxonomy.csv"))
    } else NULL
  )
  result$sediment_multiblock <- list(
    sample_metadata = read_analysis_csv(microbial_paths[["sediment_metadata"]]),
    environment = sediment_environment,
    asv_counts = sediment_asv,
    fticr_intensity = sediment_fticr,
    fticr_presence_absence = sediment_fticr_pa
  )
  result$source_pools <- list(
    paired_site_map = read_analysis_csv(microbial_paths[["source_map"]]),
    aquatic_sample_metadata = aquatic_metadata,
    aquatic_asv_counts = aquatic_asv,
    regional_soil_metadata = soil_metadata,
    regional_soil_asv_counts = soil_asv
  )
  result$dada2_ready <- TRUE

  if (prepare_views) {
    result$views$microbes_all <- make_compositional_views(
      asv_all, "16S V4 ASV count"
    )
    result$views$sediment_asv <- make_compositional_views(
      sediment_asv, "16S V4 ASV count"
    )
    result$views$sediment_fticr <- make_compositional_views(
      sediment_fticr, "FT-ICR-MS peak intensity"
    )
    result$views$source_aquatic_asv <- make_compositional_views(
      aquatic_asv, "16S V4 ASV count"
    )
    result$views$source_soil_asv <- make_compositional_views(
      soil_asv, "16S V4 ASV count"
    )
    new_log <- do.call(rbind, lapply(
      c(
        "microbes_all", "sediment_asv", "sediment_fticr",
        "source_aquatic_asv", "source_soil_asv"
      ),
      function(object) {
        data.frame(
          object = paste0(
            object,
            c("$row_relative", "$hellinger", "$presence_absence")
          ),
          input = paste0(object, "$raw"),
          operation = c(
            "x / rowSums(x)",
            "sqrt(x / rowSums(x))",
            "1 * (x > 0)"
          ),
          parameters = c("", "", "threshold=0"),
          stringsAsFactors = FALSE
        )
      }
    ))
    result$transformation_log <- rbind(result$transformation_log, new_log)
  }
  class(result) <- "hja_analysis_data"
  result
}

print.hja_analysis_data <- function(x, ...) {
  cat("HJA analysis data\n")
  cat(
    "  Sequenced aquatic sites: ",
    nrow(x$environment$aquatic_59_sites), "\n", sep = ""
  )
  cat("  FT-ICR sites: ", nrow(x$environment$fticr_60_sites), "\n", sep = "")
  cat(
    "  Primary FT-ICR features: ",
    ncol(x$fticr$intensity_60_primary_features), "\n", sep = ""
  )
  cat(
    "  Paired sediment sites declared: ",
    nrow(x$environment$sediment_44_sites), "\n", sep = ""
  )
  cat("  DADA2 integrated data ready: ", x$dada2_ready, "\n", sep = "")
  invisible(x)
}
