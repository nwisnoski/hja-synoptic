# Putative within-sample mass-difference transformations, following Danczak et al.
# (2020), https://doi.org/10.1038/s41467-020-19989-y.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Settings and plotting theme ------------------------------------------------
input_dir <- here("data", "derived", "analysis_inputs")
reference_dir <- here("reference", "danczak_2020")
results_dir <- here("results", "fticr_2016", "transformations")
figures_dir <- here("figures", "supplement", "transformations")
absolute_mass_tolerance <- 0.000010
display_transformations <- 20L
upstream_commit <- "434e28bd471c83e1eb616f4203d6ada584520d2e"
reference_md5 <- c("26c1733353bfd0e076b5fe048ca1d486",
  "9955caaf9f539785608552336cb9a5a5")
figure_width <- 6
figure_height <- 4.8
ordination_width <- 7
ordination_height <- 3.8

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

# 2. Check input identity and retain the existing sediment feature filter --------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "hja-FTICRMS.xlsx"),
  file.path(input_dir, "fticr", "fticr_site_by_peak_primary.rds"),
  file.path(input_dir, "fticr", "fticr_peak_metadata.csv.gz"),
  file.path(input_dir, "fticr", "fticr_peak_filter_status.csv.gz"),
  file.path(input_dir, "environment", "sediment_44_environment.csv"),
  file.path(reference_dir, "Transformation_Database_07-2020.csv"),
  file.path(reference_dir, "Transformation_Analysis_by_Sample.R")
)
input_md5 <- unname(tools::md5sum(input_paths))
stopifnot(!anyNA(input_md5))
stopifnot(identical(input_md5[6:7], reference_md5))

intensity <- readRDS(input_paths[2])
metadata <- read_csv(input_paths[3], show_col_types = FALSE)
feature_filter <- read_csv(input_paths[4], show_col_types = FALSE)
environment <- read_csv(input_paths[5],
  col_types = cols(site_code = col_character()))
database <- read_csv(input_paths[6], show_col_types = FALSE) %>%
  transmute(database_row = row_number(), transformation = Name, mass = Mass)
stopifnot(
  nrow(environment) == 44L,
  !anyDuplicated(environment$site_code),
  all(environment$site_code %in% rownames(intensity)),
  identical(colnames(intensity), feature_filter$peak_id[feature_filter$included_primary]),
  !anyDuplicated(database$transformation),
  all(is.finite(database$mass)),
  all(database$mass > 0),
  all(is.finite(intensity)),
  all(intensity >= 0)
)
site_audit <- tibble(site_code = rownames(intensity)) %>%
  mutate(included = site_code %in% environment$site_code,
    reason = if_else(included, "sediment_profile", "outside_sediment_44"))
write_csv(site_audit, file.path(results_dir, "site_inclusion_audit.csv"))
detected <- intensity[environment$site_code, , drop = FALSE] > 0
peaks <- metadata[match(colnames(detected), metadata$peak_id), ]
stopifnot(
  !anyNA(peaks$peak_id),
  !anyDuplicated(peaks$measured_mz),
  all(is.finite(peaks$measured_mz)),
  all(peaks$c13_indicator <= 0),
  all(peaks$measured_mz >= 200 & peaks$measured_mz <= 900),
  all(rowSums(detected) > 0)
)
# Features absent from all sediment profiles cannot form within-sample pairs.
peak_order <- order(peaks$measured_mz[colSums(detected) > 0])
peaks <- peaks[colSums(detected) > 0, ][peak_order, ]
detected <- detected[, peaks$peak_id, drop = FALSE]
database_sorted <- arrange(database, mass)

# 3. Match peak differences once, then require both peaks in each sample ---------
# This helper avoids constructing the full pairwise self-join for every sample.
# findInterval retains every database entry within the published inclusive bounds.
match_differences <- function(masses, transformation_table, tolerance) {
  matched_rows <- vector("list", length(masses))
  matched_rows[[1]] <- tibble(lower_index = integer(), upper_index = integer(),
    database_row = integer(), observed_difference = double())
  for (upper in seq_along(masses)[-1]) {
    lower <- seq_len(upper - 1L)
    difference <- masses[upper] - masses[lower]
    first <- findInterval(difference - tolerance,
      transformation_table$mass, left.open = TRUE) + 1L
    last <- findInterval(difference + tolerance, transformation_table$mass)
    number <- pmax(0L, last - first + 1L)
    if (!any(number > 0)) next
    lower_matches <- rep(lower, number)
    transformation_matches <- sequence(number[number > 0], from = first[number > 0])
    matched_rows[[upper]] <- tibble(
      lower_index = lower_matches,
      upper_index = upper,
      database_row = transformation_table$database_row[transformation_matches],
      observed_difference = masses[upper] - masses[lower_matches]
    )
  }
  bind_rows(matched_rows)
}
matching_start <- proc.time()[["elapsed"]]
assignments <- match_differences(peaks$measured_mz, database_sorted,
  absolute_mass_tolerance) %>%
  left_join(database, by = "database_row") %>%
  mutate(
    pair_id = paste(peaks$peak_id[lower_index], peaks$peak_id[upper_index], sep = ":"),
    lower_peak_id = peaks$peak_id[lower_index],
    upper_peak_id = peaks$peak_id[upper_index],
    mass_error = observed_difference - mass
  ) %>%
  group_by(pair_id) %>%
  mutate(
    candidate_assignments = n(),
    published_assignment = database_row == max(database_row)
  ) %>%
  ungroup()
stopifnot(nrow(assignments) > 0,
  all(assignments$observed_difference + absolute_mass_tolerance >= assignments$mass),
  all(assignments$observed_difference - absolute_mass_tolerance <= assignments$mass))
pairs <- assignments %>%
  filter(published_assignment) %>%
  arrange(lower_index, upper_index)
co_detected <- detected[, pairs$lower_index, drop = FALSE] &
  detected[, pairs$upper_index, drop = FALSE]
pairs$observed_sites <- colSums(co_detected)
co_detected <- co_detected[, pairs$observed_sites > 0, drop = FALSE]
pairs <- filter(pairs, observed_sites > 0)
assignments <- assignments %>%
  inner_join(select(pairs, pair_id, observed_sites), by = "pair_id")
write_csv(select(assignments, pair_id, lower_peak_id, upper_peak_id,
  observed_difference, database_row, transformation, mass, mass_error,
  candidate_assignments, published_assignment, observed_sites),
  file.path(results_dir, "peak_pair_assignments.csv.gz"))
matching_seconds <- proc.time()[["elapsed"]] - matching_start

# 4. Summarize the published last-matching-label assignment ----------------------
# Each unordered pair counts once. Alternate labels remain in the assignment file.
profiles <- matrix(0L, nrow = nrow(detected), ncol = nrow(database),
  dimnames = list(rownames(detected), database$transformation))
degree_rows <- vector("list", nrow(detected))
summary_rows <- vector("list", nrow(detected))
for (site in seq_len(nrow(detected))) {
  present_pairs <- pairs[co_detected[site, ], ]
  profiles[site, ] <- tabulate(present_pairs$database_row, nbins = nrow(database))
  degree <- tabulate(c(present_pairs$lower_index, present_pairs$upper_index),
    nbins = nrow(peaks))
  present <- detected[site, ]
  degree_rows[[site]] <- tibble(
    site_code = rownames(detected)[site],
    peak_id = peaks$peak_id[present],
    measured_mz = peaks$measured_mz[present],
    matched_pair_degree = degree[present]
  )
  n_detected <- sum(present)
  summary_rows[[site]] <- tibble(
    site_code = rownames(detected)[site],
    detected_features = n_detected,
    possible_peak_pairs = choose(n_detected, 2),
    matched_peak_pairs = nrow(present_pairs),
    all_candidate_assignments = sum(present_pairs$candidate_assignments),
    ambiguous_peak_pairs = sum(present_pairs$candidate_assignments > 1),
    transformation_labels = sum(profiles[site, ] > 0),
    connected_features = sum(degree > 0),
    mean_pair_degree = sum(degree) / n_detected,
    matched_pair_fraction = nrow(present_pairs) / choose(n_detected, 2)
  )
}
site_summary <- environment %>%
  select(site_code, site_drainage_area_ha, site_stream_order, site_segment_number) %>%
  left_join(bind_rows(summary_rows), by = "site_code")
write_csv(site_summary, file.path(results_dir, "site_summary.csv"))
write_csv(bind_rows(degree_rows), file.path(results_dir, "peak_degrees.csv.gz"))
profile_table <- bind_cols(database,
  as_tibble(t(profiles), .name_repair = "minimal"))
write_csv(profile_table, file.path(results_dir, "transformation_profiles.csv"))
stopifnot(all(rowSums(profiles) == site_summary$matched_peak_pairs),
  all(site_summary$matched_peak_pairs > 0),
  sum(bind_rows(degree_rows)$matched_pair_degree) == 2 * sum(profiles))

# 5. Verify against the original script on one complete observed sample ---------
# Run the unmodified upstream source in a temporary directory, using the smallest
# complete sediment profile to bound the memory required by its Cartesian join.
validation_site <- which.min(rowSums(detected))
validation_present <- detected[validation_site, ]
validation_masses <- peaks$measured_mz[validation_present]
validation_dir <- tempfile("hja_danczak_validation_")
dir.create(validation_dir)
write_csv(tibble(Mass = sprintf("%.17g", validation_masses), validation = 1L),
  file.path(validation_dir, "Validation_Data.csv"))
write_csv(tibble(Mass = sprintf("%.17g", validation_masses), C13 = 0L),
  file.path(validation_dir, "Validation_Mol.csv"))
invisible(file.copy(input_paths[6], file.path(validation_dir,
  "Transformation_Database_07-2020.csv"), overwrite = TRUE))
validation_environment <- new.env(parent = globalenv())
# Only the upstream hard-coded working directory is redirected; its code is intact.
validation_environment$setwd <- function(dir) base::setwd(validation_dir)
original_directory <- getwd()
reference_start <- proc.time()[["elapsed"]]
reference_output <- tryCatch(
  capture.output(sys.source(input_paths[7], envir = validation_environment)),
  finally = base::setwd(original_directory)
)
reference_seconds <- proc.time()[["elapsed"]] - reference_start
reference_pairs <- validation_environment$Distance_Results
reference_counts <- validation_environment$profiles.of.trans
reference_mass_input <- as.numeric(rownames(validation_environment$data))
stopifnot(identical(reference_mass_input, validation_masses))
expected_pairs <- pairs[co_detected[validation_site, ], ]
reference_keys <- paste(
  match(reference_pairs$peak.y, peaks$measured_mz),
  match(reference_pairs$peak.x, peaks$measured_mz),
  reference_pairs$Trans.name, sep = ":")
expected_keys <- paste(expected_pairs$lower_index, expected_pairs$upper_index,
  expected_pairs$transformation, sep = ":")
profile_difference <- reference_counts[["Sample_validation"]] -
  profiles[validation_site, match(reference_counts$Name, colnames(profiles))]
validation <- tibble(
  check = c("complete_sample_upstream_peak_pairs", "complete_sample_upstream_profiles",
    "profile_totals_equal_pairs", "degree_totals_equal_twice_pairs"),
  passed = c(setequal(reference_keys, expected_keys),
    all(profile_difference == 0),
    all(rowSums(profiles) == site_summary$matched_peak_pairs),
    sum(bind_rows(degree_rows)$matched_pair_degree) == 2 * sum(profiles)),
  site_code = c(rep(rownames(detected)[validation_site], 2), NA, NA),
  upstream_seconds = c(reference_seconds, NA, NA, NA),
  pooled_matching_seconds = c(matching_seconds, NA, NA, NA)
)

# An independent exhaustive comparison also checks overlapping labels and bounds.
check_table <- tibble(database_row = 1:3,
  mass = c(1 - absolute_mass_tolerance, 1, 1 + absolute_mass_tolerance))
check_masses <- c(100, 101, 103)
check_matches <- match_differences(check_masses, check_table, absolute_mass_tolerance)
check_pairs <- expand_grid(lower_index = 1:3, upper_index = 1:3,
  database_row = 1:3) %>%
  filter(upper_index > lower_index) %>%
  mutate(observed_difference = check_masses[upper_index] - check_masses[lower_index],
    mass = check_table$mass[database_row]) %>%
  filter(observed_difference + absolute_mass_tolerance >= mass,
    observed_difference - absolute_mass_tolerance <= mass)
validation <- bind_rows(validation, tibble(
  check = "overlapping_labels_and_inclusive_bounds",
  passed = setequal(paste(check_matches$lower_index, check_matches$upper_index,
    check_matches$database_row), paste(check_pairs$lower_index, check_pairs$upper_index,
    check_pairs$database_row))
), tibble(
  check = "sample_without_mass_matches",
  passed = nrow(match_differences(c(100, 110), check_table,
    absolute_mass_tolerance)) == 0
))
write_csv(validation, file.path(results_dir, "numerical_validation.csv"))
if (!all(validation$passed)) stop("Transformation validation failed; inspect numerical_validation.csv.")

# 6. Describe transformation profiles and their dependence on detection ----------
# PCA uses Hellinger-transformed pair-count proportions, not peak intensities.
hellinger <- sqrt(profiles / rowSums(profiles))
ordination <- prcomp(hellinger, center = TRUE, scale. = FALSE)
variance <- 100 * ordination$sdev^2 / sum(ordination$sdev^2)
scores <- site_summary %>%
  mutate(pc1 = ordination$x[, 1], pc2 = ordination$x[, 2])
write_csv(select(scores, site_code, pc1, pc2),
  file.path(results_dir, "profile_pca_scores.csv"))
write_csv(tibble(axis = seq_along(variance), percent_variance = variance),
  file.path(results_dir, "profile_pca_variance.csv"))
loadings <- bind_cols(database,
  tibble(pc1 = ordination$rotation[, 1], pc2 = ordination$rotation[, 2]))
write_csv(loadings, file.path(results_dir, "profile_pca_loadings.csv"))
ordination_plot <- ggplot(scores,
  aes(pc1, pc2, color = log10(site_drainage_area_ha))) +
  geom_point(size = 2.8) +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  coord_equal() +
  labs(x = sprintf("Transformation PC1 (%.1f%%)", variance[1]),
    y = sprintf("Transformation PC2 (%.1f%%)", variance[2])) +
  plot_theme
ggsave(file.path(figures_dir, "2016_transformation_profile_pca.pdf"),
  ordination_plot, width = ordination_width, height = ordination_height, bg = "white")

top_labels <- names(sort(colMeans(profiles / rowSums(profiles)),
  decreasing = TRUE))[seq_len(min(display_transformations, ncol(profiles)))]
top_profiles <- as_tibble(profiles[, top_labels, drop = FALSE] / rowSums(profiles)) %>%
  mutate(site_code = rownames(profiles)) %>%
  pivot_longer(-site_code, names_to = "transformation", values_to = "pair_fraction") %>%
  mutate(transformation = factor(transformation, levels = rev(top_labels)))
profile_plot <- ggplot(top_profiles, aes(pair_fraction, transformation)) +
  geom_boxplot(fill = "#56B4E9", width = 0.65, outlier.size = 1) +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(x = "Fraction of matched peak pairs", y = NULL) +
  plot_theme
ggsave(file.path(figures_dir, "2016_transformation_common_labels.pdf"),
  profile_plot, width = 7, height = 6.5, bg = "white")

diagnostic_plot <- ggplot(site_summary,
  aes(possible_peak_pairs, matched_peak_pairs)) +
  geom_point(size = 2.8, color = "#0072B2") +
  scale_x_log10(breaks = 10^seq(4, 7, by = 0.5), labels = scales::label_log()) +
  scale_y_log10(breaks = 10^seq(3, 5, by = 0.5), labels = scales::label_log()) +
  labs(x = "Possible detected peak pairs", y = "Matched peak pairs") +
  plot_theme
ggsave(file.path(figures_dir, "2016_transformation_detection_diagnostic.pdf"),
  diagnostic_plot, width = figure_width, height = figure_height, bg = "white")

# 7. Record provenance and confirm that inputs remained unchanged ---------------
integrity <- tibble(
  path = substring(input_paths, nchar(here()) + 2L),
  md5_before = input_md5,
  md5_after = unname(tools::md5sum(input_paths)),
  upstream_commit = upstream_commit
) %>%
  mutate(unchanged = md5_before == md5_after)
write_csv(integrity, file.path(results_dir, "input_integrity.csv"))
stopifnot(all(integrity$unchanged))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(results_dir, "session_info.txt"))
message("Putative transformations: ", nrow(detected), " sediment profiles; ",
  nrow(pairs), " co-detected peak pairs across sites; ",
  sum(colSums(profiles) > 0), " assigned transformation labels. Validation passed.")
