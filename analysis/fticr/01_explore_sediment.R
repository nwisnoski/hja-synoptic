# Describe sediment molecular profiles before choosing environmental models.

library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(vegan)

# 1. Settings and plotting theme ------------------------------------------------
# The primary feature filter belongs to data_prep/config.R, not this script.
input_dir <- here("data", "derived", "analysis_inputs")
results_dir <- here("results", "fticr_2016", "tables")
derived_dir <- here("data", "derived", "fticr_2016")
figures_dir <- here("figures")
presence_threshold <- 0
pcoa_correction <- "lingoes"
figure_width <- 6
figure_height <- 4.8
raster_dpi <- 600

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    legend.position = "right"
  )

# 2. Read prepared inputs and check their provenance ----------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
inventory <- read_csv(file.path(input_dir, "audit", "source_file_inventory.csv"),
  show_col_types = FALSE)
source_checks <- inventory %>%
  mutate(current_md5 = unname(tools::md5sum(here(source_path))),
    checksum_matches = current_md5 == md5)
write_csv(source_checks, file.path(results_dir, "source_checks.csv"))
stopifnot(all(source_checks$checksum_matches))

intensity_60 <- readRDS(file.path(input_dir, "fticr", "fticr_site_by_peak_primary.rds"))
qc_all <- readRDS(file.path(input_dir, "fticr", "fticr_qc_by_peak.rds"))
peaks_all <- read_csv(file.path(input_dir, "fticr", "fticr_peak_metadata.csv.gz"),
  show_col_types = FALSE)
filter_status <- read_csv(file.path(input_dir, "fticr", "fticr_peak_filter_status.csv.gz"),
  show_col_types = FALSE)
environment <- read_csv(file.path(input_dir, "environment", "sediment_44_environment.csv"),
  col_types = cols(site_code = col_character()))
crosswalk <- read_csv(file.path(input_dir, "audit", "sediment_site_crosswalk.csv"),
  col_types = cols(site_code = col_character()))
microbial_metadata <- read_csv(here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  col_types = cols(site_code = col_character()))

stopifnot(
  !anyDuplicated(environment$site_code),
  !anyDuplicated(rownames(intensity_60)),
  !anyDuplicated(peaks_all$peak_id),
  identical(colnames(intensity_60), filter_status$peak_id[filter_status$included_primary]),
  all(is.finite(intensity_60)),
  all(intensity_60 >= 0),
  setequal(environment$site_code, crosswalk$site_code[crosswalk$include_sediment_multiblock]),
  all(environment$site_code %in% rownames(intensity_60))
)

# Audit all 60 profiles plus any sediment 10K sites without a molecular profile.
sediment_10k <- microbial_metadata %>%
  filter(habitat == "sediment", included_10k) %>%
  select(site_code, sample_id)
stopifnot(!anyDuplicated(sediment_10k$site_code))
site_audit <- full_join(
  tibble(site_code = rownames(intensity_60), has_fticr_profile = TRUE),
  sediment_10k, by = "site_code"
) %>%
  mutate(
    has_fticr_profile = coalesce(has_fticr_profile, FALSE),
    in_sediment_44 = site_code %in% environment$site_code,
    has_sediment_10k = !is.na(sample_id),
    paired_eligible = in_sediment_44 & has_fticr_profile & has_sediment_10k,
    analysis_role = case_when(
      paired_eligible ~ "sediment_44_and_paired_10k",
      in_sediment_44 ~ "sediment_44_without_10k",
      has_sediment_10k ~ "sediment_10k_without_fticr",
      TRUE ~ "fticr_outside_sediment_44"
    )
  )
write_csv(site_audit, file.path(results_dir, "site_overlap_audit.csv"))

intensity <- intensity_60[environment$site_code, , drop = FALSE]
stopifnot(identical(rownames(intensity), environment$site_code), all(rowSums(intensity) > 0))
stopifnot(all(is.finite(environment$site_drainage_area_ha)), all(environment$site_drainage_area_ha > 0))
detected <- intensity > presence_threshold
presence_absence <- bind_cols(tibble(site_code = rownames(detected)), as_tibble(1L * detected))
write_csv(presence_absence, file.path(derived_dir, "sediment_44_presence_absence.csv"))
peak_rows <- match(colnames(intensity), peaks_all$peak_id)
stopifnot(!anyNA(peak_rows))
peaks <- peaks_all[peak_rows, ] %>%
  mutate(
    detected_sites_44 = colSums(detected),
    prevalence_44 = detected_sites_44 / nrow(intensity),
    hydrogen_carbon_ratio = hydrogen_count / carbon_count,
    oxygen_carbon_ratio = oxygen_count / carbon_count
  )
write_csv(peaks, file.path(results_dir, "primary_peak_properties.csv"))
stopifnot(all(is.finite(peaks$hydrogen_carbon_ratio)), all(is.finite(peaks$oxygen_carbon_ratio)))

# Intensities remain QC diagnostics; all molecular composition uses detections.
site_summary <- environment %>%
  mutate(
    detected_primary_features = rowSums(detected),
    total_primary_intensity = rowSums(intensity),
    mean_detected_h_c = as.vector(detected %*% peaks$hydrogen_carbon_ratio) / rowSums(detected),
    mean_detected_o_c = as.vector(detected %*% peaks$oxygen_carbon_ratio) / rowSums(detected),
    has_sediment_10k = site_code %in% sediment_10k$site_code
  )
write_csv(site_summary, file.path(results_dir, "sediment_site_summary.csv"))
qc_primary <- qc_all[, colnames(intensity_60), drop = FALSE]
qc_summary <- tibble(
  qc_sample_id = rownames(qc_all),
  detected_raw_features = rowSums(qc_all > presence_threshold),
  total_raw_intensity = rowSums(qc_all),
  detected_site_filtered_features = rowSums(qc_primary > presence_threshold),
  total_site_filtered_intensity = rowSums(qc_primary)
)
write_csv(qc_summary, file.path(results_dir, "laboratory_standard_summary.csv"))

# 3. Presence/absence composition ----------------------------------------------
jaccard <- vegdist(detected, method = "jaccard", binary = TRUE)
raw_pcoa <- wcmdscale(jaccard, k = 2, eig = TRUE)
pa_pcoa <- wcmdscale(jaccard, k = 2, eig = TRUE, add = pcoa_correction)
pa_percent <- 100 * pa_pcoa$eig / sum(pa_pcoa$eig[pa_pcoa$eig > 0])
pa_scores <- site_summary %>%
  mutate(pcoa1 = pa_pcoa$points[, 1], pcoa2 = pa_pcoa$points[, 2])
write_csv(pa_scores, file.path(results_dir, "jaccard_pcoa_scores.csv"))
ordination_summary <- tibble(
  method = "Jaccard incidence PCoA",
  sites = nrow(intensity),
  features_detected_in_subset = sum(colSums(detected) > 0),
  axis1_percent = pa_percent[1],
  axis2_percent = pa_percent[2],
  correction = pcoa_correction,
  additive_constant = pa_pcoa$ac,
  uncorrected_negative_eigenvalues = sum(raw_pcoa$eig < -1e-10)
)
write_csv(ordination_summary, file.path(results_dir, "ordination_summary.csv"))
signal_diagnostic <- tibble(
  comparison = "Detected primary feature count versus total primary intensity",
  site_spearman_rho = cor(site_summary$detected_primary_features,
    site_summary$total_primary_intensity, method = "spearman")
)
write_csv(signal_diagnostic, file.path(results_dir, "richness_signal_diagnostic.csv"))

# 4. Save individual figures ----------------------------------------------------
pa_plot <- ggplot(pa_scores, aes(pcoa1, pcoa2, color = log10(site_drainage_area_ha))) +
  geom_point(size = 2.8) +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  coord_equal() +
  labs(x = sprintf("Jaccard PCoA1 (%.1f%%)", pa_percent[1]),
    y = sprintf("Jaccard PCoA2 (%.1f%%)", pa_percent[2])) +
  plot_theme
ggsave(file.path(figures_dir, "2016_fticr_sediment_jaccard_pcoa.pdf"),
  pa_plot, width = figure_width, height = figure_height, bg = "white")

richness_plot <- ggplot(site_summary, aes(site_drainage_area_ha, detected_primary_features)) +
  geom_point(size = 2.8, color = "#0072B2") +
  scale_x_log10() +
  labs(x = "Drainage area (ha)", y = "Detected primary molecular features") +
  plot_theme
ggsave(file.path(figures_dir, "2016_fticr_sediment_richness_drainage.pdf"),
  richness_plot, width = figure_width, height = figure_height, bg = "white")

# Keep all detected ratios, without assigning structural or metabolite identities.
vk_plot <- ggplot(filter(peaks, detected_sites_44 > 0),
  aes(oxygen_carbon_ratio, hydrogen_carbon_ratio, color = prevalence_44)) +
  geom_point(size = 0.7, alpha = 0.65) +
  scale_color_viridis_c(name = "Proportion of\n44 sites", limits = c(0, 1)) +
  labs(x = "O:C atomic ratio", y = "H:C atomic ratio") +
  plot_theme
ggsave(file.path(figures_dir, "2016_fticr_sediment_van_krevelen.png"),
  vk_plot, width = figure_width, height = figure_height, dpi = raster_dpi, bg = "white")

# 5. Record the run -------------------------------------------------------------
session_lines <- capture.output(sessionInfo())
writeLines(trimws(session_lines, which = "right"), file.path(results_dir, "session_info.txt"))
message("FT-ICR exploration: ", nrow(intensity), " sites; ", ncol(intensity),
  " primary features; ", sum(site_audit$paired_eligible), " paired sediment 10K sites.")
