# Plot catchment-wide representation of soil-detected ASVs in aquatic samples.
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(mgcv)

# 1. Settings and plotting theme ------------------------------------------------
target_depth <- 10000
smooth_basis_dimension <- 4
reference_areas_ha <- c(100, 1000)
y_metric <- "soil_pool_read_percent"
# Alternative: "soil_pool_richness_fraction" is the fraction of each aquatic
# sample's detected ASVs that were also detected in at least one sampled soil.
habitat_levels <- c("sediment", "planktonic", "hyporheic")
habitat_labels <- c(sediment = "Sediment", planktonic = "Planktonic", hyporheic = "Hyporheic")
habitat_colors <- c(sediment = "#009E73", planktonic = "#0072B2", hyporheic = "#E69F00")
results_dir <- here("results", "diversity_2016", "soil_stream_localization")
figure_path <- here("figures", "2016_soil_asv_catchment_pattern.pdf")
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    strip.background = element_blank(), legend.position = "top")
theme_set(plot_theme)

# 2. Read the completed sharing table and validate its richness denominator ------
input_paths <- c(
  file.path(results_dir, "site_soil_asv_sharing.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
sharing <- read_csv(input_paths[1], col_types = cols(site_code = col_character()), show_col_types = FALSE) %>%
  filter(pool == "soil_detected")
count_table <- read_csv(input_paths[2], show_col_types = FALSE)
metadata <- read_csv(input_paths[3], col_types = cols(site_code = col_character()), show_col_types = FALSE)
stopifnot(!anyDuplicated(sharing$sample_id), all(sharing$sample_id %in% names(count_table)),
  all(sharing$habitat %in% habitat_levels), all(sharing$site_drainage_area_ha > 0))
aquatic_counts <- as.matrix(count_table[, sharing$sample_id])
soil_ids <- metadata$sample_id[metadata$included_10k & metadata$habitat == "soil"]
soil_detected <- rowSums(as.matrix(count_table[, soil_ids])) > 0
stopifnot(length(soil_ids) == 15, all(colSums(aquatic_counts) == target_depth),
  all(colSums(aquatic_counts[soil_detected, , drop = FALSE] > 0) == sharing$soil_pool_asvs_detected),
  max(abs(100 * colSums(aquatic_counts[soil_detected, , drop = FALSE]) / target_depth -
    sharing$soil_pool_read_percent)) < 1e-10)
plot_data <- sharing %>% mutate(aquatic_observed_asvs = colSums(aquatic_counts > 0),
  soil_pool_richness_fraction = soil_pool_asvs_detected / aquatic_observed_asvs,
  habitat = factor(habitat, levels = habitat_levels))
stopifnot(nrow(plot_data) == 85, all(table(plot_data$habitat) == c(34, 20, 31)),
  all(plot_data$soil_pool_richness_fraction >= 0 & plot_data$soil_pool_richness_fraction <= 1),
  y_metric %in% c("soil_pool_read_percent", "soil_pool_richness_fraction"))
write_csv(plot_data, file.path(results_dir, "catchment_soil_asv_plot_data.csv"))

# 3. Fit a small penalized smooth separately within each habitat -----------------
# REML can reduce each four-dimensional basis to an approximately linear fit.
# These are descriptive curves: shared segments preclude independent-sample
# inference, so neither P values nor confidence bands are presented.
prediction_rows <- list()
fit_rows <- list()
for (habitat_name in habitat_levels) {
  current <- plot_data %>% filter(habitat == habitat_name) %>%
    mutate(log_area = log10(site_drainage_area_ha), response = .data[[y_metric]])
  fit <- gam(response ~ s(log_area, k = smooth_basis_dimension), data = current, method = "REML")
  prediction <- tibble(log_area = seq(min(current$log_area), max(current$log_area), length.out = 200))
  prediction$fit <- as.numeric(predict(fit, newdata = prediction))
  prediction_rows[[habitat_name]] <- prediction %>%
    mutate(site_drainage_area_ha = 10^log_area, habitat = habitat_name)
  stopifnot(all(log10(reference_areas_ha) >= min(current$log_area)),
    all(log10(reference_areas_ha) <= max(current$log_area)), fit$converged)
  reference_fit <- as.numeric(predict(fit, newdata = tibble(log_area = log10(reference_areas_ha))))
  fit_rows[[habitat_name]] <- tibble(habitat = habitat_name, y_metric = y_metric,
    samples = nrow(current), segments = n_distinct(current$site_segment_number),
    smooth_edf = unname(summary(fit)$s.table[1, "edf"]),
    fitted_at_100_ha = reference_fit[1], fitted_at_1000_ha = reference_fit[2],
    change_100_to_1000_ha = reference_fit[2] - reference_fit[1])
}
predictions <- bind_rows(prediction_rows)
fit_summary <- bind_rows(fit_rows)
stopifnot(all(is.finite(predictions$fit)))
write_csv(predictions, file.path(results_dir, "catchment_soil_asv_fitted_curves.csv"))
write_csv(fit_summary, file.path(results_dir, "catchment_soil_asv_fit_summary.csv"))

# 4. Plot samples and fitted curves within each habitat's observed area range -----
# A logarithmic area axis keeps headwaters visible across the full catchment.
plot <- ggplot(plot_data, aes(site_drainage_area_ha, .data[[y_metric]], color = habitat)) +
  geom_line(data = predictions, aes(site_drainage_area_ha, fit, color = habitat),
    inherit.aes = FALSE, linewidth = 1) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(values = habitat_colors, labels = habitat_labels, drop = FALSE) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = scales::label_comma()) +
  labs(x = "Drainage area (ha; log scale)", color = NULL)
if (y_metric == "soil_pool_read_percent") {
  plot <- plot + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(y = "Soil-detected ASV reads (%)")
} else {
  plot <- plot + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(y = "Soil-detected fraction of aquatic ASV richness")
}
ggsave(figure_path, plot, width = 8, height = 5.5, bg = "white")

# 5. Retain plotting choices and input integrity ---------------------------------
write_csv(tibble(y_metric = y_metric, x_transform = "log10", aquatic_samples = nrow(plot_data),
  soil_samples = length(soil_ids), soil_detected_asvs = sum(soil_detected),
  fit_method = "Gaussian GAM; REML", smooth_basis_dimension = smooth_basis_dimension),
  file.path(results_dir, "catchment_soil_asv_figure_settings.csv"))
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "catchment_figure_input_checksums.csv"))
if (interactive()) {
  print(plot)
}
