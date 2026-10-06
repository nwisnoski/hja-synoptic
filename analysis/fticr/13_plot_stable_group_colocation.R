# Overlay both components of the completed stable group on one catchment map.
# Run directly in Positron; the labeled supplementary maps remain unchanged.

library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. Settings and plotting theme.
core_id <- 2L
expected_asvs <- 38L
expected_molecules <- 8L
results_dir <- here("results", "fticr_2016", "stable_group_distribution")
plot_results_dir <- file.path(results_dir, "colocation_figure")
figures_dir <- here("figures")
figure_file <- "2016_sediment_stable_group_colocation_map.pdf"
labeled_figure_file <- "2016_sediment_stable_group_colocation_map_labeled.pdf"
label_seed <- 20261044L
label_size <- 2.8
marker_size_range <- c(2, 8)
fraction_breaks <- c(0, 0.5, 1)
fraction_labels <- c("0%", "50%", "100%")
centerline_stride <- 10L
figure_width <- 8.5
figure_height <- 7
plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 2. Check that the completed site profiles still describe the current inputs.
reference_checks <- read_csv(file.path(results_dir, "input_integrity.csv"),
  show_col_types = FALSE)
reference_checks$path <- sub("^.*?/hja-synoptic/", "", reference_checks$path)
reference_checks$current_md5 <- unname(tools::md5sum(here(reference_checks$path)))
if (anyNA(reference_checks$current_md5) ||
    any(reference_checks$current_md5 != reference_checks$final_md5)) {
  stop("Stable-group inputs have changed; review the distribution analysis before plotting.")
}
preserved_paths <- c(
  file.path(results_dir, "core_site_profiles.csv"),
  file.path(results_dir, "site_network_audit.csv"),
  file.path(results_dir, "settings.csv"),
  here("data", "StreamCenterline.csv"),
  file.path(figures_dir, paste0("2016_sediment_stable_group_",
    c("asv_map", "molecular_map", "catchment_maps"), ".pdf"))
)
preserved_hashes <- unname(tools::md5sum(preserved_paths))
stopifnot(!anyNA(preserved_hashes))
site <- read_csv(preserved_paths[1], col_types = cols(site_code = col_character()))
network <- read_csv(preserved_paths[2], col_types = cols(site_code = col_character()))
settings <- read_csv(preserved_paths[3], show_col_types = FALSE)
geometry <- read_csv(preserved_paths[4], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
stopifnot(
  settings$core == core_id,
  settings$original_asvs == expected_asvs,
  settings$molecular_features == expected_molecules,
  !anyDuplicated(site$site_code),
  !anyDuplicated(network$site_code),
  setequal(site$site_code, network$site_code),
  nrow(problems(geometry)) == 0,
  all(is.finite(geometry$x)),
  all(is.finite(geometry$y))
)

# 3. Preserve nondetections and unavailable profiles as different observations.
plot_data <- site %>% mutate(
  paired_profiles_available = microbial_profile_available & molecular_profile_available
)
reference <- network[match(plot_data$site_code, network$site_code), ]
stopifnot(
  !anyNA(plot_data$paired_profiles_available),
  identical(plot_data$microbial_profile_available, !is.na(plot_data$core_asvs_detected)),
  identical(plot_data$molecular_profile_available, !is.na(plot_data$core_molecules_detected)),
  identical(plot_data$paired_profiles_available, plot_data$inference_site),
  all(abs(plot_data$site_utm_x_m - reference$site_utm_x_m) < 1e-8),
  all(abs(plot_data$site_utm_y_m - reference$site_utm_y_m) < 1e-8),
  all(plot_data$site_segment_number == reference$site_segment_number)
)
paired <- plot_data %>% filter(paired_profiles_available) %>%
  arrange(desc(core_asv_fraction))
unpaired <- plot_data %>% filter(!paired_profiles_available)
stopifnot(
  nrow(paired) == settings$paired_sites,
  all(is.finite(paired$core_asv_fraction)),
  all(is.finite(paired$core_molecular_fraction)),
  all(paired$core_asv_fraction >= 0 & paired$core_asv_fraction <= 1),
  all(paired$core_molecular_fraction >= 0 & paired$core_molecular_fraction <= 1),
  all(abs(paired$core_asv_fraction - paired$core_asvs_detected / expected_asvs) < 1e-12),
  all(abs(paired$core_molecular_fraction - paired$core_molecules_detected / expected_molecules) < 1e-12)
)
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% centerline_stride == 1L | row_number() == n()) %>%
  ungroup()

# 4. Map both original detection fractions at the same coordinates.
# A minimum circle size keeps zero-ASV observations visible, including sites
# with molecular detections but no ASVs. Crosses are missing measurements.
# This overlay does not create a joint score or a module-presence threshold.
map <- ggplot(map_lines, aes(x, y, group = segment)) +
  geom_path(color = "#C0C0C0", linewidth = 0.35) +
  geom_point(
    data = paired,
    aes(site_utm_x_m, site_utm_y_m,
      size = core_asv_fraction, fill = core_molecular_fraction),
    inherit.aes = FALSE,
    shape = 21,
    color = "black",
    stroke = 0.4
  ) +
  geom_point(
    data = unpaired,
    aes(site_utm_x_m, site_utm_y_m, shape = "One profile\nunavailable"),
    inherit.aes = FALSE,
    color = "#777777",
    size = 2.5,
    stroke = 0.7
  ) +
  scale_size_continuous(
    range = marker_size_range,
    limits = c(0, 1),
    breaks = fraction_breaks,
    labels = fraction_labels,
    name = "ASVs detected\n(% of 38)"
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    limits = c(0, 1),
    breaks = fraction_breaks,
    labels = fraction_labels,
    name = "FT-ICR features detected\n(% of 8)"
  ) +
  scale_shape_manual(values = c("One profile\nunavailable" = 4), name = NULL) +
  guides(
    size = guide_legend(order = 1, direction = "horizontal", title.position = "top",
      override.aes = list(fill = "grey75")),
    fill = guide_colorbar(order = 2, direction = "horizontal", title.position = "top"),
    shape = guide_legend(order = 3)
  ) +
  coord_equal() +
  labs(x = "UTM easting (m)", y = "UTM northing (m)") +
  plot_theme
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figures_dir, figure_file), map,
  width = figure_width, height = figure_height, bg = "white")

# Retain all site IDs with the original label seed and leader-line styling.
# Reserve room for the largest circles so labels do not cover their fills.
labeled_map <- map +
  geom_text_repel(
    data = plot_data,
    aes(site_utm_x_m, site_utm_y_m, label = site_code),
    inherit.aes = FALSE,
    seed = label_seed,
    size = label_size,
    point.size = max(marker_size_range),
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "#999999",
    segment.size = 0.2
  )
ggsave(file.path(figures_dir, labeled_figure_file), labeled_map,
  width = figure_width, height = figure_height, bg = "white")

# 5. Save the map data and verify that supplementary figures remain intact.
stopifnot(identical(preserved_hashes, unname(tools::md5sum(preserved_paths))))
dir.create(plot_results_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(plot_data, file.path(plot_results_dir, "site_plot_data.csv"))
write_csv(tibble(path = preserved_paths, md5 = preserved_hashes),
  file.path(plot_results_dir, "preserved_input_figure_checksums.csv"))
message("Saved labeled and unlabeled single-panel overlays: ", nrow(paired), " paired sites and ",
  nrow(unpaired), " sites with one profile unavailable; supplementary maps unchanged.")
