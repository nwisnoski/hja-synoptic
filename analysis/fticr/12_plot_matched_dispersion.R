# Plot the completed sediment microbial and FT-ICR dispersion comparison.
# Run directly in Positron; this script does not refit dispersion or its tests.

library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme.
results_dir <- here("results", "fticr_2016", "matched_sediment_dispersion")
plot_results_dir <- file.path(results_dir, "figure_comparison")
figures_dir <- here("figures")
figure_prefix <- "2016_sediment_microbes_fticr_dispersion"
plot_seed <- 20261017
group_levels <- c("Headwater", "Intermediate", "Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
panel_design <- tibble(
  metric = c("microbial_bray", "molecular_equal_features"),
  file_label = c("microbial", "molecular"),
  axis_label = c("Sediment microbes\nBray-Curtis", "FT-ICR molecular signatures\nEqual-feature Jaccard")
)
panel_width <- 5.5
panel_height <- 4.5
plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 12)),
    plot.tag = element_text(face = "bold")
  )

# 2. Check that the completed matched analysis still describes these inputs.
reference_checks <- read_csv(
  file.path(results_dir, "preserved_input_output_checksums.csv"),
  show_col_types = FALSE
)
# The original audit records absolute paths; resolve its project paths here.
reference_checks$path <- sub("^.*?/hja-synoptic/", "", reference_checks$path)
reference_checks$current_md5 <- unname(tools::md5sum(here(reference_checks$path)))
if (anyNA(reference_checks$current_md5) ||
    any(reference_checks$current_md5 != reference_checks$final_md5)) {
  stop("Inputs to the matched analysis have changed; review that analysis before plotting.")
}
source_paths <- file.path(results_dir, c(
  "site_dispersion.csv", "matched_site_manifest.csv", "three_group_tests.csv"
))
source_hashes <- unname(tools::md5sum(source_paths))
dispersion <- read_csv(source_paths[1], col_types = cols(site_code = col_character()))
manifest <- read_csv(source_paths[2], col_types = cols(site_code = col_character()))
tests <- read_csv(source_paths[3], show_col_types = FALSE)
stopifnot(
  !anyDuplicated(manifest$site_code),
  !anyDuplicated(manifest$sample_id),
  !anyDuplicated(tests$metric)
)

# 3. Retain the same sites and group definitions in both panels.
plot_data <- dispersion %>%
  filter(metric %in% panel_design$metric) %>%
  mutate(network_group = factor(network_group, levels = group_levels))
for (metric_name in panel_design$metric) {
  current <- plot_data %>% filter(metric == metric_name)
  expected <- manifest[match(current$site_code, manifest$site_code), ]
  stopifnot(
    !anyDuplicated(current$site_code),
    setequal(current$site_code, manifest$site_code),
    !anyNA(current$network_group),
    all(is.finite(current$distance)),
    all(current$distance >= 0),
    identical(current$sample_id, expected$sample_id),
    all(as.character(current$network_group) == expected$network_group),
    all(current$site_segment_number == expected$site_segment_number),
    all(current$site_stream_order == expected$site_stream_order)
  )
}
group_counts <- table(factor(manifest$network_group, levels = group_levels))
x_labels <- paste0(group_levels, "\n(n = ", as.integer(group_counts), ")")
names(x_labels) <- group_levels
group_means <- plot_data %>%
  group_by(metric, network_group) %>%
  summarize(n_sites = n(), mean_distance = mean(distance), .groups = "drop")
for (metric_name in panel_design$metric) {
  current_means <- group_means %>% filter(metric == metric_name)
  recorded <- tests %>% filter(metric == metric_name)
  expected_means <- as.numeric(unlist(recorded[c(
    "headwater_mean", "intermediate_mean", "mainstem_mean"
  )]))
  stopifnot(
    nrow(recorded) == 1L,
    recorded$sites == nrow(manifest),
    all(current_means$n_sites == as.integer(group_counts)),
    all(abs(current_means$mean_distance - expected_means) < 1e-12)
  )
}

# 4. Plot site distributions and group means on their original metric scales.
# Black diamonds are means; boxes describe the observed site distribution.
# No trend lines or new tests are introduced for these categorical groups.
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
panels <- list()
for (i in seq_len(nrow(panel_design))) {
  metric_name <- panel_design$metric[i]
  current <- plot_data %>% filter(metric == metric_name)
  current_means <- group_means %>% filter(metric == metric_name)
  panels[[i]] <- ggplot(current, aes(network_group, distance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.5) +
    geom_point(
      position = position_jitter(width = 0.09, height = 0, seed = plot_seed),
      size = 2,
      alpha = 0.85
    ) +
    geom_point(
      data = current_means,
      aes(network_group, mean_distance),
      inherit.aes = FALSE,
      shape = 18,
      size = 3.3,
      color = "black"
    ) +
    scale_color_manual(values = group_colors) +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    labs(
      x = panel_design$axis_label[i],
      y = "Distance to group spatial median"
    ) +
    plot_theme
  ggsave(
    file.path(figures_dir, paste0(figure_prefix, "_", panel_design$file_label[i], ".pdf")),
    panels[[i]],
    width = panel_width,
    height = panel_height,
    bg = "white"
  )
}
comparison <- wrap_plots(panels, nrow = 1) + plot_annotation(tag_levels = "A")
ggsave(
  file.path(figures_dir, paste0(figure_prefix, "_comparison.pdf")),
  comparison,
  width = 2 * panel_width,
  height = panel_height,
  bg = "white"
)

# 5. Save inspectable plot data and confirm the analysis inputs were untouched.
stopifnot(identical(source_hashes, unname(tools::md5sum(source_paths))))
dir.create(plot_results_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(plot_data, file.path(plot_results_dir, "site_plot_data.csv"))
write_csv(group_means, file.path(plot_results_dir, "group_means.csv"))
write_csv(
  tibble(path = basename(source_paths), md5 = source_hashes),
  file.path(plot_results_dir, "source_checksums.csv")
)
message("Saved matched microbial and FT-ICR dispersion panels without refitting the analysis.")
