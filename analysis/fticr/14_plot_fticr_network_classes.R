# Compare FT-ICR chemical-class fractions across established network groups.
# Run directly in Positron; retain the original drainage-ordered site figure.

library(here)
library(readr)
library(dplyr)
library(ggplot2)

# 1. Settings and plotting theme.
results_dir <- here("results", "fticr_2016", "river_composition")
plot_results_dir <- file.path(results_dir, "network_class_figures")
figures_dir <- here("figures")
group_levels <- c("Headwater", "Intermediate", "Mainstem")
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
class_order <- c("Lipid", "Protein", "AminoSugar", "Carb", "Lignin", "Tannin", "ConHC", "UnsatHC", "Other")
class_labels <- c("Lipid-like", "Protein-like", "Amino-sugar-like", "Carbohydrate-like",
  "Lignin-like", "Tannin-like", "Condensed-hydrocarbon-like", "Unsaturated-hydrocarbon-like", "Other")
class_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#882255", "#AA4499")
names(class_colors) <- class_labels
plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 12),
    legend.position = "bottom"
  )

# 2. Align the existing class fractions to the audited 44-site network groups.
input_paths <- c(
  file.path(results_dir, "site_class_fractions.csv"),
  here("results", "fticr_2016", "network_dispersion", "site_network_position_audit.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "tables", "primary_peak_properties.csv"),
  file.path(figures_dir, "2016_fticr_site_class_composition.pdf")
)
input_hashes <- unname(tools::md5sum(input_paths))
stopifnot(!anyNA(input_hashes))
class_table <- read_csv(input_paths[1], col_types = cols(site_code = col_character()))
network <- read_csv(input_paths[2], col_types = cols(site_code = col_character()))
incidence <- read_csv(input_paths[3], col_types = cols(site_code = col_character()))
properties <- read_csv(input_paths[4], show_col_types = FALSE)
stopifnot(
  !anyDuplicated(network$site_code),
  !anyDuplicated(incidence$site_code),
  !anyDuplicated(properties$peak_id),
  !anyDuplicated(class_table[c("site_code", "source_class")]),
  setequal(class_table$site_code, network$site_code),
  setequal(incidence$site_code, network$site_code),
  setequal(class_table$source_class, class_order)
)
network$expected_group <- case_when(
  network$site_stream_order %in% headwater_orders ~ "Headwater",
  network$site_stream_order %in% intermediate_orders ~ "Intermediate",
  network$site_stream_order == mainstem_order ~ "Mainstem"
)
stopifnot(!anyNA(network$expected_group), all(network$network_group == network$expected_group))
site_order <- network %>%
  mutate(network_group = factor(network_group, levels = group_levels)) %>%
  arrange(network_group, site_drainage_area_ha, site_code)
reference <- network[match(class_table$site_code, network$site_code), ]
stopifnot(
  all(class_table$site_stream_order == reference$site_stream_order),
  all(abs(class_table$site_drainage_area_ha - reference$site_drainage_area_ha) < 1e-8),
  all(is.finite(class_table$fraction_detected)),
  all(class_table$fraction_detected >= 0 & class_table$fraction_detected <= 1)
)
class_table$network_group <- factor(reference$network_group, levels = group_levels)
class_table$class_label <- factor(class_table$source_class, levels = class_order, labels = class_labels)
class_table$site_code <- factor(class_table$site_code, levels = site_order$site_code)

# 3. Verify that the saved fractions agree with current binary detections.
# Keep workbook classifications and detection rules unchanged.
pa <- as.matrix(incidence[, -1])
properties <- properties[match(colnames(pa), properties$peak_id), ]
stopifnot(
  identical(colnames(pa), properties$peak_id),
  all(pa %in% c(0, 1)),
  all(rowSums(pa) > 0),
  all(properties$molecular_class_source %in% class_order)
)
for (source_name in class_order) {
  selected <- properties$molecular_class_source == source_name
  expected <- rowSums(pa[, selected, drop = FALSE]) / rowSums(pa)
  current <- class_table %>% filter(source_class == source_name)
  expected <- expected[match(as.character(current$site_code), incidence$site_code)]
  stopifnot(nrow(current) == nrow(network), all(abs(current$fraction_detected - expected) < 1e-12))
}
site_sums <- class_table %>% group_by(site_code) %>%
  summarize(fraction_sum = sum(fraction_detected), .groups = "drop")
stopifnot(all(abs(site_sums$fraction_sum - 1) < 1e-12))

# 4. Average site fractions within groups, giving every site equal weight.
# Pooling raw feature counts would give feature-rich sites greater weight.
group_summary <- class_table %>% group_by(network_group, source_class, class_label) %>%
  summarize(
    n_sites = n(),
    mean_fraction = mean(fraction_detected),
    sd_fraction = sd(fraction_detected),
    .groups = "drop"
  )
group_sums <- group_summary %>% group_by(network_group) %>%
  summarize(fraction_sum = sum(mean_fraction), .groups = "drop")
stopifnot(all(abs(group_sums$fraction_sum - 1) < 1e-12))
group_counts <- table(factor(network$network_group, levels = group_levels))
group_labels <- paste0(group_levels, "\n(n = ", as.integer(group_counts), ")")
names(group_labels) <- group_levels

# 5. Draw the compact group comparison and the grouped site-level counterpart.
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
group_plot <- ggplot(group_summary, aes(network_group, mean_fraction, fill = class_label)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Mean fraction of detected features", fill = NULL) +
  plot_theme +
  guides(fill = guide_legend(nrow = 3))
ggsave(file.path(figures_dir, "2016_fticr_network_group_class_comparison.pdf"),
  group_plot, width = 9, height = 5.5, bg = "white")

# Retain the original site bars, ordered by drainage area inside each group.
# Bottom facet labels serve as the group labels on the x axis.
site_panels <- list()
for (group_name in c("All", group_levels)) {
  current <- class_table
  if (group_name != "All") {
    current <- class_table %>% filter(network_group == group_name)
  }
  site_panels[[group_name]] <- ggplot(current, aes(site_code, fraction_detected, fill = class_label)) +
    geom_col(width = 0.95) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25)) +
    facet_grid(. ~ network_group, scales = "free_x", space = "free_x", switch = "x",
      labeller = labeller(network_group = group_labels)) +
    labs(x = NULL, y = "Fraction of detected features", fill = NULL) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    guides(fill = guide_legend(nrow = 3))
  file_suffix <- "by_network_group"
  panel_width <- 11
  if (group_name != "All") {
    file_suffix <- tolower(group_name)
    panel_width <- 9
  }
  ggsave(file.path(figures_dir, paste0("2016_fticr_site_class_composition_", file_suffix, ".pdf")),
    site_panels[[group_name]], width = panel_width, height = 5.5, bg = "white")
}

# 6. Save plot data and confirm that the original figure and inputs are intact.
stopifnot(identical(input_hashes, unname(tools::md5sum(input_paths))))
dir.create(plot_results_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(class_table, file.path(plot_results_dir, "site_class_plot_data.csv"))
write_csv(group_summary, file.path(plot_results_dir, "network_group_class_summary.csv"))
write_csv(tibble(path = input_paths, md5 = input_hashes),
  file.path(plot_results_dir, "preserved_input_figure_checksums.csv"))
message("Saved network-group class comparisons for ", nrow(network), " FT-ICR sites.")
