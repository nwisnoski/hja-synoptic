# Reframe completed analyses for the sediment-only manuscript without new tests.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(vegan)

# 1. Settings and one shared plotting style -------------------------------------
results_dir <- here("results", "sediment_manuscript_figures_2016")
figures_dir <- here("figures", "supplement")
manuscript_dir <- here("figures", "manuscript")
figure_prefix <- "2016_sediment_ms_"
plot_seed <- 20261090
centerline_stride <- 10
group_levels <- c("Headwater", "Intermediate", "Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
profile_levels <- c("Both profiles (33)", "Microbial only (1)", "Molecular only (11)")
profile_shapes <- c(21, 24, 22)
class_order <- c("Lipid", "Protein", "AminoSugar", "Carb", "Lignin", "Tannin", "ConHC", "UnsatHC", "Other")
class_labels <- c("Lipid-like", "Protein-like", "Amino-sugar-like", "Carbohydrate-like",
  "Lignin-like", "Tannin-like", "Condensed-hydrocarbon-like", "Unsaturated-hydrocarbon-like", "Other")
class_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#882255", "#AA4499")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
    strip.background = element_blank(), legend.position = "bottom",
    plot.tag = element_text(face = "bold"), plot.margin = margin(8, 8, 8, 8))
dispersion_design <- tribble(
  ~key, ~pool, ~metric, ~axis_label,
  "dispersion_microbes_all", "microbial", "bray", "Sediment microbes\nBray-Curtis",
  "dispersion_molecules_complete_all", "molecular", "jaccard", "Sediment molecular signatures\nComplete-profile Jaccard",
  "dispersion_molecules_equal_all", "molecular", "equal_features", "Sediment molecular signatures\nEqual-feature sensitivity",
  "dispersion_microbes_paired", "paired", "microbial_bray", "Paired sediment microbes\nBray-Curtis",
  "dispersion_molecules_paired", "paired", "molecular_equal_features", "Paired sediment molecular signatures\nEqual-feature sensitivity")
panel_design <- tribble(
  ~key, ~width, ~height, ~samples, ~role,
  "sampling_map", 8, 7, "45 unique: 34 microbial, 44 molecular, 33 paired", "supporting",
  "microbial_pca", 5.8, 5, "34 microbial", "individual",
  "molecular_pcoa", 5.8, 5, "44 molecular", "supporting",
  "molecular_classes", 6.6, 5, "44 molecular", "individual",
  "composition", 12, 5.4, "34 microbial; 44 molecular", "supporting",
  "diversity", 6, 4.7, "34 microbial", "supporting",
  "dispersion_microbes_all", 5.4, 4.6, "34 microbial", "individual",
  "dispersion_molecules_complete_all", 5.4, 4.6, "44 molecular", "individual",
  "dispersion_molecules_equal_all", 5.4, 4.6, "44 molecular", "individual",
  "dispersion_all_samples", 14.5, 4.6, "34 microbial; 44 molecular", "supporting",
  "dispersion_microbes_paired", 5.4, 4.6, "33 paired", "individual",
  "dispersion_molecules_paired", 5.4, 4.6, "33 paired", "individual",
  "dispersion_paired", 10.8, 4.6, "33 identical sites", "main",
  "characteristic_associations", 8, 4.3, "Microbial 34/34/33; molecular 44/44/41", "individual",
  "paired_property_adjustment", 5.8, 4.3, "33 paired; 32 conductivity sensitivity", "individual",
  "characteristics", 14, 4.5, "All eligible samples per model", "supporting",
  "microbial_pca_organic", 6, 5, "34 microbial", "supporting",
  "unifrac_weighted", 5.8, 5, "34 microbial", "individual",
  "unifrac_unweighted", 5.8, 5, "34 microbial", "individual",
  "unifrac", 11.6, 5, "34 microbial in each panel", "supporting",
  "molecular_pcoa_equal", 5.8, 5, "44 molecular", "supporting")

# 2. Read completed results and checksum all inputs and pre-existing figures -----
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
paths <- c(
  microbial_summary = here("results", "diversity_2016", "tables", "sediment_10k_figure_summary.csv"),
  microbial_pca = here("results", "diversity_2016", "tables", "sediment_10k_pca_summary.csv"),
  site_union = here("results", "fticr_2016", "stable_group_distribution", "core_site_profiles.csv"),
  geometry = here("data", "StreamCenterline.csv"),
  metadata = here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  molecular_scores = here("results", "fticr_2016", "tables", "jaccard_pcoa_scores.csv"),
  molecular_ordination = here("results", "fticr_2016", "tables", "ordination_summary.csv"),
  classes = here("results", "fticr_2016", "river_composition", "site_class_fractions.csv"),
  microbial_dispersion = here("results", "diversity_2016", "network_dispersion", "site_dispersion.csv"),
  microbial_tests = here("results", "diversity_2016", "network_dispersion", "three_group_tests.csv"),
  molecular_dispersion = here("results", "fticr_2016", "network_dispersion", "site_dispersion.csv"),
  molecular_tests = here("results", "fticr_2016", "network_dispersion", "three_group_tests.csv"),
  paired_dispersion = here("results", "fticr_2016", "matched_sediment_dispersion", "site_dispersion.csv"),
  paired_tests = here("results", "fticr_2016", "matched_sediment_dispersion", "three_group_tests.csv"),
  paired_manifest = here("results", "fticr_2016", "matched_sediment_dispersion", "matched_site_manifest.csv"),
  characteristics = here("results", "fticr_2016", "sediment_characteristics", "composition_characteristic_tests.csv"),
  paired_properties = here("results", "fticr_2016", "sediment_characteristics", "paired_molecular_property_tests.csv"),
  weighted = here("results", "diversity_2016", "tables", "unifrac_weighted_10k.csv"),
  unweighted = here("results", "diversity_2016", "tables", "unifrac_unweighted_10k.csv"),
  equal_scores = here("results", "fticr_2016", "river_composition", "sensitivity_pcoa_scores.csv"),
  equal_ordination = here("results", "fticr_2016", "river_composition", "sensitivity_pcoa_summary.csv"))
old_figures <- list.files(here("figures"), pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE)
old_figures <- old_figures[!startsWith(basename(old_figures), figure_prefix)]
preserved_paths <- unique(c(unname(paths), old_figures))
checks <- tibble(path = preserved_paths, initial_md5 = unname(tools::md5sum(preserved_paths)))
stopifnot(!anyNA(checks$initial_md5))
tables <- list()
for (key in setdiff(names(paths), "geometry")) {
  tables[[key]] <- read_csv(paths[[key]], show_col_types = FALSE)
  if ("site_code" %in% names(tables[[key]])) tables[[key]]$site_code <- as.character(tables[[key]]$site_code)
}
site <- tables$site_union %>% mutate(network_group = factor(network_group, levels = group_levels),
  coverage = case_when(microbial_profile_available & molecular_profile_available ~ profile_levels[1],
    microbial_profile_available ~ profile_levels[2], TRUE ~ profile_levels[3]),
  coverage = factor(coverage, levels = profile_levels))
stopifnot(nrow(site) == 45, !anyDuplicated(site$site_code),
  sum(site$microbial_profile_available) == 34, sum(site$molecular_profile_available) == 44,
  sum(site$inference_site) == 33, all(site$centerline_match_error_m < 1e-5),
  all(site$site_stream_order == site$centerline_stream_order))
microbial <- tables$microbial_summary %>% left_join(site %>% select(site_code, network_group), by = "site_code")
metadata <- tables$metadata %>% filter(habitat == "sediment", included_10k)
stopifnot(nrow(microbial) == 34, setequal(microbial$sample_id, metadata$sample_id))
microbial <- microbial %>% left_join(metadata %>% select(sample_id, sediment_organic_content_percent), by = "sample_id")
molecular <- tables$molecular_scores %>% select(site_code, pcoa1, pcoa2) %>%
  left_join(site %>% select(site_code, network_group), by = "site_code")
stopifnot(nrow(molecular) == 44, setequal(molecular$site_code, site$site_code[site$molecular_profile_available]))
panels <- list()

# 3. Show all sediment sampling locations and profile availability ----------------
geometry <- read_csv(paths[["geometry"]], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
stopifnot(nrow(problems(geometry)) == 0)
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% centerline_stride == 1 | row_number() == n()) %>% ungroup()
panels$sampling_map <- ggplot(map_lines, aes(x / 1000, y / 1000, group = segment)) +
  geom_path(color = "#BFBFBF", linewidth = 0.35) +
  geom_point(data = site, aes(site_utm_x_m / 1000, site_utm_y_m / 1000,
    fill = network_group, shape = coverage), inherit.aes = FALSE,
    size = 3.3, color = "black", stroke = 0.35) +
  geom_text_repel(data = site %>% filter(site_code == "43"),
    aes(site_utm_x_m / 1000, site_utm_y_m / 1000, label = site_code),
    inherit.aes = FALSE, seed = plot_seed, size = 3, min.segment.length = 0,
    nudge_x = 0.6, nudge_y = 0.3, box.padding = 0.4) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Network position") +
  scale_shape_manual(values = setNames(profile_shapes, profile_levels), drop = FALSE,
    name = "Sediment data") + coord_equal() + labs(x = "Easting (km)", y = "Northing (km)") +
  plot_theme + guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)),
    shape = guide_legend(order = 2)) + theme(legend.box = "vertical")

# 4. Reuse the sediment PCA and molecular summaries with consistent group colors -
pca_percent <- as.numeric(tables$microbial_pca[1, c("pc1_percent", "pc2_percent")])
panels$microbial_pca <- ggplot(microbial, aes(pc1, pc2, color = network_group)) +
  geom_point(size = 2.8) + scale_color_manual(values = group_colors, drop = FALSE, name = NULL) +
  coord_equal() + labs(x = sprintf("Microbial Hellinger PC1 (%.1f%%)", pca_percent[1]),
    y = sprintf("Microbial Hellinger PC2 (%.1f%%)", pca_percent[2])) + plot_theme
panels$molecular_pcoa <- ggplot(molecular, aes(pcoa1, pcoa2, color = network_group)) +
  geom_point(size = 2.8) + scale_color_manual(values = group_colors, drop = FALSE, name = NULL) +
  coord_equal() + labs(x = sprintf("Molecular Jaccard PCo1 (%.1f%%)", tables$molecular_ordination$axis1_percent),
    y = sprintf("Molecular Jaccard PCo2 (%.1f%%)", tables$molecular_ordination$axis2_percent)) + plot_theme
class_data <- tables$classes %>% select(site_code, source_class, fraction_detected) %>%
  left_join(site %>% select(site_code, network_group), by = "site_code") %>%
  mutate(class_label = factor(source_class, levels = class_order, labels = class_labels))
stopifnot(nrow(class_data) == 44 * length(class_order),
  all(abs((class_data %>% group_by(site_code) %>% summarize(total = sum(fraction_detected)))$total - 1) < 1e-10))
class_summary <- class_data %>% group_by(network_group, class_label) %>%
  summarize(mean_fraction = mean(fraction_detected), sites = n(), .groups = "drop")
class_counts <- site %>% filter(molecular_profile_available) %>% count(network_group)
class_labels_x <- setNames(paste0(class_counts$network_group, "\n(n = ", class_counts$n, ")"), class_counts$network_group)
panels$molecular_classes <- ggplot(class_summary, aes(network_group, mean_fraction, fill = class_label)) +
  geom_col(width = 0.65) + scale_fill_manual(values = setNames(class_colors, class_labels), drop = FALSE, name = NULL) +
  scale_x_discrete(labels = class_labels_x) + scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Mean fraction of detected molecular features") + plot_theme +
  guides(fill = guide_legend(ncol = 3)) + theme(legend.text = element_text(size = 8))
panels$composition <- wrap_plots(panels$microbial_pca, panels$molecular_classes, nrow = 1) +
  plot_annotation(tag_levels = "A")
panels$diversity <- ggplot(microbial, aes(site_drainage_area_ha, q1, color = network_group)) +
  geom_point(size = 2.8) + scale_x_log10() + scale_color_manual(values = group_colors, name = NULL) +
  labs(x = "Drainage area (ha)", y = "Sediment Hill diversity (q = 1), 10K reads") + plot_theme

# 5. Replot full-sample and strictly paired dispersion without refitting tests ----
dispersion_tables <- list(microbial = tables$microbial_dispersion %>% filter(habitat == "sediment"),
  molecular = tables$molecular_dispersion, paired = tables$paired_dispersion)
dispersion_rows <- list()
for (j in seq_len(nrow(dispersion_design))) {
  pool <- dispersion_design$pool[j]
  metric <- dispersion_design$metric[j]
  current <- dispersion_tables[[pool]] %>% filter(.data$metric == .env$metric) %>%
    mutate(network_group = factor(network_group, levels = group_levels))
  expected_sites <- site$site_code[site$microbial_profile_available]
  if (pool == "molecular") expected_sites <- site$site_code[site$molecular_profile_available]
  if (pool == "paired") expected_sites <- tables$paired_manifest$site_code
  stopifnot(setequal(current$site_code, expected_sites), !anyDuplicated(current$site_code),
    all(is.finite(current$distance)), !anyNA(current$network_group))
  group_summary <- current %>% group_by(network_group) %>%
    summarize(sites = n(), mean_distance = mean(distance), .groups = "drop")
  recorded <- tables$paired_tests %>% filter(.data$metric == .env$metric)
  if (pool == "microbial") recorded <- tables$microbial_tests %>% filter(habitat == "sediment", .data$metric == .env$metric)
  if (pool == "molecular") recorded <- tables$molecular_tests %>% filter(.data$metric == .env$metric)
  expected_means <- as.numeric(recorded[1, c("headwater_mean", "intermediate_mean", "mainstem_mean")])
  stopifnot(nrow(recorded) == 1, max(abs(group_summary$mean_distance - expected_means)) < 1e-12)
  x_labels <- setNames(paste0(group_summary$network_group, "\n(n = ", group_summary$sites, ")"), group_summary$network_group)
  panels[[dispersion_design$key[j]]] <- ggplot(current, aes(network_group, distance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.45) +
    geom_point(position = position_jitter(width = 0.09, height = 0, seed = plot_seed), size = 2, alpha = 0.85) +
    geom_point(data = group_summary, aes(network_group, mean_distance), inherit.aes = FALSE,
      shape = 18, color = "black", size = 3) +
    scale_color_manual(values = group_colors, guide = "none") + scale_x_discrete(labels = x_labels) +
    labs(x = dispersion_design$axis_label[j], y = "Distance to group spatial median") + plot_theme
  dispersion_rows[[j]] <- current %>% mutate(figure_panel = dispersion_design$key[j])
}
panels$dispersion_all_samples <- wrap_plots(panels$dispersion_microbes_all,
  panels$dispersion_molecules_complete_all, panels$dispersion_molecules_equal_all, nrow = 1) +
  plot_annotation(tag_levels = "A")
panels$dispersion_paired <- wrap_plots(panels$dispersion_microbes_paired,
  panels$dispersion_molecules_paired, nrow = 1) + plot_annotation(tag_levels = "A")

# 6. Display sediment-condition models and the paired adjustment explicitly -------
associations <- tables$characteristics %>% filter(role == "primary") %>%
  mutate(characteristic = factor(comparison, levels = c("organic", "hydraulic", "enzyme"),
    labels = c("Organic content", "Hydraulic conductivity", "Enzyme profile")),
    response_label = factor(response, levels = c("microbial_hellinger", "molecular_jaccard"),
      labels = c("Microbial composition", "Molecular composition")))
stopifnot(nrow(associations) == 6,
  all(associations$sites[associations$response == "microbial_hellinger" & associations$comparison != "hydraulic"] == 34))
panels$characteristic_associations <- ggplot(associations,
  aes(100 * conditional_adjusted_fraction, characteristic, color = response_label)) +
  geom_vline(xintercept = 0, color = "gray65", linewidth = 0.3) + geom_point(size = 2.8) +
  geom_text(aes(label = paste0("n = ", sites)), nudge_y = 0.2, size = 3, show.legend = FALSE) +
  facet_wrap(~response_label, nrow = 1) +
  scale_color_manual(values = c("Microbial composition" = "#0072B2", "Molecular composition" = "#D55E00"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  labs(x = "Adjusted variation beyond conditions (%)", y = NULL) + plot_theme
properties <- tables$paired_properties %>% mutate(model_label = factor(model,
  levels = rev(c("reference", "sediment_adjusted", "hydraulic_subset")),
  labels = rev(c("Drainage + detection", "+ Organic content + enzymes", "+ Conductivity (subset)"))))
stopifnot(identical(properties$sites, c(33, 33, 32)))
panels$paired_property_adjustment <- ggplot(properties,
  aes(100 * conditional_adjusted_fraction, model_label)) +
  geom_vline(xintercept = 0, color = "gray65", linewidth = 0.3) + geom_point(color = "#0072B2", size = 2.8) +
  geom_text(aes(label = paste0("n = ", sites)), nudge_y = 0.2, size = 3) +
  scale_x_continuous(limits = c(0, 5.2)) +
  labs(x = "Adjusted microbial variation associated\nwith molecular-property PCs (%)", y = NULL) + plot_theme
panels$characteristics <- wrap_plots(panels$characteristic_associations,
  panels$paired_property_adjustment, widths = c(1.7, 1), nrow = 1) + plot_annotation(tag_levels = "A")
panels$microbial_pca_organic <- ggplot(microbial, aes(pc1, pc2, color = sediment_organic_content_percent)) +
  geom_point(size = 2.8) + scale_color_viridis_c(name = "Sediment organic content (%)") + coord_equal() +
  labs(x = sprintf("Microbial Hellinger PC1 (%.1f%%)", pca_percent[1]),
    y = sprintf("Microbial Hellinger PC2 (%.1f%%)", pca_percent[2])) + plot_theme

# 7. Refit descriptive UniFrac coordinates on all 34 sediment samples only --------
# Reuse saved pairwise distances; no tree rebuilding or new hypothesis tests.
unifrac_rows <- list()
unifrac_summary <- list()
for (metric in c("weighted", "unweighted")) {
  saved <- tables[[metric]]
  distance_matrix <- as.matrix(saved[, metadata$sample_id])
  rownames(distance_matrix) <- saved$sample_id
  sediment_distance <- distance_matrix[metadata$sample_id, , drop = FALSE]
  stopifnot(all(is.finite(sediment_distance)), all(diag(sediment_distance) == 0),
    max(abs(sediment_distance - t(sediment_distance))) < 1e-12)
  fit <- wcmdscale(as.dist(sediment_distance), eig = TRUE, add = "lingoes")
  corrected <- sqrt(sediment_distance^2 + 2 * fit$ac)
  diag(corrected) <- 0
  reconstruction_error <- max(abs(as.matrix(dist(fit$points)) - corrected))
  stopifnot(reconstruction_error < 1e-8)
  percent <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
  scores <- metadata %>% select(sample_id, site_code) %>%
    mutate(axis1 = fit$points[, 1], axis2 = fit$points[, 2], metric) %>%
    left_join(site %>% select(site_code, network_group), by = "site_code")
  metric_label <- paste0(toupper(substr(metric, 1, 1)), substring(metric, 2), " UniFrac")
  panels[[paste0("unifrac_", metric)]] <- ggplot(scores, aes(axis1, axis2, color = network_group)) +
    geom_point(size = 2.8) + scale_color_manual(values = group_colors, name = NULL) + coord_equal() +
    labs(x = sprintf("%s PCo1 (%.1f%%)", metric_label, percent[1]),
      y = sprintf("PCo2 (%.1f%%)", percent[2])) + plot_theme
  unifrac_rows[[metric]] <- scores
  unifrac_summary[[metric]] <- tibble(metric, samples = nrow(scores), lingoes_constant = fit$ac,
    axis1_percent = percent[1], axis2_percent = percent[2], reconstruction_error)
}
panels$unifrac <- (wrap_plots(panels$unifrac_weighted, panels$unifrac_unweighted,
  nrow = 1, guides = "collect") + plot_annotation(tag_levels = "A")) &
  theme(legend.position = "bottom")
equal_scores <- tables$equal_scores %>% filter(view == "equal_feature") %>%
  select(site_code, axis1, axis2) %>% left_join(site %>% select(site_code, network_group), by = "site_code")
equal_percent <- tables$equal_ordination %>% filter(view == "equal_feature")
stopifnot(nrow(equal_scores) == 44)
panels$molecular_pcoa_equal <- ggplot(equal_scores, aes(axis1, axis2, color = network_group)) +
  geom_point(size = 2.8) + scale_color_manual(values = group_colors, name = NULL) + coord_equal() +
  labs(x = sprintf("Equal-feature molecular PCo1 (%.1f%%)", equal_percent$axis1_percent),
    y = sprintf("PCo2 (%.1f%%)", equal_percent$axis2_percent)) + plot_theme

# 8. Save every panel and composite, then verify source preservation --------------
stopifnot(setequal(names(panels), panel_design$key))
panel_design <- panel_design %>% mutate(file = paste0(figure_prefix, key, ".pdf"))
for (j in seq_len(nrow(panel_design))) {
  output_dir <- figures_dir
  if (panel_design$key[j] == "dispersion_paired") {
    output_dir <- manuscript_dir
  }
  ggsave(file.path(output_dir, panel_design$file[j]), panels[[panel_design$key[j]]],
    width = panel_design$width[j], height = panel_design$height[j], bg = "white")
}
write_csv(panel_design, file.path(results_dir, "figure_manifest.csv"))
write_csv(site, file.path(results_dir, "sampling_map_sites.csv"))
write_csv(microbial, file.path(results_dir, "microbial_composition_plot_data.csv"))
write_csv(molecular, file.path(results_dir, "molecular_composition_plot_data.csv"))
write_csv(class_summary, file.path(results_dir, "molecular_class_plot_data.csv"))
write_csv(bind_rows(dispersion_rows), file.path(results_dir, "dispersion_plot_data.csv"))
write_csv(associations, file.path(results_dir, "characteristic_association_plot_data.csv"))
write_csv(properties, file.path(results_dir, "paired_property_plot_data.csv"))
write_csv(bind_rows(unifrac_rows), file.path(results_dir, "sediment_unifrac_scores.csv"))
write_csv(bind_rows(unifrac_summary), file.path(results_dir, "sediment_unifrac_summary.csv"))
checks$final_md5 <- unname(tools::md5sum(checks$path))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "preserved_source_checksums.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
message("Saved ", nrow(panel_design), " sediment-focused PDFs; all original figures and inputs unchanged.")
