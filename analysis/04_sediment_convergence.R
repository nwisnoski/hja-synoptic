# Figure 3: sediment dispersion and molecular classes on identical sites.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "matched_sediment_dispersion")
class_results_dir <- here("results", "fticr_2016", "river_composition")
figures_dir <- here("figures", "supplement")
manuscript_dir <- here("figures", "manuscript")
target_depth <- 10000
target_features <- 391
n_feature_draws <- 1000
n_permutations <- 9999
n_balanced_draws <- 1000
permutation_seed <- 20261015
resampling_seed <- 20261016
plot_seed <- 20261017
point_size <- 3
point_alpha <- 0.5
panel_aspect <- 0.8
panel_width <- 5.8
figure_height <- 5.2
composite_height <- 11.5
group_levels <- c("Headwater", "Intermediate", "Mainstem")
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
contrast_names <- c("Headwater-Intermediate", "Headwater-Mainstem", "Intermediate-Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
test_design <- tribble(
  ~metric, ~role, ~label,
  "molecular_equal_features", "primary", "Molecular signatures (equal-feature Jaccard)",
  "microbial_hellinger", "primary", "Sediment microbes (Hellinger)",
  "molecular_jaccard", "sensitivity", "Molecular signatures (complete Jaccard)",
  "microbial_jaccard", "sensitivity", "Sediment microbes (binary Jaccard)") %>%
  mutate(seed = permutation_seed)
draw_design <- tibble(run = seq_len(n_balanced_draws), seed = resampling_seed + run)
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "none")

site_panel_aspect <- 0.8
site_figure_width <- 7.5
site_figure_height <- 5.8
class_order <- c("Lipid", "Protein", "AminoSugar", "Carb", "Lignin", "Tannin", "ConHC", "UnsatHC", "Other")
class_labels <- c("Lipid-like", "Protein-like", "Amino-sugar-like", "Carbohydrate-like",
  "Lignin-like", "Tannin-like", "Condensed-hydrocarbon-like", "Unsaturated-hydrocarbon-like", "Other")
class_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#882255", "#AA4499")
names(class_colors) <- class_labels
class_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 12),
    legend.position = "bottom"
  )

# 2. Audit the intersection without translating or renaming sample identities -----
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "network_dispersion", "site_network_position_audit.csv"),
  here("results", "fticr_2016", "river_composition", "mean_subsampled_jaccard.csv"),
  here("results", "fticr_2016", "river_composition", "subsampling_summary.csv"))
preserved_paths <- c(input_paths,
  list.files(here("results", "fticr_2016", "network_dispersion"), full.names = TRUE),
  here("figures", "supplement", "2016_fticr_dispersion_equal_features.pdf"))
preserved_paths <- unique(preserved_paths[file.exists(preserved_paths)])
checks <- tibble(path = preserved_paths, initial_md5 = unname(tools::md5sum(preserved_paths)))
metadata <- read_csv(input_paths[1], col_types = cols(site_code = col_character())) %>%
  filter(included_10k, habitat == "sediment") %>%
  select(sample_id, site_code, site_stream_order, site_segment_number, site_utm_x_m, site_utm_y_m)
molecular <- read_csv(input_paths[3], col_types = cols(site_code = col_character()))
network <- read_csv(input_paths[4], col_types = cols(site_code = col_character()))
stopifnot(!anyDuplicated(metadata$site_code), !anyDuplicated(metadata$sample_id),
  !anyDuplicated(molecular$site_code), !anyDuplicated(network$site_code),
  setequal(molecular$site_code, network$site_code))
overlap <- full_join(metadata %>% transmute(site_code, sample_id, microbial_10k = TRUE),
  molecular %>% transmute(site_code, molecular_profile = TRUE), by = "site_code") %>%
  mutate(microbial_10k = replace_na(microbial_10k, FALSE),
    molecular_profile = replace_na(molecular_profile, FALSE),
    included = microbial_10k & molecular_profile,
    reason = case_when(included ~ "matched sediment microbial and molecular site",
      !molecular_profile ~ "no molecular profile", !microbial_10k ~ "no sediment sample retained at 10K"))
write_csv(overlap, file.path(results_dir, "site_inclusion_audit.csv"))
site <- metadata %>% filter(site_code %in% molecular$site_code) %>%
  arrange(site_code) %>% mutate(network_group = factor(case_when(
    site_stream_order %in% headwater_orders ~ "Headwater",
    site_stream_order %in% intermediate_orders ~ "Intermediate",
    site_stream_order == mainstem_order ~ "Mainstem"), levels = group_levels))
reference <- network[match(site$site_code, network$site_code), ]
stopifnot(!anyNA(site), identical(reference$site_code, site$site_code),
  all(as.character(site$network_group) == reference$network_group),
  all(site$site_stream_order == reference$site_stream_order),
  all(site$site_segment_number == reference$site_segment_number),
  max(abs(site$site_utm_x_m - reference$site_utm_x_m)) < 1e-8,
  max(abs(site$site_utm_y_m - reference$site_utm_y_m)) < 1e-8)
site$detected_molecular_features <- reference$detected_primary_features
site$basin_outlet_stream_distance_m <- reference$basin_outlet_stream_distance_m
write_csv(site, file.path(results_dir, "matched_site_manifest.csv"))

# 3. Retain preparation settings and construct distances on the shared sites -------
count_table <- read_csv(input_paths[2], show_col_types = FALSE)
stopifnot(!anyDuplicated(count_table$asv_id), all(site$sample_id %in% names(count_table)))
community <- t(as.matrix(count_table[, site$sample_id]))
colnames(community) <- count_table$asv_id
rownames(community) <- site$site_code
community <- community[, colSums(community) > 0, drop = FALSE]
pa <- as.matrix(molecular[match(site$site_code, molecular$site_code), -1])
rownames(pa) <- site$site_code
pa <- pa[, colSums(pa) > 0, drop = FALSE]
stopifnot(all(is.finite(community)), all(community >= 0), all(community == floor(community)),
  all(rowSums(community) == target_depth), all(pa %in% c(0, 1)),
  all(rowSums(pa) == site$detected_molecular_features), min(rowSums(pa)) == target_features)
feature_summary <- read_csv(input_paths[6], show_col_types = FALSE)
stopifnot(feature_summary$target_features == target_features,
  feature_summary$completed_runs == n_feature_draws)
thinned <- read_csv(input_paths[5], col_types = cols(site_code = col_character()))
stopifnot(!anyDuplicated(thinned$site_code), all(site$site_code %in% thinned$site_code))
# Individual-site thinning does not change when other sites are excluded.
# Reuse the completed mean pairwise distances at the unchanged 391-feature target.
equal_features <- as.matrix(thinned[match(site$site_code, thinned$site_code), site$site_code])
rownames(equal_features) <- site$site_code
distances <- list(molecular_equal_features = equal_features,
  microbial_hellinger = as.matrix(dist(decostand(community, method = "hellinger"))),
  molecular_jaccard = as.matrix(vegdist(pa, method = "jaccard", binary = TRUE)),
  microbial_jaccard = as.matrix(vegdist(community, method = "jaccard", binary = TRUE)))

# 4. Refit spatial medians and use the same residual permutations for each view ----
set.seed(permutation_seed)
permutations <- shuffleSet(nrow(site), nset = n_permutations)
dispersion_rows <- list()
test_rows <- list()
contrast_rows <- list()
correction_rows <- list()
corrected_distances <- list()
for (test in seq_len(nrow(test_design))) {
  metric <- test_design$metric[test]
  raw <- distances[[metric]]
  stopifnot(identical(rownames(raw), site$site_code), identical(colnames(raw), site$site_code),
    all(is.finite(raw)), max(abs(raw - t(raw))) < 1e-12, all(diag(raw) == 0))
  ordination <- wcmdscale(as.dist(raw), eig = TRUE, add = "lingoes")
  correction <- ordination$ac
  if (is.null(correction)) {
    correction <- 0
  }
  corrected <- sqrt(raw^2 + 2 * correction)
  diag(corrected) <- 0
  corrected_distances[[metric]] <- corrected
  write_csv(tibble(site_code = site$site_code, as_tibble(raw, .name_repair = "minimal")),
    file.path(results_dir, paste0(metric, "_distance_matrix.csv")))
  correction_rows[[test]] <- tibble(metric, lingoes_constant = correction,
    max_distance_change = max(abs(corrected - raw)))
  fit <- betadisper(as.dist(corrected), site$network_group, type = "median", bias.adjust = TRUE)
  tested <- permutest(fit, permutations = permutations, pairwise = TRUE, parallel = 1)
  means <- tapply(fit$distances, site$network_group, mean)
  dispersion_rows[[test]] <- site %>% mutate(metric, distance = fit$distances)
  test_rows[[test]] <- tibble(metric, role = test_design$role[test], sites = nrow(site),
    headwater_mean = means["Headwater"], intermediate_mean = means["Intermediate"],
    mainstem_mean = means["Mainstem"], headwater_mainstem_ratio = means["Headwater"] / means["Mainstem"],
    statistic = tested$tab[1, "F"], p_value = tested$tab[1, "Pr(>F)"])
  for (contrast in contrast_names) {
    groups <- strsplit(contrast, "-", fixed = TRUE)[[1]]
    contrast_rows[[length(contrast_rows) + 1]] <- tibble(metric, role = test_design$role[test], contrast,
      difference = means[groups[1]] - means[groups[2]], ratio = means[groups[1]] / means[groups[2]],
      p_value = tested$pairwise$permuted[contrast])
  }
}
dispersion <- bind_rows(dispersion_rows)
tests <- bind_rows(test_rows) %>% group_by(role) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
contrasts <- bind_rows(contrast_rows) %>% group_by(role) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
stopifnot(all(is.finite(dispersion$distance)), !anyNA(tests$p_value), !anyNA(contrasts$p_value))
write_csv(dispersion, file.path(results_dir, "site_dispersion.csv"))
write_csv(tests, file.path(results_dir, "three_group_tests.csv"))
write_csv(contrasts, file.path(results_dir, "pairwise_group_contrasts.csv"))
write_csv(bind_rows(correction_rows), file.path(results_dir, "distance_correction_audit.csv"))

# 5. Select identical segment-balanced site sets for molecules and microbes -------
mainstem_indices <- which(site$network_group == "Mainstem")
balanced_n <- length(mainstem_indices)
stopifnot(n_distinct(site$site_segment_number[mainstem_indices]) == balanced_n, balanced_n >= 3)
balanced_rows <- list()
selection_rows <- list()
start_time <- proc.time()[["elapsed"]]
for (run in seq_len(n_balanced_draws)) {
  set.seed(draw_design$seed[run])
  selected <- integer(0)
  for (group in group_levels) {
    candidates <- which(site$network_group == group)
    segments <- unique(site$site_segment_number[candidates])
    stopifnot(length(segments) >= balanced_n)
    if (group == "Mainstem") {
      selected <- c(selected, mainstem_indices)
    } else {
      selected_segments <- segments[sample.int(length(segments), balanced_n)]
      for (segment in selected_segments) {
        indices <- candidates[site$site_segment_number[candidates] == segment]
        selected <- c(selected, indices[sample.int(length(indices), 1)])
      }
    }
  }
  labels <- site$network_group[selected]
  selection_rows[[run]] <- site[selected, ] %>%
    transmute(run, seed = draw_design$seed[run], site_code, sample_id, network_group,
      segment = site_segment_number)
  for (metric in test_design$metric) {
    diagnostic <- character(0)
    raw <- corrected_distances[[metric]][selected, selected]
    fit <- tryCatch(withCallingHandlers(
      betadisper(as.dist(raw), labels, type = "median", bias.adjust = TRUE),
      warning = function(w) {
        diagnostic <<- c(diagnostic, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    succeeded <- !inherits(fit, "error")
    means <- setNames(rep(NA_real_, length(group_levels)), group_levels)
    if (succeeded) {
      means <- tapply(fit$distances, labels, mean)
      succeeded <- all(is.finite(means)) && all(means > 0)
    } else {
      diagnostic <- c(diagnostic, conditionMessage(fit))
    }
    balanced_rows[[length(balanced_rows) + 1]] <- tibble(run, seed = draw_design$seed[run], metric,
      samples_per_group = balanced_n, succeeded, diagnostic = paste(diagnostic, collapse = " | "),
      headwater_mean = means["Headwater"], intermediate_mean = means["Intermediate"],
      mainstem_mean = means["Mainstem"], headwater_mainstem_ratio = means["Headwater"] / means["Mainstem"],
      intermediate_mainstem_ratio = means["Intermediate"] / means["Mainstem"],
      headwater_intermediate_ratio = means["Headwater"] / means["Intermediate"])
  }
  if (run == 1) {
    message("Representative four-view balanced draw: ", round(proc.time()[["elapsed"]] - start_time, 3), " s")
  }
}
balanced_runtime_seconds <- proc.time()[["elapsed"]] - start_time
balanced <- bind_rows(balanced_rows)
selections <- bind_rows(selection_rows)
selection_audit <- selections %>% group_by(run, network_group) %>%
  summarize(samples = n(), segments = n_distinct(segment), .groups = "drop")
stopifnot(nrow(balanced) == n_balanced_draws * nrow(test_design),
  all(selection_audit$samples == balanced_n), all(selection_audit$segments == balanced_n))
write_csv(balanced, file.path(results_dir, "balanced_segment_draws.csv"))
write_csv(selections, file.path(results_dir, "balanced_site_selection.csv"))
balanced_summary <- balanced %>% group_by(metric) %>%
  summarize(draws = n(), successful_draws = sum(succeeded), failed_draws = sum(!succeeded),
    draws_with_diagnostics = sum(nchar(diagnostic) > 0),
    median_ratio = median(headwater_mainstem_ratio[succeeded]),
    min_ratio = min(headwater_mainstem_ratio[succeeded]), max_ratio = max(headwater_mainstem_ratio[succeeded]),
    fraction_headwater_greater_mainstem = mean(headwater_mainstem_ratio[succeeded] > 1),
    fraction_intermediate_greater_mainstem = mean(intermediate_mainstem_ratio[succeeded] > 1),
    fraction_headwater_greater_intermediate = mean(headwater_intermediate_ratio[succeeded] > 1), .groups = "drop")
write_csv(balanced_summary, file.path(results_dir, "balanced_segment_summary.csv"))

# 6. Save primary panels and their combination without altering full-site figures --
panels <- list()
for (test in which(test_design$role == "primary")) {
  metric <- test_design$metric[test]
  panel <- dispersion %>% filter(.data$metric == .env$metric)
  sizes <- table(panel$network_group)
  x_labels <- paste0(group_levels, "\n(n = ", as.integer(sizes[group_levels]), ")")
  set.seed(plot_seed)
  p <- ggplot(panel, aes(network_group, distance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.5) +
    geom_jitter(width = 0.1, height = 0, size = point_size, alpha = point_alpha) +
    facet_wrap(~metric, labeller = as_labeller(setNames(test_design$label, test_design$metric))) +
    scale_color_manual(values = group_colors) + scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = "Distance to group spatial median") +
    plot_theme + theme(aspect.ratio = panel_aspect)
  panels[[metric]] <- p
  ggsave(file.path(figures_dir, paste0("2016_matched_sediment_dispersion_", metric, ".pdf")),
    p, width = panel_width, height = figure_height, bg = "white")
}
ggsave(file.path(figures_dir, "2016_matched_sediment_dispersion_comparison.pdf"),
  wrap_plots(panels, ncol = 2) + plot_annotation(tag_levels = "A"),
  width = 2 * panel_width, height = figure_height, bg = "white")
# 7. Align the existing class fractions to the audited 44-site network groups.
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
class_input_paths <- c(
  file.path(class_results_dir, "site_class_fractions.csv"),
  here("results", "fticr_2016", "network_dispersion", "site_network_position_audit.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "tables", "primary_peak_properties.csv"),
  file.path(figures_dir, "2016_fticr_site_class_composition.pdf")
)
class_input_hashes <- unname(tools::md5sum(class_input_paths))
stopifnot(!anyNA(class_input_hashes))
class_table <- read_csv(class_input_paths[1], col_types = cols(site_code = col_character()))
class_network <- read_csv(class_input_paths[2], col_types = cols(site_code = col_character()))
class_incidence <- read_csv(class_input_paths[3], col_types = cols(site_code = col_character()))
class_properties <- read_csv(class_input_paths[4], show_col_types = FALSE)
stopifnot(
  !anyDuplicated(class_network$site_code),
  !anyDuplicated(class_incidence$site_code),
  !anyDuplicated(class_properties$peak_id),
  !anyDuplicated(class_table[c("site_code", "source_class")]),
  setequal(class_table$site_code, class_network$site_code),
  setequal(class_incidence$site_code, class_network$site_code),
  setequal(class_table$source_class, class_order)
)
class_network$expected_group <- case_when(
  class_network$site_stream_order %in% headwater_orders ~ "Headwater",
  class_network$site_stream_order %in% intermediate_orders ~ "Intermediate",
  class_network$site_stream_order == mainstem_order ~ "Mainstem"
)
stopifnot(!anyNA(class_network$expected_group), all(class_network$network_group == class_network$expected_group))
class_site_order <- class_network %>%
  mutate(network_group = factor(network_group, levels = group_levels)) %>%
  arrange(network_group, site_drainage_area_ha, site_code)
class_reference <- class_network[match(class_table$site_code, class_network$site_code), ]
stopifnot(
  all(class_table$site_stream_order == class_reference$site_stream_order),
  all(abs(class_table$site_drainage_area_ha - class_reference$site_drainage_area_ha) < 1e-8),
  all(is.finite(class_table$fraction_detected)),
  all(class_table$fraction_detected >= 0 & class_table$fraction_detected <= 1)
)
class_table$network_group <- factor(class_reference$network_group, levels = group_levels)
class_table$class_label <- factor(class_table$source_class, levels = class_order, labels = class_labels)
class_table$site_code <- factor(class_table$site_code, levels = class_site_order$site_code)

# 8. Verify that the saved fractions agree with current binary detections.
# Keep workbook classifications and detection rules unchanged.
class_pa <- as.matrix(class_incidence[, -1])
class_properties <- class_properties[match(colnames(class_pa), class_properties$peak_id), ]
stopifnot(
  identical(colnames(class_pa), class_properties$peak_id),
  all(class_pa %in% c(0, 1)),
  all(rowSums(class_pa) > 0),
  all(class_properties$molecular_class_source %in% class_order)
)
for (source_name in class_order) {
  selected <- class_properties$molecular_class_source == source_name
  expected <- rowSums(class_pa[, selected, drop = FALSE]) / rowSums(class_pa)
  current <- class_table %>% filter(source_class == source_name)
  expected <- expected[match(as.character(current$site_code), class_incidence$site_code)]
  stopifnot(nrow(current) == nrow(class_network), all(abs(current$fraction_detected - expected) < 1e-12))
}
site_sums <- class_table %>% group_by(site_code) %>%
  summarize(fraction_sum = sum(fraction_detected), .groups = "drop")
stopifnot(all(abs(site_sums$fraction_sum - 1) < 1e-12))

# 9. Average site fractions within groups, giving every site equal weight.
# Pooling raw feature counts would give feature-rich sites greater weight.
class_group_summary <- class_table %>% group_by(network_group, source_class, class_label) %>%
  summarize(
    n_sites = n(),
    mean_fraction = mean(fraction_detected),
    sd_fraction = sd(fraction_detected),
    .groups = "drop"
  )
group_sums <- class_group_summary %>% group_by(network_group) %>%
  summarize(fraction_sum = sum(mean_fraction), .groups = "drop")
stopifnot(all(abs(group_sums$fraction_sum - 1) < 1e-12))
class_group_counts <- table(factor(class_network$network_group, levels = group_levels))
class_group_labels <- paste0(group_levels, "\n(n = ", as.integer(class_group_counts), ")")
names(class_group_labels) <- group_levels

# 10. Retain full-site class panels and match the manuscript to dispersion.
group_plot <- ggplot(class_group_summary, aes(network_group, mean_fraction, fill = class_label)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_x_discrete(labels = class_group_labels) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Mean fraction of detected features", fill = NULL) +
  class_theme +
  guides(fill = guide_legend(nrow = 3))
ggsave(file.path(figures_dir, "2016_fticr_network_group_class_means.pdf"),
  group_plot, width = 9, height = 5.5, bg = "white")

# Retain the original site bars, ordered by drainage area inside each group.
# Bottom facet labels serve as the group labels on the x axis.
# Facet widths follow site counts so every site's bar has the same width.
site_panels <- list()
for (group_name in c("All", "Matched", group_levels)) {
  current <- class_table
  if (group_name == "Matched") {
    current <- class_table %>% filter(as.character(site_code) %in% site$site_code)
    stopifnot(setequal(as.character(current$site_code), site$site_code),
      n_distinct(current$site_code) == 33)
    matched_reference <- site[match(as.character(current$site_code), site$site_code), ]
    stopifnot(all(as.character(current$network_group) == as.character(matched_reference$network_group)))
  } else if (group_name != "All") {
    current <- class_table %>% filter(network_group == group_name)
  }
  current_sizes <- current %>% distinct(site_code, network_group) %>% count(network_group)
  current_labels <- setNames(paste0(current_sizes$network_group, "\n(n = ",
    current_sizes$n, ")"), as.character(current_sizes$network_group))
  site_panels[[group_name]] <- ggplot(current, aes(site_code, fraction_detected, fill = class_label)) +
    geom_col(width = 0.95) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25)) +
    facet_grid(. ~ network_group, scales = "free_x",
      space = if (group_name %in% c("All", "Matched")) "free_x" else "fixed", switch = "x",
      labeller = labeller(network_group = current_labels)) +
    labs(x = NULL, y = "Fraction of detected features", fill = NULL) +
    class_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    guides(fill = guide_legend(nrow = 3))
  if (group_name %in% c("All", "Matched")) {
    site_panels[[group_name]] <- site_panels[[group_name]] +
      guides(fill = guide_legend(ncol = 2, byrow = TRUE))
  } else {
    site_panels[[group_name]] <- site_panels[[group_name]] +
      theme(aspect.ratio = site_panel_aspect)
  }
  file_suffix <- "by_network_group"
  class_panel_width <- site_figure_width
  class_panel_height <- site_figure_height
  if (group_name == "Matched") {
    file_suffix <- "matched_sites"
    class_panel_width <- 2 * panel_width
  } else if (group_name != "All") {
    file_suffix <- tolower(group_name)
    class_panel_width <- 5.8
    class_panel_height <- 5.2
  }
  ggsave(file.path(figures_dir, paste0("2016_fticr_site_class_composition_", file_suffix, ".pdf")),
    site_panels[[group_name]], width = class_panel_width, height = class_panel_height, bg = "white")
}
# The lower panel describes broad class fractions, not feature-level turnover.
# All three panels use exactly the same 33 sites (18, 10, and 5 by network group).
combined <- wrap_plots(panels[["molecular_equal_features"]], panels[["microbial_hellinger"]],
  site_panels$Matched + guides(fill = guide_legend(nrow = 3, byrow = TRUE)),
  design = "AB\nCC", heights = c(1, 1.25)) + plot_annotation(tag_levels = "A")
ggsave(file.path(manuscript_dir, "2016_matched_sediment_dispersion_comparison.pdf"),
  combined, width = 2 * panel_width, height = composite_height, bg = "white")
stopifnot(identical(class_input_hashes, unname(tools::md5sum(class_input_paths))))

# 11. Confirm that the reused inputs and full-site dispersion results are intact.
checks$final_md5 <- unname(tools::md5sum(checks$path))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "preserved_input_output_checksums.csv"))
write_csv(tibble(n_permutations, n_balanced_draws, target_depth, target_features,
  n_feature_draws, permutation_seed, resampling_seed, plot_seed, balanced_runtime_seconds),
  file.path(results_dir, "run_settings.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
print(tests)
print(balanced_summary)
