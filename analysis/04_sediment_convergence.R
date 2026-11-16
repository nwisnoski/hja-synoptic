# Figure 3: downstream decline of sediment molecular richness, and the share of
# molecular and microbial dissimilarity that comes from richness differences
# along flow paths. Dispersion and molecular class views of the same sites are
# supporting panels.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
library(adespatial)
library(igraph)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "matched_sediment_dispersion")
class_results_dir <- here("results", "fticr_2016", "river_composition")
figures_dir <- here("figures", "supplement", "sediment_comparison")
manuscript_dir <- here("figures", "manuscript")
target_depth <- 10000
target_features <- 391
n_feature_draws <- 1000
n_permutations <- 9999
n_balanced_draws <- 1000
permutation_seed <- 20261015
significance_threshold <- 0.05
resampling_seed <- 20261016
plot_seed <- 20261017
point_size <- 3
point_alpha <- 0.5
panel_aspect <- 0.8
panel_width <- 5.8
figure_height <- 5.2
# The Figure 3 panels are built at printed size: each panel has a fixed size and
# the canvas is measured around it. Dispersion and class panels keep the earlier
# canvas-based sizes until the whole figure set is restyled together.
print_panel_width <- 2.6
print_panel_height <- 2.0
facet_panel_width <- 2.1
facet_panel_height <- 1.6
print_base_size <- 9
site_point_size <- 1.6
pair_point_size <- 0.8
pair_point_alpha <- 0.5
group_levels <- c("Headwater", "Intermediate", "Mainstem")
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
contrast_names <- c("Headwater-Intermediate", "Headwater-Mainstem", "Intermediate-Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
# Dispersion views. The primary views are the standard dissimilarity for each
# data type on the complete data. Replacement removes the part of binary
# dissimilarity that comes from richness differences. The equal-feature view is
# kept only to document that thinning binary detections restates richness.
test_design <- tribble(
  ~metric, ~role, ~label,
  "molecular_equal_features", "thinning_diagnostic", "Molecular signatures (equal-feature Jaccard)",
  "microbial_hellinger", "primary", "Sediment microbes (Hellinger)",
  "molecular_jaccard", "primary", "Molecular signatures (complete Jaccard)",
  "microbial_jaccard", "sensitivity", "Sediment microbes (binary Jaccard)",
  "molecular_replacement", "replacement", "Molecular signatures (Jaccard replacement)",
  "microbial_replacement", "replacement", "Sediment microbes (Jaccard replacement)") %>%
  mutate(seed = permutation_seed)
# Each view is compared with the richness of its own data type.
molecular_metrics <- c("molecular_equal_features", "molecular_jaccard", "molecular_replacement")
# Site-level models of the downstream decline. All variables are log10, so
# slopes are scaling exponents. The joint model asks whether drainage area
# explains feature richness beyond sediment organic content.
model_design <- tribble(
  ~model, ~role, ~response, ~predictors,
  "richness_drainage", "simple", "log10_features", "log10_drainage",
  "organic_drainage", "simple", "log10_organic", "log10_drainage",
  "richness_organic", "simple", "log10_features", "log10_organic",
  "richness_joint", "joint", "log10_features", "log10_organic + log10_drainage",
  "asv_richness_drainage", "microbial_contrast", "log10_asvs", "log10_drainage")
# Binary Jaccard dissimilarity between two sites is the sum of replacement and
# richness difference (Podani and Schmera 2011, Legendre 2014). The richness-
# difference share is richness difference divided by Jaccard dissimilarity.
# Jaccard dissimilarity itself is tested in supporting/10_distance_decay.R.
component_design <- tribble(
  ~component, ~label,
  "replacement", "Replacement",
  "richness_difference", "Richness difference",
  "richness_share", "Richness-difference share")
data_type_labels <- c(molecular = "FT-ICR features", microbial = "Microbial ASVs")
# Two separations between sites: channel distance and the ratio of their
# drainage areas. Each is tested within flow-connected and unconnected pairs.
separation_design <- tribble(
  ~separation, ~label,
  "watercourse_km", "Watercourse distance (km)",
  "log10_drainage_ratio", "Drainage area ratio (log10)")
connection_levels <- c("Flow-connected", "Not flow-connected")
connection_colors <- c("Flow-connected" = "#000000", "Not flow-connected" = "#BBBBBB")
# Fitted lines for unconnected pairs use a darker gray so they stay visible over
# the light gray points.
connection_line_colors <- c("Flow-connected" = "#000000", "Not flow-connected" = "#707070")
# Tolerance for matching graph distances to the centerline's own outlet distances.
distance_tolerance_m <- 0.01
draw_design <- tibble(run = seq_len(n_balanced_draws), seed = resampling_seed + run)
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "none",
    axis.text = element_text(size = 11), strip.text = element_text(size = 12),
    plot.tag = element_text(size = 14, face = "bold"))
print_theme <- theme_bw(base_size = print_base_size) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "none",
    panel.widths = grid::unit(print_panel_width, "in"),
    panel.heights = grid::unit(print_panel_height, "in"),
    plot.tag = element_text(size = 11, face = "bold"))
# Save a figure whose panels have fixed sizes on a canvas measured around them.
# Measuring text needs an open graphics device, so a null device is opened.
save_at_panel_size <- function(plot, path) {
  pdf(NULL)
  if (inherits(plot, "patchwork")) {
    grob <- patchworkGrob(plot)
  } else {
    grob <- ggplotGrob(plot)
  }
  width <- grid::convertWidth(sum(grob$widths), "in", valueOnly = TRUE)
  height <- grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE)
  dev.off()
  ggsave(path, plot, width = width, height = height, bg = "white")
}

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
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = grid::unit(0.45, "cm"),
    plot.tag = element_text(size = 14, face = "bold")
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
  here("results", "fticr_2016", "river_composition", "subsampling_summary.csv"),
  here("data", "derived", "analysis_inputs", "environment", "sediment_44_environment.csv"))
preserved_paths <- c(input_paths,
  list.files(here("results", "fticr_2016", "network_dispersion"), full.names = TRUE),
  here("figures", "supplement", "molecular", "2016_fticr_dispersion_equal_features.pdf"))
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
# Drainage area and sediment organic content for all 44 molecular sites.
environment <- read_csv(input_paths[7], col_types = cols(site_code = col_character())) %>%
  select(site_code, site_drainage_area_ha, sediment_organic_content_percent)
stopifnot(!anyDuplicated(environment$site_code), setequal(environment$site_code, network$site_code),
  all(is.finite(environment$site_drainage_area_ha)), all(environment$site_drainage_area_ha > 0),
  all(is.finite(environment$sediment_organic_content_percent)),
  all(environment$sediment_organic_content_percent > 0))

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
# Binary Jaccard dissimilarity divides into replacement, where each site holds
# features the other lacks, and richness difference, where one site holds more
# features than the other (Podani and Schmera 2011; Carvalho et al. 2012).
# Both data types use presence/absence on the complete data, with no thinning.
# The abundance-weighted form would instead measure differences in total
# abundance, which are zero for samples rarefied to equal depth.
community_pa <- (community > 0) * 1
site$detected_asvs <- as.integer(rowSums(community_pa))
partitions <- list(molecular = beta.div.comp(pa, coef = "J", quant = FALSE),
  microbial = beta.div.comp(community_pa, coef = "J", quant = FALSE))
distances <- list(molecular_equal_features = equal_features,
  microbial_hellinger = as.matrix(dist(decostand(community, method = "hellinger"))),
  molecular_jaccard = as.matrix(vegdist(pa, method = "jaccard", binary = TRUE)),
  microbial_jaccard = as.matrix(vegdist(community, method = "jaccard", binary = TRUE)),
  molecular_replacement = as.matrix(partitions$molecular$repl),
  microbial_replacement = as.matrix(partitions$microbial$repl))
for (data_type in names(partitions)) {
  dimnames(distances[[paste0(data_type, "_replacement")]]) <- list(site$site_code, site$site_code)
}
stopifnot(identical(names(distances), test_design$metric),
  max(abs(as.matrix(partitions$molecular$D) - distances$molecular_jaccard)) < 1e-12,
  max(abs(as.matrix(partitions$microbial$D) - distances$microbial_jaccard)) < 1e-12)

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
    message("Representative balanced draw across all views: ", round(proc.time()[["elapsed"]] - start_time, 3), " s")
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

# 6. Relate each site's distance to its richness -----------------------------------
# A measure suited to comparing dispersion among groups should not rank sites by
# the number of features or ASVs detected in them. Correlations near one mean
# that the view mostly restates richness.
richness_rows <- list()
for (test in seq_len(nrow(test_design))) {
  metric <- test_design$metric[test]
  current <- dispersion %>% filter(.data$metric == .env$metric)
  richness_measure <- "detected_asvs"
  if (metric %in% molecular_metrics) {
    richness_measure <- "detected_molecular_features"
  }
  current$richness <- current[[richness_measure]]
  within_group <- current %>% group_by(network_group) %>%
    summarize(rho = cor(distance, richness, method = "spearman"), .groups = "drop")
  richness_rows[[test]] <- tibble(metric, role = test_design$role[test], richness_measure,
    spearman_all_sites = cor(current$distance, current$richness, method = "spearman"),
    spearman_headwater = within_group$rho[within_group$network_group == "Headwater"],
    spearman_intermediate = within_group$rho[within_group$network_group == "Intermediate"],
    spearman_mainstem = within_group$rho[within_group$network_group == "Mainstem"])
}
richness_dependence <- bind_rows(richness_rows)
write_csv(richness_dependence, file.path(results_dir, "richness_dependence.csv"))
richness_panel <- dispersion %>% filter(metric %in% molecular_metrics) %>%
  mutate(metric = factor(metric, levels = molecular_metrics))
richness_plot <- ggplot(richness_panel, aes(detected_molecular_features, distance, color = network_group)) +
  geom_point(size = site_point_size, alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y", ncol = 3,
    labeller = as_labeller(setNames(test_design$label, test_design$metric))) +
  scale_color_manual(values = group_colors) +
  labs(x = "Detected FT-ICR features", y = "Distance to group median", color = NULL) +
  print_theme + theme(legend.position = "bottom",
    panel.widths = grid::unit(rep(facet_panel_width, length(molecular_metrics)), "in"),
    panel.heights = grid::unit(facet_panel_height, "in"))
save_at_panel_size(richness_plot,
  file.path(figures_dir, "2016_matched_sediment_dispersion_richness_dependence.pdf"))

# 7. Model the downstream decline in molecular richness ----------------------------
# Three simple regressions describe the chain from network position to sediment
# organic content to detected features. The joint model then asks whether
# drainage area still matters once organic content is known. Extracts used a
# fixed sediment mass, so detected features are counted per unit sediment.
# The matched 33 sites are the primary set; all 44 molecular sites are a check.
# ASV richness at equal sequencing depth is the microbial contrast.
site_sets <- list(
  matched_33 = site %>% select(site_code, detected_molecular_features, detected_asvs),
  molecular_44 = network %>% transmute(site_code, detected_molecular_features = detected_primary_features,
    detected_asvs = NA_integer_))
model_rows <- list()
model_data <- list()
for (site_set in names(site_sets)) {
  current <- site_sets[[site_set]] %>% left_join(environment, by = "site_code") %>%
    mutate(log10_features = log10(detected_molecular_features), log10_asvs = log10(detected_asvs),
      log10_drainage = log10(site_drainage_area_ha), log10_organic = log10(sediment_organic_content_percent))
  stopifnot(!anyNA(current$log10_features), !anyNA(current$log10_drainage), !anyNA(current$log10_organic))
  model_data[[site_set]] <- current
  for (i in seq_len(nrow(model_design))) {
    if (model_design$response[i] == "log10_asvs" && site_set != "matched_33") {
      next
    }
    fit <- lm(as.formula(paste(model_design$response[i], "~", model_design$predictors[i])), data = current)
    coefficients <- summary(fit)$coefficients
    model_rows[[length(model_rows) + 1]] <- tibble(site_set, model = model_design$model[i],
      role = model_design$role[i], response = model_design$response[i], term = rownames(coefficients),
      estimate = coefficients[, "Estimate"], std_error = coefficients[, "Std. Error"],
      statistic = coefficients[, "t value"], p_value = coefficients[, "Pr(>|t|)"],
      sites = nrow(current), r_squared = summary(fit)$r.squared)
  }
}
# Intercepts are kept for drawing lines; only slopes enter the BH families.
richness_models <- bind_rows(model_rows) %>% group_by(site_set, role) %>%
  mutate(p_bh = if_else(term == "(Intercept)", NA_real_,
    p.adjust(if_else(term == "(Intercept)", NA_real_, p_value), method = "BH"))) %>% ungroup()
write_csv(richness_models, file.path(results_dir, "molecular_richness_models.csv"))

# 8. Measure the separation of every pair of sites along the stream network --------
# Each centerline row is one point; D is its channel distance to the basin outlet.
# The raw header has a trailing comma but data rows have seven fields.
geometry <- read_csv(here("data", "StreamCenterline.csv"), skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
geometry$node <- seq_len(nrow(geometry))
# Rows run down each segment, so consecutive rows of a segment are joined by the
# change in their outlet distance.
same_segment <- geometry$segment[-1] == geometry$segment[-nrow(geometry)]
channel_edges <- tibble(from = geometry$node[-nrow(geometry)][same_segment],
  to = geometry$node[-1][same_segment])
channel_edges$weight <- abs(geometry$outlet_distance[channel_edges$from] -
  geometry$outlet_distance[channel_edges$to])
# Segments meet at confluences, where their end points have identical coordinates.
shared_points <- split(geometry$node, paste(geometry$x, geometry$y))
shared_points <- shared_points[lengths(shared_points) > 1]
confluence_rows <- list()
for (i in seq_along(shared_points)) {
  nodes <- shared_points[[i]]
  confluence_rows[[i]] <- tibble(from = nodes[1], to = nodes[-1], weight = 0)
}
confluence_edges <- bind_rows(confluence_rows)
stream_graph <- graph_from_data_frame(bind_rows(channel_edges, confluence_edges),
  directed = FALSE, vertices = tibble(name = geometry$node))
# Check the graph against the source: along the connected network, the graph
# distance from every point to the outlet must equal the centerline's own D.
outlet_node <- which.min(geometry$outlet_distance)
graph_outlet_distance <- as.numeric(distances(stream_graph, v = outlet_node,
  weights = E(stream_graph)$weight))
membership <- components(stream_graph)$membership
reaches_outlet <- membership == membership[outlet_node]
graph_checks <- tibble(
  centerline_points = nrow(geometry),
  segments = n_distinct(geometry$segment),
  confluences = length(shared_points),
  components = max(membership),
  segments_not_reaching_outlet = paste(sort(unique(geometry$segment[!reaches_outlet])), collapse = ";"),
  points_not_reaching_outlet = sum(!reaches_outlet),
  max_outlet_distance_error_m = max(abs(graph_outlet_distance[reaches_outlet] -
    geometry$outlet_distance[reaches_outlet])))
stopifnot(graph_checks$max_outlet_distance_error_m < distance_tolerance_m)
# Each site is placed at the nearest centerline point of its own segment.
site$centerline_node <- NA_integer_
site$centerline_match_error_m <- NA_real_
for (i in seq_len(nrow(site))) {
  candidates <- which(geometry$segment == site$site_segment_number[i])
  stopifnot(length(candidates) > 0)
  offset <- sqrt((geometry$x[candidates] - site$site_utm_x_m[i])^2 +
    (geometry$y[candidates] - site$site_utm_y_m[i])^2)
  site$centerline_node[i] <- candidates[which.min(offset)]
  site$centerline_match_error_m[i] <- min(offset)
}
site_outlet_distance <- geometry$outlet_distance[site$centerline_node]
stopifnot(all(reaches_outlet[site$centerline_node]), all(site$centerline_match_error_m < 1),
  max(abs(site_outlet_distance - site$basin_outlet_stream_distance_m)) < distance_tolerance_m)
watercourse_m <- distances(stream_graph, v = site$centerline_node, to = site$centerline_node,
  weights = E(stream_graph)$weight)
# Two sites are flow-connected when water passes from one to the other, so the
# channel distance between them equals the difference of their outlet distances.
outlet_difference_m <- abs(outer(site_outlet_distance, site_outlet_distance, "-"))
connected_matrix <- abs(watercourse_m - outlet_difference_m) < distance_tolerance_m
straight_line_m <- as.matrix(dist(cbind(site$site_utm_x_m, site$site_utm_y_m)))
stopifnot(all(watercourse_m + distance_tolerance_m >= outlet_difference_m),
  all(watercourse_m + 1e-3 >= straight_line_m))
drainage <- model_data$matched_33$site_drainage_area_ha
stopifnot(identical(model_data$matched_33$site_code, site$site_code))
separations <- list(watercourse_km = watercourse_m / 1000,
  log10_drainage_ratio = abs(outer(log10(drainage), log10(drainage), "-")))
# One row per unordered pair of sites. For flow-connected pairs, the downstream
# site is the one nearer the outlet.
lower <- which(lower.tri(watercourse_m), arr.ind = TRUE)
ia <- lower[, "col"]
ib <- lower[, "row"]
pairs <- tibble(site_a = site$site_code[ia], site_b = site$site_code[ib],
  group_a = site$network_group[ia], group_b = site$network_group[ib],
  flow_connected = connected_matrix[lower],
  watercourse_km = separations$watercourse_km[lower], straight_line_km = straight_line_m[lower] / 1000,
  log10_drainage_ratio = separations$log10_drainage_ratio[lower],
  downstream_site = if_else(connected_matrix[lower],
    if_else(site_outlet_distance[ia] < site_outlet_distance[ib], site$site_code[ia], site$site_code[ib]),
    NA_character_))
write_csv(graph_checks %>% mutate(flow_connected_pairs = sum(pairs$flow_connected),
  site_pairs = nrow(pairs), max_site_match_error_m = max(site$centerline_match_error_m)),
  file.path(results_dir, "stream_graph_checks.csv"))

# 9. Divide Jaccard dissimilarity into replacement and richness difference ---------
# Total beta diversity and its two shares summarize each data type over all site
# pairs (Legendre 2014). For each pair, the richness-difference share says how
# much of the dissimilarity between the two sites comes from one holding more
# features or ASVs than the other.
share_rows <- list()
components <- list()
for (data_type in names(partitions)) {
  part <- partitions[[data_type]]$part
  share_rows[[data_type]] <- tibble(data_type, sites = nrow(site), site_pairs = nrow(pairs),
    total_beta_diversity = part[["BDtotal"]], replacement = part[["Repl"]],
    richness_difference = part[["RichDif"]], replacement_share = part[["Repl/BDtotal"]],
    richness_difference_share = part[["RichDif/BDtotal"]])
  jaccard <- as.matrix(partitions[[data_type]]$D)
  richness_difference <- as.matrix(partitions[[data_type]]$rich)
  stopifnot(all(jaccard[lower] > 0))
  richness_share <- richness_difference / jaccard
  diag(richness_share) <- 0
  components[[data_type]] <- list(replacement = as.matrix(partitions[[data_type]]$repl),
    richness_difference = richness_difference, richness_share = richness_share)
  pairs[[paste0(data_type, "_jaccard")]] <- jaccard[lower]
  for (component in component_design$component) {
    pairs[[paste(data_type, component, sep = "_")]] <- components[[data_type]][[component]][lower]
  }
}
partition_shares <- bind_rows(share_rows)
write_csv(partition_shares, file.path(results_dir, "jaccard_partition_shares.csv"))
write_csv(pairs, file.path(results_dir, "site_pair_components.csv"))

# 10. Test each component against separation within flow-connected and other pairs
# Water carries material between flow-connected sites, so composition may change
# along a flow path and not between separate branches. The correlation with each
# separation is computed within each set of pairs. Its null distribution permutes
# whole sites in the component matrix, with the permutations used for dispersion,
# while separations and flow connections stay in place. The last test compares
# the slopes of the two pair sets. All tests are two-sided.
connected <- pairs$flow_connected
flow_rows <- list()
for (data_type in names(partitions)) {
  for (component in component_design$component) {
    current <- components[[data_type]][[component]]
    observed <- current[lower]
    for (separation in separation_design$separation) {
      x <- pairs[[separation]]
      connected_line <- lm(observed[connected] ~ x[connected])
      unconnected_line <- lm(observed[!connected] ~ x[!connected])
      observed_r <- c(cor(observed[connected], x[connected]), cor(observed[!connected], x[!connected]))
      observed_difference <- coef(connected_line)[[2]] - coef(unconnected_line)[[2]]
      null_r <- matrix(NA_real_, n_permutations, 2)
      null_difference <- rep(NA_real_, n_permutations)
      for (k in seq_len(n_permutations)) {
        order <- permutations[k, ]
        shuffled <- current[order, order][lower]
        null_r[k, 1] <- cor(shuffled[connected], x[connected])
        null_r[k, 2] <- cor(shuffled[!connected], x[!connected])
        null_difference[k] <- cov(shuffled[connected], x[connected]) / var(x[connected]) -
          cov(shuffled[!connected], x[!connected]) / var(x[!connected])
      }
      flow_rows[[length(flow_rows) + 1]] <- tibble(data_type, component, separation,
        pair_set = c(connection_levels, "Slope difference"),
        site_pairs = c(sum(connected), sum(!connected), length(connected)),
        mean_value = c(mean(observed[connected]), mean(observed[!connected]), NA_real_),
        pearson_r = c(observed_r, NA_real_),
        intercept = c(coef(connected_line)[[1]], coef(unconnected_line)[[1]], NA_real_),
        slope = c(coef(connected_line)[[2]], coef(unconnected_line)[[2]], observed_difference),
        p_value = c((sum(abs(null_r[, 1]) >= abs(observed_r[1])) + 1) / (n_permutations + 1),
          (sum(abs(null_r[, 2]) >= abs(observed_r[2])) + 1) / (n_permutations + 1),
          (sum(abs(null_difference) >= abs(observed_difference)) + 1) / (n_permutations + 1)))
    }
  }
}
# Each separation forms two BH families: slopes within pair sets, and slope
# differences.
flow_tests <- bind_rows(flow_rows) %>%
  mutate(family = if_else(pair_set == "Slope difference", "slope_difference", "within_pair_set")) %>%
  group_by(separation, family) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup() %>% select(-family)
write_csv(flow_tests, file.path(results_dir, "flow_path_component_tests.csv"))

# 11. Draw the manuscript figure: organic content, richness, and flow paths -------
# Lines are drawn only for slopes whose adjusted test passes the threshold.
matched <- model_data$matched_33
simple <- richness_models %>% filter(site_set == "matched_33", role == "simple")
figure_panels <- list()
site_panel_design <- tribble(
  ~panel, ~x, ~y, ~x_label, ~y_label,
  "richness_drainage", "site_drainage_area_ha", "detected_molecular_features",
    "Drainage area (ha)", "Detected FT-ICR features",
  "organic_drainage", "site_drainage_area_ha", "sediment_organic_content_percent",
    "Drainage area (ha)", "Sediment organic content (%)",
  "richness_organic", "sediment_organic_content_percent", "detected_molecular_features",
    "Sediment organic content (%)", "Detected FT-ICR features")
for (i in seq_len(nrow(site_panel_design))) {
  design <- site_panel_design[i, ]
  fitted <- simple %>% filter(model == design$panel)
  slope <- fitted %>% filter(term != "(Intercept)")
  p <- ggplot(matched, aes(.data[[design$x]], .data[[design$y]])) +
    geom_point(size = site_point_size, alpha = 0.7) +
    labs(x = design$x_label, y = design$y_label) + print_theme
  # Drainage area spans three orders of magnitude and is labeled in powers of
  # ten. Feature count and organic content span less than one, so their log
  # axes carry plain numbers.
  if (design$x == "site_drainage_area_ha") {
    p <- p + scale_x_log10(limits = c(10, 10^4), breaks = 10^(1:4), labels = scales::label_log())
  } else {
    p <- p + scale_x_log10(breaks = c(3, 5, 10, 15))
  }
  if (design$y == "detected_molecular_features") {
    p <- p + scale_y_log10(breaks = c(500, 1000, 2000))
  } else {
    p <- p + scale_y_log10(breaks = c(3, 5, 10, 15))
  }
  if (slope$p_bh <= significance_threshold) {
    p <- p + geom_abline(intercept = fitted$estimate[fitted$term == "(Intercept)"],
      slope = slope$estimate, linewidth = 0.6)
  }
  figure_panels[[design$panel]] <- p
  save_at_panel_size(p, file.path(figures_dir, paste0("2016_sediment_molecular_", design$panel, ".pdf")))
}
# The richness-difference share of each pair against watercourse distance, with
# the same axes for both data types so their levels can be compared directly.
# Unconnected pairs are drawn first so the fewer flow-connected pairs stay visible.
pairs$connection <- factor(if_else(pairs$flow_connected, connection_levels[1], connection_levels[2]),
  levels = connection_levels)
for (data_type in names(partitions)) {
  panel <- pairs %>% arrange(desc(connection)) %>%
    transmute(watercourse_km, connection, share = .data[[paste0(data_type, "_richness_share")]],
      data_type = data_type_labels[[data_type]])
  p <- ggplot(panel, aes(watercourse_km, share, color = connection)) +
    geom_point(size = pair_point_size, alpha = pair_point_alpha) +
    facet_wrap(~data_type) +
    scale_color_manual(values = connection_colors) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Watercourse distance (km)", y = "Richness-difference share of Jaccard", color = NULL) +
    print_theme + theme(legend.position = "bottom")
  # Every fitted slope is drawn so the four can be compared: solid where the
  # adjusted test passes the threshold, dashed where it does not. Each line spans
  # only the distances of the pairs it was fitted to. Unconnected pairs are drawn
  # first so the flow-connected line stays on top.
  for (pair_set in rev(connection_levels)) {
    tested <- flow_tests %>% filter(.data$data_type == .env$data_type, component == "richness_share",
      separation == "watercourse_km", .data$pair_set == .env$pair_set)
    fitted_range <- range(panel$watercourse_km[panel$connection == pair_set])
    line <- tibble(watercourse_km = fitted_range, share = tested$intercept + tested$slope * fitted_range)
    line_type <- if_else(tested$p_bh <= significance_threshold, "solid", "dashed")
    p <- p + geom_line(data = line, aes(watercourse_km, share), inherit.aes = FALSE,
      color = connection_line_colors[[pair_set]], linetype = line_type, linewidth = 0.6)
  }
  figure_panels[[data_type]] <- p
  save_at_panel_size(p, file.path(figures_dir, paste0("2016_sediment_", data_type, "_richness_share_distance.pdf")))
}
figure <- wrap_plots(figure_panels[c("organic_drainage", "richness_organic", "molecular", "microbial")],
  ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
save_at_panel_size(figure, file.path(manuscript_dir, "2016_sediment_richness_difference.pdf"))

# 12. Save dispersion panels for every view ----------------------------------------
panels <- list()
for (test in seq_len(nrow(test_design))) {
  metric <- test_design$metric[test]
  panel <- dispersion %>% filter(.data$metric == .env$metric)
  sizes <- table(panel$network_group)
  x_labels <- paste0(group_levels, "\n(n = ", as.integer(sizes[group_levels]), ")")
  set.seed(plot_seed)
  p <- ggplot(panel, aes(network_group, distance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.5) +
    geom_point(position = position_jitter(width = 0.1, height = 0, seed = plot_seed),
      size = point_size, alpha = point_alpha) +
    facet_wrap(~metric, labeller = as_labeller(setNames(test_design$label, test_design$metric))) +
    scale_color_manual(values = group_colors) + scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = "Distance to group median") +
    plot_theme + theme(aspect.ratio = panel_aspect)
  panels[[metric]] <- p
  ggsave(file.path(figures_dir, paste0("2016_matched_sediment_dispersion_", metric, ".pdf")),
    p, width = panel_width, height = figure_height, bg = "white")
}
ggsave(file.path(figures_dir, "2016_matched_sediment_dispersion_comparison.pdf"),
  wrap_plots(panels[c("molecular_jaccard", "microbial_hellinger")], ncol = 2) +
    plot_annotation(tag_levels = "A"),
  width = 2 * panel_width, height = figure_height, bg = "white")
ggsave(file.path(figures_dir, "2016_matched_sediment_dispersion_replacement_comparison.pdf"),
  wrap_plots(panels[c("molecular_replacement", "microbial_replacement")], ncol = 2) +
    plot_annotation(tag_levels = "A"),
  width = 2 * panel_width, height = figure_height, bg = "white")

# 13. Align the existing class fractions to the audited 44-site network groups.
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
class_input_paths <- c(
  file.path(class_results_dir, "site_class_fractions.csv"),
  here("results", "fticr_2016", "network_dispersion", "site_network_position_audit.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "tables", "primary_peak_properties.csv"),
  here("figures", "supplement", "molecular", "2016_fticr_site_class_composition.pdf")
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

# 14. Verify that the saved fractions agree with current binary detections.
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

# 15. Average site fractions within groups, giving every site equal weight.
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

# 16. Retain class panels for all sites and for the matched sites.
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
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
  if (!group_name %in% c("All", "Matched")) {
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
stopifnot(identical(class_input_hashes, unname(tools::md5sum(class_input_paths))))

# 17. Confirm that the reused inputs and full-site dispersion results are intact.
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
print(richness_dependence)
print(richness_models, n = Inf)
print(partition_shares)
print(flow_tests, n = Inf)
