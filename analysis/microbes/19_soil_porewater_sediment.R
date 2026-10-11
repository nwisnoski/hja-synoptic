# Analyze the reviewed three-habitat cohort; run 18_prepare_manuscript_cohort.R first.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(permute)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme -----------------------------------------------
data_dir <- here("data", "derived", "soil_porewater_sediment_2016")
results_dir <- here("results", "diversity_2016", "soil_porewater_sediment")
figures_dir <- here("figures")
figure_prefix <- "2016_soil_porewater_sediment_"
target_depth <- 10000L
n_permutations <- 9999L
base_seed <- 20261010L
soil_inventory_draws <- 1000L
soil_inventory_size <- 3L
centerline_stride <- 10L
alpha <- 0.05
habitat_levels <- c("soil", "sediment", "hyporheic")
habitat_labels <- c(soil = "Soil", sediment = "Sediment", hyporheic = "Porewater")
habitat_colors <- c(soil = "#D55E00", sediment = "#009E73", hyporheic = "#0072B2")
group_levels <- c("Headwater", "Intermediate", "Mainstem")
plot_theme <- theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  strip.background = element_blank(), legend.position = "bottom")
habitat_design <- tibble::tribble(
  ~comparison, ~first, ~second, ~seed,
  "All three habitats", NA_character_, NA_character_, base_seed + 1L,
  "Soil-sediment", "soil", "sediment", base_seed + 2L,
  "Soil-porewater", "soil", "hyporheic", base_seed + 3L,
  "Sediment-porewater", "sediment", "hyporheic", base_seed + 4L
)
pair_design <- expand_grid(
  metric = c("Bray-Curtis", "Weighted UniFrac", "Unweighted UniFrac", "Sediment ASV read share"),
  null = c("All paired sites", "Within stream order")
) %>% mutate(seed = base_seed + 100L + row_number())
environment_design <- tibble::tribble(
  ~variable, ~transform,
  "site_drainage_area_ha", "log10",
  "hyporheic_npoc_mg_l", "log10",
  "hyporheic_suva254", "identity",
  "hyporheic_fluorescence_index", "identity",
  "hyporheic_do_mg_l", "identity",
  "hyporheic_water_temperature_c", "identity",
  "hyporheic_tdn_mg_l", "log10",
  "sediment_hydraulic_conductivity_geomean_m_s", "log10"
) %>% mutate(seed = base_seed + 300L + row_number())

# 2. Align the cohort and record checksums for all reused inputs -----------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  file.path(data_dir, "asv_counts_10k.csv"),
  file.path(data_dir, "sample_metadata.csv"),
  file.path(data_dir, "sediment_porewater_pair_audit.csv"),
  here("results", "diversity_2016", "tables", "unifrac_weighted_10k.csv"),
  here("results", "diversity_2016", "tables", "unifrac_unweighted_10k.csv"),
  here("data", "StreamCenterline.csv"),
  here("data", "derived", "soil_drainage_2016", "soil_coordinate_audit.csv"),
  here("results", "fticr_2016", "stable_group_distribution", "core_site_profiles.csv")
)
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
stopifnot(!anyNA(checks$initial_md5))
table <- read_csv(input_paths[1], show_col_types = FALSE)
counts <- t(as.matrix(table[, -1]))
colnames(counts) <- table$asv_id
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()))
stopifnot(identical(metadata$sample_id, rownames(counts)),
  all(metadata$included_manuscript), !any(metadata$excluded_suspect_library),
  nrow(counts) == 76, all(rowSums(counts) == target_depth))
metadata$habitat <- factor(metadata$habitat, levels = habitat_levels)
metadata$network_group <- factor(metadata$network_group, levels = group_levels)
bray <- vegdist(counts, method = "bray")
bray_matrix <- as.matrix(bray)
pair_audit <- read_csv(input_paths[3], col_types = cols(site_code = col_character()))
pairs <- pair_audit %>% filter(paired_manuscript) %>% arrange(site_code)
pairs <- pairs %>% left_join(metadata %>% filter(habitat == "hyporheic") %>%
  select(site_code, site_stream_order, site_segment_number, network_group), by = "site_code")
stopifnot(nrow(pairs) == 18, !anyNA(pairs$site_stream_order))
write_csv(pairs, file.path(results_dir, "paired_site_manifest.csv"))
write_csv(habitat_design, file.path(results_dir, "habitat_test_design.csv"))
write_csv(pair_design, file.path(results_dir, "pairing_test_design.csv"))
write_csv(environment_design, file.path(results_dir, "porewater_environment_design.csv"))

# 3. Refit habitat ordination, diversity, and exploratory distribution tests -----
hellinger <- decostand(counts, method = "hellinger")
pca <- prcomp(hellinger, center = TRUE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
sample_summary <- metadata %>% mutate(
  pc1 = pca$x[, 1], pc2 = pca$x[, 2],
  q0 = specnumber(counts), q1 = exp(diversity(counts, "shannon")),
  q2 = diversity(counts, "invsimpson"))
write_csv(sample_summary, file.path(results_dir, "sample_summary.csv"))
write_csv(tibble(axis = seq_along(variance), percent = variance),
  file.path(results_dir, "hellinger_pca_variance.csv"))
write_csv(sample_summary %>% group_by(habitat) %>% summarize(
  samples = n(), median_q0 = median(q0), median_q1 = median(q1), median_q2 = median(q2),
  .groups = "drop"), file.path(results_dir, "habitat_diversity_summary.csv"))
habitat_rows <- list()
dispersion_rows <- list()
for (i in seq_len(nrow(habitat_design))) {
  current <- metadata
  if (i > 1) current <- metadata %>% filter(habitat %in% c(habitat_design$first[i], habitat_design$second[i]))
  current$habitat <- droplevels(current$habitat)
  d <- as.dist(bray_matrix[current$sample_id, current$sample_id])
  set.seed(habitat_design$seed[i])
  fit <- adonis2(d ~ habitat, data = current, permutations = n_permutations)
  dispersion <- betadisper(d, current$habitat, type = "median", bias.adjust = TRUE)
  set.seed(habitat_design$seed[i] + 1000L)
  dispersion_test <- permutest(dispersion, permutations = n_permutations)
  habitat_rows[[i]] <- tibble(comparison = habitat_design$comparison[i], samples = nrow(current),
    r_squared = fit$R2[1], f = fit$F[1], p = fit$`Pr(>F)`[1],
    dispersion_f = dispersion_test$tab$F[1], dispersion_p = dispersion_test$tab$`Pr(>F)`[1],
    null = "Unrestricted labels; exploratory, spatial and repeated-site dependence not controlled")
  if (i == 1) dispersion_rows[[1]] <- tibble(sample_id = current$sample_id,
    habitat = current$habitat, distance_to_median = dispersion$distances)
}
habitat_tests <- bind_rows(habitat_rows) %>% mutate(
  p_bh_family4 = p.adjust(p, "BH"), dispersion_p_bh_family4 = p.adjust(dispersion_p, "BH"))
write_csv(habitat_tests, file.path(results_dir, "habitat_distribution_tests.csv"))
write_csv(bind_rows(dispersion_rows), file.path(results_dir, "habitat_dispersion.csv"))
# This paired contrast conditions on site and exchanges habitat only within site.
paired_metadata <- metadata %>% filter(sample_id %in% c(pairs$sediment_id, pairs$porewater_id))
paired_metadata$pair_site <- factor(paired_metadata$site_code)
paired_d <- as.dist(bray_matrix[paired_metadata$sample_id, paired_metadata$sample_id])
set.seed(base_seed + 20L)
paired_fit <- adonis2(paired_d ~ pair_site + habitat, data = paired_metadata, by = "margin",
  permutations = how(nperm = n_permutations, blocks = paired_metadata$pair_site))
# Site is a conditioning term and cannot be tested by within-site exchanges.
write_csv(tibble(term = rownames(paired_fit), as.data.frame(paired_fit),
  null = "Habitat labels exchangeable within each of 18 paired sites; site conditioned") %>%
    filter(term == "habitat"),
  file.path(results_dir, "paired_habitat_permanova.csv"))

# 4. Define the soil inventory and quantify coverage sensitivity ----------------
soil_counts <- counts[metadata$habitat == "soil", , drop = FALSE]
receiving <- sample_summary %>% filter(habitat != "soil")
receiving_counts <- counts[receiving$sample_id, , drop = FALSE]
soil_detected <- colSums(soil_counts > 0) > 0
sharing <- receiving %>% mutate(
  soil_read_percent = 100 * rowSums(receiving_counts[, soil_detected, drop = FALSE]) / target_depth,
  soil_asv_percent = 100 * rowSums(receiving_counts[, soil_detected, drop = FALSE] > 0) / q0)
write_csv(tibble(asv_id = colnames(counts), soil_samples_detected = colSums(soil_counts > 0)),
  file.path(results_dir, "soil_inventory.csv"))
# Hold receiving communities fixed and draw three of 15 soils. This isolates
# inventory coverage; it does not reproduce the earlier ten-cell balancing.
draw_rows <- list()
selection_rows <- list()
for (draw_id in seq_len(soil_inventory_draws)) {
  set.seed(base_seed + 2000L + draw_id)
  selected <- sample(seq_len(nrow(soil_counts)), soil_inventory_size)
  detected <- colSums(soil_counts[selected, , drop = FALSE] > 0) > 0
  selection_rows[[draw_id]] <- tibble(draw_id, seed = base_seed + 2000L + draw_id,
    soil_sample_id = rownames(soil_counts)[selected])
  current <- tibble(sample_id = receiving$sample_id, habitat = receiving$habitat,
    read_percent = 100 * rowSums(receiving_counts[, detected, drop = FALSE]) / target_depth,
    asv_percent = 100 * rowSums(receiving_counts[, detected, drop = FALSE] > 0) / receiving$q0)
  draw_rows[[draw_id]] <- current %>% group_by(habitat) %>% summarize(
    median_read_percent = median(read_percent), median_asv_percent = median(asv_percent),
    .groups = "drop") %>% mutate(draw_id, seed = base_seed + 2000L + draw_id,
      soil_samples = soil_inventory_size, status = "ok")
}
draws <- bind_rows(draw_rows)
write_csv(draws, file.path(results_dir, "soil_inventory_draws.csv"))
write_csv(bind_rows(selection_rows), file.path(results_dir, "soil_inventory_selections.csv"))
write_csv(draws %>% group_by(habitat) %>% summarize(
  read_q025 = quantile(median_read_percent, 0.025), read_q975 = quantile(median_read_percent, 0.975),
  asv_q025 = quantile(median_asv_percent, 0.025), asv_q975 = quantile(median_asv_percent, 0.975),
  median_read_percent = median(median_read_percent), median_asv_percent = median(median_asv_percent),
  .groups = "drop"), file.path(results_dir, "soil_inventory_sensitivity_summary.csv"))
write_csv(sharing, file.path(results_dir, "soil_sharing_by_sample.csv"))
write_csv(sharing %>% group_by(habitat) %>% summarize(samples = n(),
  median_read_percent = median(soil_read_percent), median_asv_percent = median(soil_asv_percent),
  .groups = "drop"), file.path(results_dir, "soil_sharing_summary.csv"))
# Average within stream segments before testing drainage slopes. Four slopes
# form one BH family; segment means retain all sites without treating them as replicates.
segment_sharing <- sharing %>% group_by(habitat, site_segment_number) %>% summarize(
  log_drainage = mean(log10(site_drainage_area_ha)),
  soil_read_percent = mean(soil_read_percent), soil_asv_percent = mean(soil_asv_percent),
  sites = n(), .groups = "drop")
slope_rows <- list()
for (habitat_name in c("sediment", "hyporheic")) {
  for (metric in c("soil_read_percent", "soil_asv_percent")) {
    current <- segment_sharing %>% filter(habitat == habitat_name)
    fit <- lm(reformulate("log_drainage", response = metric), data = current)
    slope_rows[[length(slope_rows) + 1L]] <- tibble(habitat = habitat_name, metric,
      segments = nrow(current), sites = sum(current$sites),
      intercept = coef(fit)[1], slope = coef(fit)[2],
      p = summary(fit)$coefficients[2, 4], r_squared = summary(fit)$r.squared,
      status = "Exploratory segment-mean regression; residual spatial dependence possible")
  }
}
slopes <- bind_rows(slope_rows) %>% mutate(p_bh_family4 = p.adjust(p, "BH"))
write_csv(segment_sharing, file.path(results_dir, "soil_sharing_segment_means.csv"))
write_csv(slopes, file.path(results_dir, "soil_sharing_drainage_slopes.csv"))

# 5. Test same-site resemblance against reassigned porewater identities ----------
distance_matrices <- list("Bray-Curtis" = bray_matrix)
for (i in 4:5) {
  saved <- read_csv(input_paths[i], show_col_types = FALSE)
  saved_matrix <- as.matrix(saved[, -1])
  rownames(saved_matrix) <- saved$sample_id
  stopifnot(identical(rownames(saved_matrix), colnames(saved_matrix)),
    all(metadata$sample_id %in% rownames(saved_matrix)),
    isTRUE(all.equal(saved_matrix, t(saved_matrix))), all(diag(saved_matrix) == 0))
  distance_matrices[[c("Weighted UniFrac", "Unweighted UniFrac")[i - 3L]]] <- saved_matrix
}
# Rows are sediment profiles; columns are porewater profiles. A cell records
# the fraction of that porewater library detected in that sediment profile.
read_share <- 100 * (counts[pairs$sediment_id, , drop = FALSE] > 0) %*%
  t(counts[pairs$porewater_id, , drop = FALSE]) / target_depth
distance_matrices[["Sediment ASV read share"]] <- read_share
pair_results <- list()
null_rows <- list()
pair_sample_rows <- list()
for (i in seq_len(nrow(pair_design))) {
  metric <- pair_design$metric[i]
  if (metric == "Sediment ASV read share") {
    cross <- read_share
  } else {
    cross <- distance_matrices[[metric]][pairs$sediment_id, pairs$porewater_id, drop = FALSE]
  }
  observed <- mean(diag(cross))
  set.seed(pair_design$seed[i])
  if (pair_design$null[i] == "Within stream order") {
    permutations <- shuffleSet(nrow(pairs), nset = n_permutations,
      control = how(blocks = factor(pairs$site_stream_order)))
  } else {
    permutations <- shuffleSet(nrow(pairs), nset = n_permutations)
  }
  null_values <- numeric(nrow(permutations))
  for (j in seq_len(nrow(permutations))) {
    null_values[j] <- mean(cross[cbind(seq_len(nrow(pairs)), permutations[j, ])])
  }
  if (metric == "Sediment ASV read share") {
    extreme <- sum(null_values >= observed)
  } else {
    extreme <- sum(null_values <= observed)
  }
  pair_results[[i]] <- tibble(metric, null = pair_design$null[i], pairs = nrow(pairs),
    observed_mean = observed, observed_median = median(diag(cross)),
    offsite_mean = mean(cross[row(cross) != col(cross)]),
    random_assignment_mean = mean(null_values),
    p = (extreme + 1) / (length(null_values) + 1), permutations = length(null_values),
    alternative = ifelse(metric == "Sediment ASV read share", "Greater same-site sharing", "Lower same-site distance"),
    status = "ok")
  null_rows[[i]] <- tibble(metric, null = pair_design$null[i],
    draw_id = seq_along(null_values), seed = pair_design$seed[i], value = null_values)
  if (pair_design$null[i] == "All paired sites") pair_sample_rows[[i]] <- pairs %>% mutate(
    metric, same_site = diag(cross), other_site_mean = (colSums(cross) - diag(cross)) / (nrow(pairs) - 1))
}
pair_tests <- bind_rows(pair_results) %>% group_by(null) %>%
  mutate(p_bh_family4 = p.adjust(p, "BH")) %>% ungroup()
pair_samples <- bind_rows(pair_sample_rows)
write_csv(pair_tests, file.path(results_dir, "same_site_tests.csv"))
write_csv(bind_rows(null_rows), file.path(results_dir, "same_site_null_draws.csv"))
write_csv(pair_samples, file.path(results_dir, "same_site_sample_comparisons.csv"))
# Compare soil sharing within the same 18 sites using sign exchanges, not
# independent sample tests on overlapping or unpaired sites.
soil_difference <- sharing$soil_read_percent[match(pairs$sediment_id, sharing$sample_id)] -
  sharing$soil_read_percent[match(pairs$porewater_id, sharing$sample_id)]
set.seed(base_seed + 200L)
soil_null <- numeric(n_permutations)
for (i in seq_len(n_permutations)) soil_null[i] <- mean(soil_difference * sample(c(-1, 1), nrow(pairs), replace = TRUE))
write_csv(tibble(pairs = nrow(pairs), mean_difference_percentage_points = mean(soil_difference),
  median_difference_percentage_points = median(soil_difference),
  p_two_sided = (1 + sum(abs(soil_null) >= abs(mean(soil_difference)))) / (1 + n_permutations),
  seed = base_seed + 200L, null = "Exchange sediment and porewater sharing within each paired site"),
  file.path(results_dir, "paired_soil_sharing_test.csv"))

# 6. Refit porewater network and chemistry checks after the four exclusions ------
porewater <- metadata %>% filter(habitat == "hyporheic")
porewater_bray <- as.dist(bray_matrix[porewater$sample_id, porewater$sample_id])
set.seed(base_seed + 250L)
network_fit <- adonis2(porewater_bray ~ network_group, data = porewater, permutations = n_permutations)
network_dispersion <- betadisper(porewater_bray, porewater$network_group, type = "median", bias.adjust = TRUE)
set.seed(base_seed + 251L)
network_dispersion_test <- permutest(network_dispersion, permutations = n_permutations)
write_csv(tibble(samples = nrow(porewater), r_squared = network_fit$R2[1],
  p = network_fit$`Pr(>F)`[1], dispersion_p = network_dispersion_test$tab$`Pr(>F)`[1],
  null = "Unrestricted sites; exploratory network association"),
  file.path(results_dir, "porewater_network_tests.csv"))
write_csv(porewater %>% count(network_group, name = "samples"),
  file.path(results_dir, "porewater_network_sample_sizes.csv"))
environment_rows <- list()
environment_audit <- list()
for (i in seq_len(nrow(environment_design))) {
  variable <- environment_design$variable[i]
  x <- porewater[[variable]]
  usable <- is.finite(x)
  if (environment_design$transform[i] == "log10") usable <- usable & x > 0
  environment_audit[[i]] <- tibble(sample_id = porewater$sample_id, variable, included = usable,
    reason = ifelse(usable, "Included", "Missing, nonfinite, or nonpositive for log10"))
  current <- porewater[usable, ]
  value <- x[usable]
  if (environment_design$transform[i] == "log10") value <- log10(value)
  for (restriction in c("All sites", "Within stream order")) {
    if (nrow(current) < 4 || length(unique(value)) < 2) {
      environment_rows[[length(environment_rows) + 1L]] <- tibble(variable, restriction,
        samples = nrow(current), rho = NA_real_, p = NA_real_, status = "Insufficient variable data")
      next
    }
    d <- as.dist(bray_matrix[current$sample_id, current$sample_id])
    set.seed(environment_design$seed[i])
    control <- n_permutations
    if (restriction == "Within stream order") control <- how(nperm = n_permutations,
      blocks = factor(current$site_stream_order))
    fit <- mantel(d, dist(value), method = "spearman", permutations = control)
    environment_rows[[length(environment_rows) + 1L]] <- tibble(variable, restriction,
      samples = nrow(current), rho = unname(fit$statistic), p = fit$signif,
      status = "ok; positive-association Mantel test")
  }
}
environment_tests <- bind_rows(environment_rows) %>% group_by(restriction) %>%
  mutate(p_bh_family8 = p.adjust(p, "BH")) %>% ungroup()
write_csv(environment_tests, file.path(results_dir, "porewater_environment_tests.csv"))
write_csv(bind_rows(environment_audit), file.path(results_dir, "porewater_environment_inclusion.csv"))

# 7. Save individual vector panels and manuscript combinations ------------------
plots <- list()
# Separate habitat panels retain co-located measurements without symbol overlap.
geometry <- read_csv(input_paths[6], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double())) %>% group_by(segment) %>%
  filter(row_number() %% centerline_stride == 1 | row_number() == n()) %>% ungroup()
soil_coordinates <- read_csv(input_paths[7], show_col_types = FALSE)
molecular_sites <- read_csv(input_paths[8], col_types = cols(site_code = col_character())) %>%
  filter(molecular_profile_available)
map_points <- bind_rows(
  metadata %>% filter(habitat != "soil") %>% transmute(sample_id, site_code,
    layer = recode(as.character(habitat), sediment = "Sediment microbes (34)", hyporheic = "Porewater microbes (27)"),
    x = site_utm_x_m, y = site_utm_y_m,
    molecular_pair = habitat == "sediment" & site_code %in% molecular_sites$site_code),
  soil_coordinates %>% filter(sample_id %in% metadata$sample_id, spatial_included) %>%
    transmute(sample_id, site_code, layer = "Soil microbes (13 mapped / 15)",
      x = soil_utm_x_m, y = soil_utm_y_m, molecular_pair = FALSE),
  molecular_sites %>% transmute(sample_id = NA_character_, site_code,
    layer = "Sediment molecular profiles (44)", x = site_utm_x_m, y = site_utm_y_m,
    molecular_pair = site_code %in% metadata$site_code[metadata$habitat == "sediment"])
)
stopifnot(nrow(map_points) == 118, sum(map_points$molecular_pair) == 66,
  all(is.finite(map_points$x)), all(is.finite(map_points$y)))
map_points$layer <- factor(map_points$layer, levels = c("Soil microbes (13 mapped / 15)",
  "Porewater microbes (27)", "Sediment microbes (34)", "Sediment molecular profiles (44)"))
write_csv(map_points, file.path(results_dir, "sampling_map_points.csv"))
plots$sampling_map <- ggplot(geometry, aes(x / 1000, y / 1000, group = segment)) +
  geom_path(color = "gray75", linewidth = 0.3) +
  geom_point(data = map_points, aes(x / 1000, y / 1000, color = layer),
    inherit.aes = FALSE, size = 2.2) +
  geom_point(data = map_points %>% filter(molecular_pair), aes(x / 1000, y / 1000),
    inherit.aes = FALSE, shape = 1, size = 3.3, color = "black", stroke = 0.4) +
  facet_wrap(~layer, ncol = 2, drop = FALSE) + coord_equal() +
  scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")) +
  labs(x = "Easting (km, UTM zone 10N)", y = "Northing (km)") +
  plot_theme + theme(legend.position = "none")
plots$hellinger_pca <- ggplot(sample_summary, aes(pc1, pc2, color = habitat, shape = habitat)) +
  geom_point(size = 2.3, alpha = 0.85) + coord_equal() +
  scale_color_manual(values = habitat_colors, labels = habitat_labels) +
  scale_shape_manual(values = c(soil = 17, sediment = 16, hyporheic = 15), labels = habitat_labels) +
  labs(x = sprintf("Hellinger PC1 (%.1f%%)", variance[1]),
    y = sprintf("Hellinger PC2 (%.1f%%)", variance[2]), color = "Habitat", shape = "Habitat") + plot_theme
for (order in c("q0", "q1", "q2")) {
  plots[[paste0("diversity_", order)]] <- ggplot(sample_summary, aes(habitat, .data[[order]], color = habitat)) +
    geom_boxplot(outlier.shape = NA, color = "gray40") +
    geom_point(position = position_jitter(width = 0.12, seed = base_seed), size = 1.5) +
    scale_color_manual(values = habitat_colors) + scale_x_discrete(labels = habitat_labels) +
    labs(x = NULL, y = paste0("Hill diversity (", order, "), 10K reads")) +
    plot_theme + theme(legend.position = "none")
}
for (metric in c("soil_read_percent", "soil_asv_percent")) {
  plot <- ggplot(sharing, aes(site_drainage_area_ha, .data[[metric]], color = habitat)) +
    geom_point(size = 2.2, alpha = 0.85) + scale_x_log10() +
    scale_color_manual(values = habitat_colors, labels = habitat_labels) +
    facet_wrap(~habitat, labeller = as_labeller(habitat_labels)) +
    labs(x = "Drainage area (ha)", y = ifelse(metric == "soil_read_percent",
      "Reads shared with soil inventory (%)", "ASVs shared with soil inventory (%)"), color = "Habitat") +
    plot_theme + theme(legend.position = "none")
  significant <- slopes %>% filter(.data$metric == .env$metric, p_bh_family4 <= alpha)
  if (nrow(significant) > 0) plot <- plot + geom_abline(data = significant,
    aes(intercept = intercept, slope = slope, color = habitat), inherit.aes = FALSE)
  plots[[metric]] <- plot
}
distance_data <- pair_samples %>% filter(metric != "Sediment ASV read share") %>%
  pivot_longer(c(same_site, other_site_mean), names_to = "comparison", values_to = "distance") %>%
  mutate(comparison = factor(comparison, levels = c("same_site", "other_site_mean")))
plots$same_site_distances <- ggplot(distance_data, aes(comparison, distance, group = site_code)) +
  geom_line(color = "gray70") + geom_point(size = 1.8, color = "#0072B2") +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_discrete(labels = c(same_site = "Same site", other_site_mean = "Other sites\n(mean)")) +
  labs(x = NULL, y = "Sediment-porewater dissimilarity") + plot_theme
read_data <- pair_samples %>% filter(metric == "Sediment ASV read share") %>%
  pivot_longer(c(same_site, other_site_mean), names_to = "comparison", values_to = "read_percent") %>%
  mutate(comparison = factor(comparison, levels = c("same_site", "other_site_mean")))
plots$same_site_read_sharing <- ggplot(read_data, aes(comparison, read_percent, group = site_code)) +
  geom_line(color = "gray70") + geom_point(size = 2, color = "#0072B2") +
  scale_x_discrete(labels = c(same_site = "Same site", other_site_mean = "Other sites\n(mean)")) +
  labs(x = NULL, y = "Porewater reads detected in sediment (%)") + plot_theme
plots$habitat_context <- (plots$hellinger_pca | plots$diversity_q1) + plot_layout(widths = c(1.4, 1)) +
  plot_annotation(tag_levels = "A")
plots$soil_sharing <- (plots$soil_read_percent / plots$soil_asv_percent) + plot_annotation(tag_levels = "A")
plots$same_site_comparison <- (plots$same_site_distances / plots$same_site_read_sharing) +
  plot_annotation(tag_levels = "A")
figure_design <- tibble::tribble(
  ~panel, ~width, ~height,
  "sampling_map", 9, 9,
  "hellinger_pca", 7, 5.5,
  "diversity_q0", 5.5, 4.6,
  "diversity_q1", 5.5, 4.6,
  "diversity_q2", 5.5, 4.6,
  "soil_read_percent", 9, 4.5,
  "soil_asv_percent", 9, 4.5,
  "same_site_distances", 10, 4.5,
  "same_site_read_sharing", 5.5, 4.5,
  "habitat_context", 12, 5.5,
  "soil_sharing", 9, 8,
  "same_site_comparison", 10, 8
)
for (i in seq_len(nrow(figure_design))) {
  ggsave(file.path(figures_dir, paste0(figure_prefix, figure_design$panel[i], ".pdf")),
    plots[[figure_design$panel[i]]], width = figure_design$width[i], height = figure_design$height[i], bg = "white")
}
write_csv(figure_design, file.path(results_dir, "figure_manifest.csv"))

# 8. Check provenance and save settings/session information ---------------------
checks$final_md5 <- unname(tools::md5sum(checks$path))
stopifnot(identical(checks$initial_md5, checks$final_md5))
write_csv(checks, file.path(results_dir, "input_checksums.csv"))
write_csv(tibble(target_depth, n_permutations, base_seed, soil_inventory_draws, soil_inventory_size, centerline_stride, alpha),
  file.path(results_dir, "analysis_settings.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
print(habitat_tests)
print(pair_tests)
print(environment_tests)
message("Three-habitat analysis complete: 76 samples, 18 pairs, 12 vector figures.")
