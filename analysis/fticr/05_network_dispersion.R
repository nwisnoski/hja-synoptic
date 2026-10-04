# Molecular heterogeneity across the 2016 stream network; no intensity weighting.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "network_dispersion")
figures_dir <- here("figures")
null_seed <- 20261007
permutation_seed <- 20261008
resampling_seed <- 20261009
plot_seed <- 20261010
reuse_matching_null <- TRUE
n_null <- 1999
null_batch_size <- 20
n_permutations <- 9999
n_balanced_draws <- 1000
headwater_orders <- c(1, 2)
mainstem_order <- 5
chemistry_variables <- c("mean_detected_h_c", "mean_detected_o_c",
  "nitrogen_fraction", "protein_like_fraction", "tannin_like_fraction")
group_levels <- c("Headwater", "Intermediate", "Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
metric_labels <- c(jaccard = "Binary Jaccard", replacement = "Jaccard replacement",
  equal_features = "Equal-feature Jaccard", chemical_properties = "Formula properties")
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "bottom")

# 2. Match site identities and derive audited network positions ------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
pa_table <- read_csv(here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  col_types = cols(site_code = col_character()))
pa <- as.matrix(pa_table[, -1])
rownames(pa) <- pa_table$site_code
pa <- pa[, colSums(pa) > 0, drop = FALSE]
site <- read_csv(here("results", "fticr_2016", "river_composition", "site_formula_composition.csv"),
  col_types = cols(site_code = col_character()))
site <- site[match(rownames(pa), site$site_code), ]
stopifnot(nrow(site) == 44, ncol(pa) == 4718, identical(site$site_code, rownames(pa)),
  !anyDuplicated(site$site_code), all(pa %in% c(0, 1)),
  identical(as.numeric(rowSums(pa)), site$detected_primary_features))
site$network_group <- factor(ifelse(site$site_stream_order %in% headwater_orders,
  "Headwater", ifelse(site$site_stream_order == mainstem_order, "Mainstem", "Intermediate")),
  levels = group_levels)

# The raw header has a trailing comma but data rows have seven fields.
# D is basin-outlet stream distance; the site's original distance is valley-local.
geometry <- read_csv(here("data", "StreamCenterline.csv"), skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
# The outlet's elevation is NaN in the source; elevation is not used here.
stopifnot(nrow(problems(geometry)) == 0,
  all(is.finite(as.matrix(geometry[, c("x", "y", "outlet_distance", "upstream_area", "order", "segment")]))))
write_csv(tibble(variable = names(geometry),
  nonfinite_values = colSums(!is.finite(as.matrix(geometry)))),
  file.path(results_dir, "geometry_missingness_audit.csv"))
network_audit <- site %>% select(site_code, network_group, site_stream_order,
  site_segment_number, site_utm_x_m, site_utm_y_m, site_drainage_area_ha,
  detected_primary_features, site_distance_to_outlet_m)
network_audit$centerline_match_error_m <- NA_real_
network_audit$basin_outlet_stream_distance_m <- NA_real_
network_audit$centerline_stream_order <- NA_real_
network_audit$centerline_upstream_area_ha <- NA_real_
for (i in seq_len(nrow(site))) {
  candidates <- which(geometry$segment == site$site_segment_number[i])
  stopifnot(length(candidates) > 0)
  distance_squared <- (geometry$x[candidates] - site$site_utm_x_m[i])^2 +
    (geometry$y[candidates] - site$site_utm_y_m[i])^2
  nearest <- candidates[which.min(distance_squared)]
  network_audit$centerline_match_error_m[i] <- sqrt(min(distance_squared))
  network_audit$basin_outlet_stream_distance_m[i] <- geometry$outlet_distance[nearest]
  network_audit$centerline_stream_order[i] <- geometry$order[nearest]
  network_audit$centerline_upstream_area_ha[i] <- geometry$upstream_area[nearest] / 10000
}
network_audit$order_matches <- network_audit$site_stream_order == network_audit$centerline_stream_order
write_csv(network_audit, file.path(results_dir, "site_network_position_audit.csv"))
write_csv(site %>% group_by(network_group) %>%
  summarize(sites = n(), segments = n_distinct(site_segment_number),
    minimum_features = min(detected_primary_features), maximum_features = max(detected_primary_features),
    .groups = "drop"), file.path(results_dir, "group_sampling_summary.csv"))

# 3. Molecular incidence and standardized chemical-property distances -----------
jaccard <- as.matrix(vegdist(pa, method = "jaccard", binary = TRUE))
components <- read_csv(here("results", "fticr_2016", "river_composition",
  "molecular_pairwise_components.csv"), col_types = cols(site_1 = col_character(), site_2 = col_character()))
replacement <- matrix(0, nrow(site), nrow(site), dimnames = dimnames(jaccard))
for (i in seq_len(nrow(components))) {
  a <- match(components$site_1[i], site$site_code)
  b <- match(components$site_2[i], site$site_code)
  replacement[a, b] <- replacement[b, a] <- components$replacement[i]
  stopifnot(abs(jaccard[a, b] - components$jaccard[i]) < 1e-12)
}
thinned <- read_csv(here("results", "fticr_2016", "river_composition",
  "mean_subsampled_jaccard.csv"), col_types = cols(site_code = col_character()))
thinned <- thinned[match(site$site_code, thinned$site_code), ]
equal_features <- as.matrix(thinned[, site$site_code])
dimnames(equal_features) <- dimnames(jaccard)
chemical_matrix <- scale(as.matrix(site[, chemistry_variables]))
rownames(chemical_matrix) <- site$site_code
write_csv(tibble(variable = chemistry_variables, center = attr(chemical_matrix, "scaled:center"),
  scale = attr(chemical_matrix, "scaled:scale")), file.path(results_dir, "chemical_property_scaling.csv"))
distances <- list(jaccard = jaccard, replacement = replacement,
  equal_features = equal_features, chemical_properties = as.matrix(dist(chemical_matrix)))

# 4. Richness-constrained Raup-Crick using one regional 44-site pool --------------
# r1 keeps each site's observed feature count and draws without replacement,
# weighted by feature frequency across ALL 44 sites. Groups use the same pool.
# Half of ties is explicit, matching the PNNL presence/absence RC workflow.
# Shared counts decrease as dissimilarity increases: null shared > observed
# contributes to a positive RC. Values run from -1 (similar) to +1 (different).
triangle <- lower.tri(jaccard)
key_path <- file.path(results_dir, "null_input_checksum.csv")
current_null_key <- tibble(input_md5 = unname(tools::md5sum(here("data", "derived",
  "fticr_2016", "sediment_44_presence_absence.csv"))), null_model = "r1",
  tie_rule = "midpoint", simulations = n_null, seed = null_seed)
null_files <- file.path(results_dir, c("null_input_checksum.csv", "null_model_audit.csv",
  "raup_crick_signed_matrix.csv", "raup_crick_pairwise_diagnostics.csv"))
reuse_null <- FALSE
if (reuse_matching_null && all(file.exists(null_files))) {
  saved_key <- read_csv(key_path, show_col_types = FALSE)
  reuse_null <- isTRUE(all.equal(as.data.frame(saved_key), as.data.frame(current_null_key)))
}
if (reuse_null) {
  rc_table <- read_csv(file.path(results_dir, "raup_crick_signed_matrix.csv"),
    col_types = cols(site_code = col_character()))
  stopifnot(identical(rc_table$site_code, site$site_code))
  rc <- as.matrix(rc_table[, site$site_code])
  dimnames(rc) <- dimnames(jaccard)
  null_diagnostics <- read_csv(file.path(results_dir, "raup_crick_pairwise_diagnostics.csv"),
    col_types = cols(site_1 = col_character(), site_2 = col_character()))
  stopifnot(nrow(null_diagnostics) == sum(triangle), all(is.finite(rc)),
    max(abs(rc - t(rc))) < 1e-12, all(rc >= -1 & rc <= 1),
    max(abs(as.vector(rc[triangle]) - null_diagnostics$rc_signed)) < 1e-12)
  message("Reusing checksum-, seed-, and simulation-matched null outputs.")
} else {
  set.seed(null_seed)
  null_check <- simulate(nullmodel(pa, "r1"), nsim = 3)
  stopifnot(all(apply(null_check, c(1, 3), sum) == rowSums(pa)), all(null_check %in% c(0, 1)))
  set.seed(null_seed)
  null_timing <- system.time(null_fit <- oecosimu(pa,
    function(x) tcrossprod(x)[triangle], method = "r1", nsimul = n_null,
    batchsize = null_batch_size, parallel = 1, alternative = "greater"))
  observed_shared <- as.vector(tcrossprod(pa)[triangle])
  simulated_shared <- null_fit$oecosimu$simulated
  stopifnot(identical(as.numeric(null_fit$statistic), observed_shared),
    all(dim(simulated_shared) == c(length(observed_shared), n_null)))
  greater_fraction <- rowMeans(simulated_shared > observed_shared)
  tie_fraction <- rowMeans(simulated_shared == observed_shared)
  rc_signed <- 2 * (greater_fraction + 0.5 * tie_fraction) - 1
  rc <- matrix(0, nrow(site), nrow(site), dimnames = dimnames(jaccard))
  rc[triangle] <- rc_signed
  rc <- rc + t(rc)
  write_csv(tibble(site_code = rownames(rc), as_tibble(rc, .name_repair = "minimal")),
    file.path(results_dir, "raup_crick_signed_matrix.csv"))
  null_diagnostics <- tibble(site_1 = rownames(pa)[row(jaccard)[triangle]],
    site_2 = rownames(pa)[col(jaccard)[triangle]], observed_shared = observed_shared,
    null_mean_shared = rowMeans(simulated_shared), rc_signed = rc_signed,
    tie_fraction = tie_fraction,
    monte_carlo_se = 2 * apply((simulated_shared > observed_shared) +
      0.5 * (simulated_shared == observed_shared), 1, sd) / sqrt(n_null))
  write_csv(null_diagnostics, file.path(results_dir, "raup_crick_pairwise_diagnostics.csv"))
  write_csv(tibble(null_model = "r1", sites = nrow(pa), regional_features = ncol(pa),
    simulations = n_null, seed = null_seed, elapsed_seconds = unname(null_timing["elapsed"]),
    saturated_minus_one = sum(rc_signed == -1), saturated_plus_one = sum(rc_signed == 1)),
    file.path(results_dir, "null_model_audit.csv"))
  write_csv(current_null_key, key_path)
  rm(null_fit, simulated_shared, null_check)
}

# 5. Three-group dispersion and site-label RC tests ------------------------------
# All 44 sites enter primary three-group tests. Pairwise contrasts are secondary.
# Spatial median, bias adjustment, and all dimensions are retained throughout.
set.seed(permutation_seed)
permutations <- shuffleSet(nrow(site), nset = n_permutations)
contrast_names <- c("Headwater-Intermediate", "Headwater-Mainstem", "Intermediate-Mainstem")
dispersion_rows <- list()
omnibus_rows <- list()
contrast_rows <- list()
correction_rows <- list()
corrected_distances <- list()
for (metric in names(distances)) {
  raw <- distances[[metric]]
  ordination <- wcmdscale(as.dist(raw), eig = TRUE, add = "lingoes")
  correction <- ordination$ac
  corrected <- sqrt(raw^2 + 2 * correction)
  diag(corrected) <- 0
  corrected_distances[[metric]] <- corrected
  correction_rows[[metric]] <- tibble(metric = metric, lingoes_constant = correction)
  fit <- betadisper(as.dist(corrected), site$network_group,
    type = "median", bias.adjust = TRUE)
  dispersion_rows[[metric]] <- tibble(site_code = site$site_code,
    network_group = site$network_group, metric = metric, distance = fit$distances)
  permutation_result <- permutest(fit, permutations = permutations, pairwise = TRUE, parallel = 1)
  means <- tapply(fit$distances, site$network_group, mean)
  omnibus_rows[[metric]] <- tibble(metric = metric, test = "Three-group PERMDISP",
    sites = nrow(site), headwater_mean = means["Headwater"],
    intermediate_mean = means["Intermediate"], mainstem_mean = means["Mainstem"],
    statistic = permutation_result$tab[1, "F"], p_value = permutation_result$tab[1, "Pr(>F)"])
  for (contrast in contrast_names) {
    groups <- strsplit(contrast, "-", fixed = TRUE)[[1]]
    contrast_rows[[length(contrast_rows) + 1]] <- tibble(metric = metric, contrast = contrast,
      first_mean = means[groups[1]], second_mean = means[groups[2]],
      difference = means[groups[1]] - means[groups[2]], ratio = means[groups[1]] / means[groups[2]],
      p_two_sided = permutation_result$pairwise$permuted[contrast], p_first_greater = NA_real_)
  }
}
dispersion <- bind_rows(dispersion_rows)
write_csv(dispersion, file.path(results_dir, "site_dispersion.csv"))
write_csv(bind_rows(correction_rows), file.path(results_dir, "distance_correction_audit.csv"))

# RC group means compare WITHIN-group pairs. Permute sites across all three
# groups, preserving 24/14/6 labels. No pair is an independent replicate.
observed_rc_means <- numeric(length(group_levels))
names(observed_rc_means) <- group_levels
for (group in group_levels) {
  indices <- which(site$network_group == group)
  observed_rc_means[group] <- mean(as.dist(rc[indices, indices]))
}
observed_omnibus <- sum((observed_rc_means - mean(observed_rc_means))^2)
rc_permutation_means <- matrix(NA_real_, n_permutations, length(group_levels),
  dimnames = list(NULL, group_levels))
for (i in seq_len(n_permutations)) {
  labels <- site$network_group[permutations[i, ]]
  for (group in group_levels) {
    indices <- which(labels == group)
    rc_permutation_means[i, group] <- mean(as.dist(rc[indices, indices]))
  }
}
rc_permutation_omnibus <- rowSums((rc_permutation_means - rowMeans(rc_permutation_means))^2)
omnibus_rows$raup_crick <- tibble(metric = "raup_crick", test = "Three-group site-label RC spread",
  sites = nrow(site), headwater_mean = observed_rc_means["Headwater"],
  intermediate_mean = observed_rc_means["Intermediate"], mainstem_mean = observed_rc_means["Mainstem"],
  statistic = observed_omnibus,
  p_value = (1 + sum(rc_permutation_omnibus >= observed_omnibus)) / (1 + n_permutations))
for (contrast in contrast_names) {
  groups <- strsplit(contrast, "-", fixed = TRUE)[[1]]
  observed <- observed_rc_means[groups[1]] - observed_rc_means[groups[2]]
  simulated <- rc_permutation_means[, groups[1]] - rc_permutation_means[, groups[2]]
  contrast_rows[[length(contrast_rows) + 1]] <- tibble(metric = "raup_crick", contrast = contrast,
    first_mean = observed_rc_means[groups[1]], second_mean = observed_rc_means[groups[2]],
    difference = observed, ratio = NA_real_,
    p_two_sided = (1 + sum(abs(simulated) >= abs(observed))) / (1 + n_permutations),
    p_first_greater = (1 + sum(simulated >= observed)) / (1 + n_permutations))
}
tests <- bind_rows(omnibus_rows) %>% mutate(p_bh_five_tests = p.adjust(p_value, method = "BH"))
contrasts <- bind_rows(contrast_rows) %>% mutate(p_bh_fifteen_contrasts = p.adjust(p_two_sided, method = "BH"))
write_csv(tests, file.path(results_dir, "three_group_tests.csv"))
write_csv(contrasts, file.path(results_dir, "pairwise_group_contrasts.csv"))
write_csv(filter(contrasts, contrast == "Headwater-Mainstem"),
  file.path(results_dir, "headwater_mainstem_tests.csv"))
write_csv(tibble(permutation = seq_len(n_permutations),
  as_tibble(rc_permutation_means), omnibus_spread = rc_permutation_omnibus),
  file.path(results_dir, "raup_crick_site_label_permutations.csv"))

# Pair tables include geography to make differences in spatial extent inspectable.
pairwise <- null_diagnostics %>%
  mutate(group_1 = site$network_group[match(site_1, site$site_code)],
    group_2 = site$network_group[match(site_2, site$site_code)],
    within_group = ifelse(group_1 == group_2, as.character(group_1), NA_character_),
    geographic_distance_m = as.vector(dist(as.matrix(site[, c("site_utm_x_m", "site_utm_y_m")]))),
    jaccard = as.vector(as.dist(jaccard)), replacement = as.vector(as.dist(replacement)))
write_csv(pairwise, file.path(results_dir, "molecular_network_pairs.csv"))
write_csv(pairwise %>% filter(!is.na(within_group)) %>% group_by(within_group) %>%
  summarize(pairs = n(), mean_rc = mean(rc_signed), median_rc = median(rc_signed),
    mean_geographic_distance_m = mean(geographic_distance_m), .groups = "drop"),
  file.path(results_dir, "within_group_pair_summary.csv"))

# 6. Balance site numbers and avoid sampling the same headwater segment twice ----
# Each draw selects six distinct headwater AND intermediate segments, then
# one site per segment. All three groups contribute six sites.
# All six mainstem sites (each on a distinct segment) are retained. These are
# sensitivity effect distributions, not confidence intervals or independent tests.
headwater_indices <- which(site$network_group == "Headwater")
headwater_segments <- unique(site$site_segment_number[headwater_indices])
intermediate_indices <- which(site$network_group == "Intermediate")
intermediate_segments <- unique(site$site_segment_number[intermediate_indices])
mainstem_indices <- which(site$network_group == "Mainstem")
stopifnot(length(mainstem_indices) == 6,
  length(unique(site$site_segment_number[mainstem_indices])) == 6,
  length(headwater_segments) == 16, length(intermediate_segments) == 13)
balanced_rows <- list()
selection_rows <- list()
for (run in seq_len(n_balanced_draws)) {
  set.seed(resampling_seed + run)
  segments <- sample(headwater_segments, length(mainstem_indices))
  selected_head <- integer(length(segments))
  for (j in seq_along(segments)) {
    candidates <- headwater_indices[site$site_segment_number[headwater_indices] == segments[j]]
    selected_head[j] <- candidates[sample.int(length(candidates), 1)]
  }
  middle_segments <- sample(intermediate_segments, length(mainstem_indices))
  selected_middle <- integer(length(middle_segments))
  for (j in seq_along(middle_segments)) {
    candidates <- intermediate_indices[site$site_segment_number[intermediate_indices] == middle_segments[j]]
    selected_middle[j] <- candidates[sample.int(length(candidates), 1)]
  }
  selected <- c(selected_head, selected_middle, mainstem_indices)
  labels <- droplevels(site$network_group[selected])
  selection_rows[[run]] <- tibble(run = run, seed = resampling_seed + run,
    site_code = site$site_code[selected], network_group = labels,
    segment = site$site_segment_number[selected])
  for (metric in names(corrected_distances)) {
    fit <- betadisper(as.dist(corrected_distances[[metric]][selected, selected]),
      labels, type = "median", bias.adjust = TRUE)
    means <- tapply(fit$distances, labels, mean)
    balanced_rows[[length(balanced_rows) + 1]] <- tibble(run = run, metric = metric,
      headwater_mean = means["Headwater"], intermediate_mean = means["Intermediate"],
      mainstem_mean = means["Mainstem"], intermediate_mainstem_ratio = means["Intermediate"] / means["Mainstem"],
      headwater_intermediate_ratio = means["Headwater"] / means["Intermediate"],
      difference = means["Headwater"] - means["Mainstem"],
      ratio = means["Headwater"] / means["Mainstem"])
  }
  h <- mean(as.dist(rc[selected_head, selected_head]))
  middle <- mean(as.dist(rc[selected_middle, selected_middle]))
  m <- mean(as.dist(rc[mainstem_indices, mainstem_indices]))
  balanced_rows[[length(balanced_rows) + 1]] <- tibble(run = run, metric = "raup_crick",
    headwater_mean = h, intermediate_mean = middle, mainstem_mean = m,
    intermediate_mainstem_ratio = NA_real_, headwater_intermediate_ratio = NA_real_,
    difference = h - m, ratio = NA_real_)
}
balanced <- bind_rows(balanced_rows)
stopifnot(nrow(balanced) == n_balanced_draws * 5, all(is.finite(balanced$difference)))
write_csv(balanced, file.path(results_dir, "balanced_segment_draws.csv"))
write_csv(bind_rows(selection_rows), file.path(results_dir, "balanced_site_selection.csv"))
write_csv(balanced %>% group_by(metric) %>% summarize(draws = n(),
  median_difference = median(difference), min_difference = min(difference),
  max_difference = max(difference), fraction_headwater_greater = mean(difference > 0),
  median_ratio = median(ratio),
  median_intermediate_mainstem_ratio = median(intermediate_mainstem_ratio),
  median_headwater_intermediate_ratio = median(headwater_intermediate_ratio),
  fraction_intermediate_greater_mainstem = mean(intermediate_mean > mainstem_mean),
  fraction_headwater_greater_intermediate = mean(headwater_mean > intermediate_mean), .groups = "drop"),
  file.path(results_dir, "balanced_segment_summary.csv"))

# 7. Individual vector figures and reproducibility record -----------------------
# Only the plotted centerline is reduced to approximately 10-m vertex spacing.
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% 10 == 1 | row_number() == n()) %>% ungroup()
p_map <- ggplot(map_lines, aes(x, y, group = segment)) +
  geom_path(color = "#B0B0B0", linewidth = 0.35) +
  geom_point(data = site, aes(site_utm_x_m, site_utm_y_m, color = network_group),
    inherit.aes = FALSE, size = 2.2) +
  scale_color_manual(values = group_colors, name = NULL) + coord_equal() +
  labs(x = "UTM easting (m)", y = "UTM northing (m)") + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_network_groups.pdf"), p_map,
  width = 8, height = 6, bg = "white")
set.seed(plot_seed)
for (metric in c("jaccard", "equal_features", "chemical_properties")) {
  panel <- filter(dispersion, .data$metric == .env$metric)
  p <- ggplot(panel, aes(network_group, distance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.5) +
    geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.85) +
    scale_color_manual(values = group_colors, guide = "none") +
    labs(x = NULL, y = "Distance to group spatial median") + plot_theme
  ggsave(file.path(figures_dir, paste0("2016_fticr_dispersion_", metric, ".pdf")),
    p, width = 6, height = 4.5, bg = "white")
}
within_pairs <- pairwise %>% filter(!is.na(within_group)) %>%
  mutate(within_group = factor(within_group, levels = group_levels))
p_rc <- ggplot(within_pairs, aes(within_group, rc_signed, color = within_group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#777777") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.4) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = NULL, y = "Within-group Raup-Crick (signed)") + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_network_raup_crick.pdf"), p_rc,
  width = 6, height = 4.5, bg = "white")
balanced_plot <- balanced %>% filter(metric != "raup_crick") %>%
  select(run, metric, ratio, intermediate_mainstem_ratio, headwater_intermediate_ratio) %>%
  pivot_longer(c(ratio, intermediate_mainstem_ratio, headwater_intermediate_ratio),
    names_to = "comparison", values_to = "ratio_value") %>%
  mutate(comparison = factor(comparison,
    levels = c("ratio", "intermediate_mainstem_ratio", "headwater_intermediate_ratio"),
    labels = c("Headwater / mainstem", "Intermediate / mainstem", "Headwater / intermediate")))
p_balanced <- ggplot(balanced_plot, aes(factor(metric,
  levels = names(metric_labels)), ratio_value, fill = comparison)) +
  geom_hline(yintercept = 1, color = "#777777", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.3, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73"), name = NULL) +
  scale_x_discrete(labels = metric_labels) + coord_flip() +
  labs(x = NULL, y = "Dispersion ratio (6 sites per group)") + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_balanced_dispersion.pdf"), p_balanced,
  width = 8, height = 5.5, bg = "white")
capture.output(sessionInfo(), file = file.path(results_dir, "session_info.txt"))
print(tests)
print(read_csv(file.path(results_dir, "balanced_segment_summary.csv"), show_col_types = FALSE))
