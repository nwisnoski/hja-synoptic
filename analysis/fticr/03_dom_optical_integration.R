# Compare water-column DOM optical signatures with sediment molecular incidence.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "dom_optical_integration")
derived_dir <- here("data", "derived", "fticr_2016")
figures_dir <- here("figures")
n_permutations <- 9999
permutation_seed <- 20261005
optical_variables <- c("suva254", "fluorescence_index", "log_t_c")
formula_traits <- c("mean_detected_h_c", "mean_detected_o_c", "protein_like_fraction",
  "tannin_like_fraction", "lignin_like_fraction")
water_compartments <- c("surface", "hyporheic")
compartment_colors <- c(surface = "#0072B2", hyporheic = "#D55E00")
compartment_labels <- c(surface = "Surface water", hyporheic = "Hyporheic water")
variable_labels <- c(suva254 = "SUVA254", fluorescence_index = "Fluorescence index",
  log_t_c = "log10 peak T:C")
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank())

# 2. Audit source measurements and construct optical summaries ------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
environment <- read_csv(here("data", "derived", "analysis_inputs", "environment",
  "aquatic_59_site_environment.csv"), col_types = cols(site_code = col_character()))
master <- read_csv(here("data", "hja-env_data_clean.csv"),
  col_types = cols(`Site Code` = col_character()), show_col_types = FALSE)
dictionary <- read_csv(here("data", "derived", "analysis_inputs", "audit",
  "environment_variable_dictionary.csv"), show_col_types = FALSE)
stopifnot(nrow(environment) == 59, !anyDuplicated(environment$site_code),
  all(is.finite(environment$site_drainage_area_ha)), all(environment$site_drainage_area_ha > 0))
source_rules <- dictionary %>%
  filter(analysis_column %in% names(environment),
    grepl("^(surface|hyporheic)_(npoc|suva|fluorescence|eem|spectral)", analysis_column))
source_audit <- list()
raw_key <- paste(master$`Site Code`, master$`Sample Type`, sep = "::")
for (j in seq_len(nrow(source_rules))) {
  rule <- source_rules[j, ]
  matched <- match(paste(environment$site_code, rule$source_sample_type, sep = "::"), raw_key)
  expected <- as.numeric(master[[rule$source_column]][matched])
  observed <- environment[[rule$analysis_column]]
  stopifnot(identical(is.na(expected), is.na(observed)),
    all(abs(expected - observed) < 1e-12, na.rm = TRUE))
  source_audit[[j]] <- tibble(analysis_column = rule$analysis_column,
    source_column = rule$source_column, source_sample_type = rule$source_sample_type,
    available_sites = sum(is.finite(observed)), source_values_match = TRUE)
}
write_csv(bind_rows(source_audit), file.path(results_dir, "optical_source_checks.csv"))

optical_rows <- list()
for (compartment in water_compartments) {
  optical_rows[[compartment]] <- environment %>%
    transmute(site_code, compartment = compartment, site_drainage_area_ha,
      site_stream_order, site_distance_to_outlet_m,
      npoc_mg_l = .data[[paste0(compartment, "_npoc_mg_l")]],
      suva254 = .data[[paste0(compartment, "_suva254")]],
      fluorescence_index = .data[[paste0(compartment, "_fluorescence_index")]],
      peak_a_ru = .data[[paste0(compartment, "_eem_peak_a_ru")]],
      peak_c_ru = .data[[paste0(compartment, "_eem_peak_c_ru")]],
      peak_t_ru = .data[[paste0(compartment, "_eem_peak_t_ru")]],
      source_spectral_slope_ratio = .data[[paste0(compartment, "_spectral_slope_ratio")]])
}
optical <- bind_rows(optical_rows)
optical$log_area <- log10(optical$site_drainage_area_ha)
optical$source_slope_ratio_275_295_to_300_350 <- NA_real_
optical$source_slope_ratio_275_295_to_300_350[optical$compartment == "surface"] <-
  environment$surface_spectral_slope_ratio_275_295_to_300_350
# Ratios are optical contrasts, not fractions of carbon or identified compounds.
optical$t_c_ratio <- optical$peak_t_ru / optical$peak_c_ru
optical$t_a_ratio <- optical$peak_t_ru / optical$peak_a_ru
optical$log_t_c <- log10(optical$t_c_ratio)
optical$log_t_a <- log10(optical$t_a_ratio)
optical$peak_a_per_npoc <- optical$peak_a_ru / optical$npoc_mg_l
optical$peak_c_per_npoc <- optical$peak_c_ru / optical$npoc_mg_l
optical$peak_t_per_npoc <- optical$peak_t_ru / optical$npoc_mg_l
optical$complete_primary_optics <- complete.cases(optical[optical_variables])
audit <- optical %>% select(site_code, compartment, complete_primary_optics)
audit$missing_primary_variables <- ""
for (i in seq_len(nrow(optical))) {
  missing <- optical_variables[!is.finite(as.numeric(optical[i, optical_variables]))]
  audit$missing_primary_variables[i] <- paste(missing, collapse = ";")
}
write_csv(audit, file.path(results_dir, "optical_site_inclusion_audit.csv"))
stopifnot(sum(optical$complete_primary_optics & optical$compartment == "surface") == 59,
  sum(optical$complete_primary_optics & optical$compartment == "hyporheic") == 56)
for (variable in c("npoc_mg_l", "peak_a_ru", "peak_c_ru", "peak_t_ru", "suva254")) {
  stopifnot(all(optical[[variable]][is.finite(optical[[variable]])] > 0))
}

# 3. Optical covariance and a descriptive PCA ----------------------------------
# Surface-water centering/scaling defines one reference for both compartments.
# The axes are not assigned terrestrial/algal names in advance.
surface <- optical %>% filter(compartment == "surface", complete_primary_optics)
surface_matrix <- as.matrix(surface[optical_variables])
reference_pca <- prcomp(surface_matrix, center = TRUE, scale. = TRUE)
for (axis in seq_len(ncol(reference_pca$rotation))) {
  if (reference_pca$rotation["suva254", axis] < 0) {
    reference_pca$rotation[, axis] <- -reference_pca$rotation[, axis]
    reference_pca$x[, axis] <- -reference_pca$x[, axis]
  }
}
optical$optical_pc1 <- NA_real_
optical$optical_pc2 <- NA_real_
complete <- optical$complete_primary_optics
projected <- predict(reference_pca, as.matrix(optical[complete, optical_variables]))
optical$optical_pc1[complete] <- projected[, 1]
optical$optical_pc2[complete] <- projected[, 2]
for (variable in optical_variables) {
  optical[[paste0(variable, "_z")]] <-
    (optical[[variable]] - reference_pca$center[variable]) / reference_pca$scale[variable]
}
loadings <- as_tibble(reference_pca$rotation, rownames = "variable")
write_csv(loadings, file.path(results_dir, "surface_optical_pca_loadings.csv"))
write_csv(tibble(variable = optical_variables, center = reference_pca$center,
  standard_deviation = reference_pca$scale), file.path(results_dir, "surface_optical_pca_scaling.csv"))
percent <- 100 * reference_pca$sdev^2 / sum(reference_pca$sdev^2)
write_csv(tibble(axis = seq_along(percent), percent_variance = percent),
  file.path(results_dir, "surface_optical_pca_variance.csv"))

# Rank-based and leave-one-site-out views assess sensitivity to unusual indices.
rank_matrix <- surface_matrix
for (j in seq_len(ncol(rank_matrix))) {
  rank_matrix[, j] <- rank(rank_matrix[, j], ties.method = "average")
}
rank_pca <- prcomp(rank_matrix, scale. = TRUE)
if (rank_pca$rotation["suva254", 1] < 0) {
  rank_pca$rotation[, 1] <- -rank_pca$rotation[, 1]
  rank_pca$x[, 1] <- -rank_pca$x[, 1]
}
write_csv(as_tibble(rank_pca$rotation, rownames = "variable"),
  file.path(results_dir, "rank_optical_pca_loadings.csv"))
leave_one_out <- list()
for (i in seq_len(nrow(surface))) {
  fit <- prcomp(surface_matrix[-i, , drop = FALSE], scale. = TRUE)
  if (fit$rotation["suva254", 1] < 0) {
    fit$rotation[, 1] <- -fit$rotation[, 1]
    fit$x[, 1] <- -fit$x[, 1]
  }
  leave_one_out[[i]] <- tibble(omitted_site = surface$site_code[i],
    pc1_percent = 100 * fit$sdev[1]^2 / sum(fit$sdev^2),
    pc1_score_correlation = cor(fit$x[, 1], reference_pca$x[-i, 1]),
    suva_loading = fit$rotation["suva254", 1],
    fi_loading = fit$rotation["fluorescence_index", 1],
    log_t_c_loading = fit$rotation["log_t_c", 1])
}
write_csv(bind_rows(leave_one_out), file.path(results_dir, "surface_optical_pca_leave_one_out.csv"))
write_csv(tibble(rank_pc1_spearman_correlation = cor(rank_pca$x[, 1],
  reference_pca$x[, 1], method = "spearman")), file.path(results_dir, "pca_rank_sensitivity.csv"))
write_csv(optical, file.path(derived_dir, "dom_optical_site_profiles.csv"))

# 4. Test coherence and river-position associations -----------------------------
covariance_rows <- list()
gradient_rows <- list()
row_id <- 0
secondary_variables <- c("log_t_a", "peak_a_per_npoc", "peak_c_per_npoc", "peak_t_per_npoc",
  "source_spectral_slope_ratio", "source_slope_ratio_275_295_to_300_350", "optical_pc1")
for (compartment in water_compartments) {
  sample <- optical %>% filter(.data$compartment == .env$compartment)
  stopifnot(!anyDuplicated(sample$site_code))
  pairs <- combn(c(optical_variables, "log_t_a", "peak_a_ru", "peak_c_ru", "peak_t_ru"), 2)
  for (j in seq_len(ncol(pairs))) {
    keep <- is.finite(sample[[pairs[1, j]]]) & is.finite(sample[[pairs[2, j]]])
    test <- cor.test(sample[[pairs[1, j]]][keep], sample[[pairs[2, j]]][keep],
      method = "spearman", exact = FALSE)
    row_id <- row_id + 1
    covariance_rows[[row_id]] <- tibble(compartment, variable_1 = pairs[1, j],
      variable_2 = pairs[2, j], n = sum(keep), rho = unname(test$estimate), p_value = test$p.value)
  }
  for (variable in c(optical_variables, secondary_variables)) {
    keep <- is.finite(sample[[variable]]) & is.finite(sample$log_area)
    if (sum(keep) < 3) {
      next
    }
    test <- cor.test(sample[[variable]][keep], sample$log_area[keep],
      method = "spearman", exact = FALSE)
    gradient_rows[[paste(compartment, variable)]] <- tibble(compartment, variable,
      test_family = ifelse(variable %in% optical_variables, "primary", "secondary"),
      n = sum(keep), rho = unname(test$estimate), p_value = test$p.value)
  }
}
write_csv(bind_rows(covariance_rows) %>% group_by(compartment) %>%
  mutate(p_fdr = p.adjust(p_value, "BH")) %>% ungroup(),
  file.path(results_dir, "optical_covariance_correlations.csv"))
write_csv(bind_rows(gradient_rows) %>% group_by(test_family) %>%
  mutate(p_fdr = p.adjust(p_value, "BH")) %>% ungroup(),
  file.path(results_dir, "optical_drainage_correlations.csv"))
write_csv(optical %>% group_by(compartment) %>%
  summarize(n = sum(complete_primary_optics), across(all_of(c(optical_variables, "t_c_ratio")),
    list(min = ~min(.x, na.rm = TRUE), median = ~median(.x, na.rm = TRUE),
      max = ~max(.x, na.rm = TRUE))), .groups = "drop"),
  file.path(results_dir, "optical_ranges.csv"))

# Surface-hyporheic contrasts use paired sites, not two independent site pools.
paired_optics <- optical %>% filter(complete_primary_optics) %>%
  select(site_code, compartment, all_of(c(optical_variables, "optical_pc1", "optical_pc2"))) %>%
  pivot_wider(names_from = compartment,
    values_from = all_of(c(optical_variables, "optical_pc1", "optical_pc2"))) %>%
  filter(is.finite(suva254_surface), is.finite(suva254_hyporheic))
paired_tests <- list()
for (variable in optical_variables) {
  differences <- paired_optics[[paste0(variable, "_hyporheic")]] -
    paired_optics[[paste0(variable, "_surface")]]
  test <- wilcox.test(differences, mu = 0, exact = FALSE)
  paired_tests[[variable]] <- tibble(variable, n_pairs = nrow(paired_optics),
    median_hyporheic_minus_surface = median(differences), p_value = test$p.value)
}
write_csv(bind_rows(paired_tests) %>% mutate(p_fdr = p.adjust(p_value, "BH")),
  file.path(results_dir, "paired_water_compartment_tests.csv"))
write_csv(paired_optics, file.path(results_dir, "paired_optical_site_coordinates.csv"))

# 5. Cross-method comparisons at sites with sediment molecular data -------------
formula <- read_csv(here("results", "fticr_2016", "river_composition",
  "site_formula_composition.csv"), col_types = cols(site_code = col_character()))
binary <- read_csv(here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  col_types = cols(site_code = col_character()))
pa <- as.matrix(binary[, -1])
rownames(pa) <- binary$site_code
stopifnot(all(pa %in% c(0, 1)), identical(binary$site_code, formula$site_code))
jaccard <- vegdist(pa, "jaccard", binary = TRUE)
pair_components <- read_csv(here("results", "fticr_2016", "river_composition",
  "molecular_pairwise_components.csv"), col_types = cols(site_1 = col_character(), site_2 = col_character()))
replacement <- matrix(0, nrow(pa), nrow(pa), dimnames = list(rownames(pa), rownames(pa)))
row_index <- match(pair_components$site_1, rownames(pa))
column_index <- match(pair_components$site_2, rownames(pa))
replacement[cbind(row_index, column_index)] <- pair_components$replacement
replacement[cbind(column_index, row_index)] <- pair_components$replacement
equal_features <- read_csv(here("results", "fticr_2016", "river_composition",
  "mean_subsampled_jaccard.csv"), col_types = cols(site_code = col_character()))
stopifnot(identical(equal_features$site_code, rownames(pa)),
  identical(names(equal_features)[-1], rownames(pa)))
equal_distance <- as.matrix(equal_features[, -1])
rownames(equal_distance) <- equal_features$site_code
distance_views <- list(jaccard = as.matrix(jaccard), replacement = replacement,
  equal_feature_mean = equal_distance)
matched_optics <- optical %>%
  inner_join(formula %>% select(site_code, all_of(formula_traits), detected_primary_features,
    pcoa1, pcoa2), by = "site_code")
write_csv(matched_optics, file.path(results_dir, "matched_optical_molecular_profiles.csv"))
chemical_rows <- list()
distance_rows <- list()
joint_rows <- list()
for (compartment in water_compartments) {
  matched <- matched_optics %>% filter(.data$compartment == .env$compartment, complete_primary_optics)
  expected_sites <- c(surface = 44, hyporheic = 42)
  stopifnot(!anyDuplicated(matched$site_code), nrow(matched) == expected_sites[compartment])
  indices <- match(matched$site_code, rownames(pa))
  molecular_distance <- as.dist(as.matrix(jaccard)[indices, indices])
  z_columns <- paste0(optical_variables, "_z")
  water_matrix <- as.matrix(matched[z_columns])
  rownames(water_matrix) <- matched$site_code
  optical_distance <- dist(water_matrix)
  stopifnot(identical(attr(molecular_distance, "Labels"), attr(optical_distance, "Labels")))
  set.seed(permutation_seed + match(compartment, water_compartments))
  permutations <- shuffleSet(nrow(matched), nset = n_permutations)
  for (view in names(distance_views)) {
    molecular_view <- as.dist(distance_views[[view]][indices, indices])
    test <- mantel(molecular_view, optical_distance, method = "spearman", permutations = permutations)
    distance_rows[[paste(compartment, view)]] <- tibble(compartment, molecular_view = view,
      n_sites = nrow(matched), rho = unname(test$statistic), p_value = test$signif)
  }
  # Condition on detection count and drainage before testing the three optical variables jointly.
  fit <- dbrda(molecular_distance ~ suva254_z + fluorescence_index_z + log_t_c_z +
    Condition(log_area + detected_primary_features), data = matched, add = "lingoes")
  test <- anova(fit, permutations = permutations)
  joint_rows[[compartment]] <- tibble(compartment, n_sites = nrow(matched), optical_df = test$Df[1],
    total_inertia = fit$tot.chi, conditioned_inertia = fit$pCCA$tot.chi,
    optical_constrained_inertia = fit$CCA$tot.chi,
    fraction_total_inertia = fit$CCA$tot.chi / fit$tot.chi,
    fraction_after_conditioning = fit$CCA$tot.chi / (fit$tot.chi - fit$pCCA$tot.chi),
    pseudo_f = test$F[1], p_value = test$`Pr(>F)`[1])
  for (optical_variable in optical_variables) {
    for (trait in formula_traits) {
      test <- cor.test(matched[[optical_variable]], matched[[trait]], method = "spearman", exact = FALSE)
      chemical_rows[[paste(compartment, optical_variable, trait)]] <- tibble(compartment,
        optical_variable, formula_trait = trait, n_sites = nrow(matched),
        rho = unname(test$estimate), p_value = test$p.value)
    }
  }
}
write_csv(bind_rows(distance_rows) %>% mutate(p_fdr = p.adjust(p_value, "BH")),
  file.path(results_dir, "optical_molecular_distance_mantel.csv"))
write_csv(bind_rows(joint_rows) %>% mutate(p_fdr = p.adjust(p_value, "BH")),
  file.path(results_dir, "optical_molecular_conditional_dbrda.csv"))
write_csv(bind_rows(chemical_rows) %>% group_by(compartment) %>%
  mutate(p_fdr = p.adjust(p_value, "BH")) %>% ungroup(),
  file.path(results_dir, "optical_molecular_trait_correlations.csv"))

# 6. Save seven single-panel figures -------------------------------------------
surface <- optical %>% filter(compartment == "surface", complete_primary_optics)
arrows <- loadings %>% mutate(axis1 = PC1 * 2, axis2 = PC2 * 2,
  label = unname(variable_labels[variable]), label_vjust = ifelse(axis2 > 0, -0.7, 1.6),
  label_hjust = ifelse(axis1 < 0, 1.1, 0.5), label_x = axis1 + ifelse(axis1 < 0, -0.8, 0))
pca_plot <- ggplot(surface, aes(optical_pc1, optical_pc2)) +
  geom_point(aes(color = log_area), size = 2.7) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = axis1, yend = axis2),
    inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.1, "inches"))) +
  geom_text(data = arrows, aes(label_x, axis2, label = label, vjust = label_vjust, hjust = label_hjust),
    inherit.aes = FALSE, size = 3.4) +
  geom_text(data = surface %>% filter(site_code == "CC-1"), aes(label = site_code),
    vjust = 1.5, size = 3.4) +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") + coord_equal() +
  labs(x = sprintf("Optical PC1 (%.1f%%)", percent[1]),
    y = sprintf("Optical PC2 (%.1f%%)", percent[2])) + plot_theme
ggsave(file.path(figures_dir, "2016_dom_surface_optical_pca.pdf"),
  pca_plot, width = 7, height = 5.5, bg = "white")

contrast_plot <- ggplot(surface, aes(suva254, t_c_ratio)) +
  geom_point(aes(color = log_area), size = 2.7) +
  geom_text(data = surface %>% filter(site_code == "CC-1"), aes(label = site_code),
    hjust = -0.2, size = 3.4) + scale_y_log10() +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  labs(x = "Surface-water SUVA254", y = "Surface-water EEM peak T:C") + plot_theme
ggsave(file.path(figures_dir, "2016_dom_surface_suva_peak_ratio.pdf"),
  contrast_plot, width = 6, height = 4.8, bg = "white")

plot_variables <- c("suva254", "fluorescence_index", "t_c_ratio")
plot_labels <- c("SUVA254", "Fluorescence index", "EEM peak T:C")
for (j in seq_along(plot_variables)) {
  variable <- plot_variables[j]
  drainage_plot <- ggplot(optical %>% filter(complete_primary_optics),
    aes(site_drainage_area_ha, .data[[variable]], color = compartment, shape = compartment)) +
    geom_point(size = 2.6, alpha = 0.8) + scale_x_log10() +
    scale_color_manual(name = NULL, values = compartment_colors, labels = compartment_labels) +
    scale_shape_manual(name = NULL, values = c(surface = 16, hyporheic = 2), labels = compartment_labels) +
    labs(x = "Drainage area (ha)", y = plot_labels[j]) + plot_theme +
    theme(legend.position = "bottom")
  if (variable == "t_c_ratio") {
    drainage_plot <- drainage_plot + scale_y_log10()
  }
  ggsave(file.path(figures_dir, paste0("2016_dom_", variable, "_drainage.pdf")),
    drainage_plot, width = 6, height = 4.8, bg = "white")
}

chemical_plot <- ggplot(matched_optics %>% filter(complete_primary_optics),
  aes(t_c_ratio, protein_like_fraction, color = compartment, shape = compartment)) +
  geom_point(size = 2.6, alpha = 0.8) + scale_x_log10() +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_manual(name = NULL, values = compartment_colors, labels = compartment_labels) +
  scale_shape_manual(name = NULL, values = c(surface = 16, hyporheic = 2), labels = compartment_labels) +
  labs(x = "Water-column EEM peak T:C", y = "Sediment protein-like feature fraction") +
  plot_theme + theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "2016_dom_optical_molecular_protein_comparison.pdf"),
  chemical_plot, width = 6, height = 4.8, bg = "white")

paired_plot <- ggplot(paired_optics) +
  geom_segment(aes(x = optical_pc1_surface, y = optical_pc2_surface,
    xend = optical_pc1_hyporheic, yend = optical_pc2_hyporheic), color = "gray70") +
  geom_point(aes(optical_pc1_surface, optical_pc2_surface, color = "surface", shape = "surface"), size = 2.5) +
  geom_point(aes(optical_pc1_hyporheic, optical_pc2_hyporheic, color = "hyporheic", shape = "hyporheic"), size = 2.5) +
  scale_color_manual(name = NULL, values = compartment_colors, labels = compartment_labels) +
  scale_shape_manual(name = NULL, values = c(surface = 16, hyporheic = 2), labels = compartment_labels) +
  coord_equal() + labs(x = "Optical PC1", y = "Optical PC2") + plot_theme +
  theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "2016_dom_paired_water_optical_positions.pdf"),
  paired_plot, width = 6, height = 5.5, bg = "white")
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
message("DOM optical integration completed: 59 surface sites, 56 hyporheic sites; 44/42 molecular matches.")
