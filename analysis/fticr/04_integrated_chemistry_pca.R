# Summarize complementary chemistry blocks and test microbial composition.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "integrated_chemistry")
derived_dir <- here("data", "derived", "fticr_2016")
figures_dir <- here("figures")
n_permutations <- 9999
permutation_seed <- 20261006
retained_axes <- 2
variable_design <- tribble(
  ~variable, ~block, ~label,
  "mean_detected_h_c", "fticr", "FT-ICR mean H:C",
  "mean_detected_o_c", "fticr", "FT-ICR mean O:C",
  "nitrogen_fraction", "fticr", "FT-ICR N-bearing fraction",
  "protein_like_fraction", "fticr", "FT-ICR protein-like fraction",
  "tannin_like_fraction", "fticr", "FT-ICR tannin-like fraction",
  "surface_suva254", "surface", "Surface SUVA254",
  "surface_fluorescence_index", "surface", "Surface fluorescence index",
  "surface_log_t_c", "surface", "Surface log10 peak T:C",
  "hyporheic_suva254", "hyporheic", "Hyporheic SUVA254",
  "hyporheic_fluorescence_index", "hyporheic", "Hyporheic fluorescence index",
  "hyporheic_log_t_c", "hyporheic", "Hyporheic log10 peak T:C")
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank())

# 2. Construct and audit the common site table ----------------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
formula <- read_csv(here("results", "fticr_2016", "river_composition",
  "site_formula_composition.csv"), col_types = cols(site_code = col_character()))
optical <- read_csv(here("data", "derived", "fticr_2016", "dom_optical_site_profiles.csv"),
  col_types = cols(site_code = col_character()))
water <- optical %>% select(site_code, compartment, suva254, fluorescence_index, log_t_c) %>%
  pivot_wider(names_from = compartment, values_from = c(suva254, fluorescence_index, log_t_c),
    names_glue = "{compartment}_{.value}")
site_table <- formula %>%
  select(site_code, sediment_sample_id, site_drainage_area_ha, site_stream_order,
    detected_primary_features, all_of(variable_design$variable[variable_design$block == "fticr"])) %>%
  left_join(water, by = "site_code")
stopifnot(nrow(site_table) == 44, !anyDuplicated(site_table$site_code))
audit <- site_table %>% select(site_code)
audit$complete_integrated_chemistry <- complete.cases(site_table[variable_design$variable])
audit$missing_variables <- ""
for (i in seq_len(nrow(site_table))) {
  missing <- variable_design$variable[!is.finite(as.numeric(site_table[i, variable_design$variable]))]
  audit$missing_variables[i] <- paste(missing, collapse = ";")
}
write_csv(audit, file.path(results_dir, "integrated_site_inclusion_audit.csv"))
site <- site_table[audit$complete_integrated_chemistry, ]
stopifnot(nrow(site) == 42, all(is.finite(as.matrix(site[variable_design$variable]))))
site$log_area <- log10(site$site_drainage_area_ha)
write_csv(site, file.path(derived_dir, "integrated_chemistry_42_site_inputs.csv"))

# 3. PCA within a common site set, with equal total variance per method block ----
# Standardization occurs before weighting. Each block contributes variance one,
# so the five FT-ICR properties do not outweigh a three-variable optical block.
variable_design <- variable_design %>% group_by(block) %>%
  mutate(block_variables = n(), block_weight = 1 / sqrt(block_variables)) %>% ungroup()
raw_matrix <- as.matrix(site[variable_design$variable])
rownames(raw_matrix) <- site$site_code
rank_matrix <- raw_matrix
for (j in seq_len(ncol(rank_matrix))) {
  rank_matrix[, j] <- rank(rank_matrix[, j], ties.method = "average")
}
raw_z <- scale(raw_matrix)
rank_z <- scale(rank_matrix)
weighted_z <- sweep(raw_z, 2, variable_design$block_weight, "*")
weighted_rank_z <- sweep(rank_z, 2, variable_design$block_weight, "*")
block_variance <- tibble(variable = colnames(weighted_z), variance = apply(weighted_z, 2, var)) %>%
  left_join(variable_design, by = "variable") %>% group_by(block) %>%
  summarize(total_variance = sum(variance), .groups = "drop")
stopifnot(all(abs(block_variance$total_variance - 1) < 1e-12))
write_csv(block_variance, file.path(results_dir, "block_variance_audit.csv"))
write_csv(variable_design %>% mutate(center = as.numeric(attr(raw_z, "scaled:center")),
  standard_deviation = as.numeric(attr(raw_z, "scaled:scale"))),
  file.path(results_dir, "integrated_variable_scaling.csv"))
fticr_columns <- which(variable_design$block == "fticr")
matrix_views <- list(fticr_only = raw_z[, fticr_columns, drop = FALSE],
  integrated = weighted_z, integrated_rank = weighted_rank_z)
scores_list <- list()
loading_rows <- list()
variance_rows <- list()
for (view in names(matrix_views)) {
  input_matrix <- matrix_views[[view]]
  fit <- prcomp(input_matrix, center = FALSE, scale. = FALSE)
  # Signs orient higher H:C toward positive scores; they carry no source attribution.
  for (axis in seq_len(ncol(fit$rotation))) {
    if (fit$rotation["mean_detected_h_c", axis] < 0) {
      fit$rotation[, axis] <- -fit$rotation[, axis]
      fit$x[, axis] <- -fit$x[, axis]
    }
  }
  scores_list[[view]] <- fit$x[, seq_len(retained_axes), drop = FALSE]
  correlations <- cor(raw_matrix[, colnames(input_matrix), drop = FALSE], fit$x)
  if (view == "integrated_rank") {
    correlations <- cor(rank_matrix[, colnames(input_matrix), drop = FALSE], fit$x)
  }
  loading_rows[[view]] <- tibble(view, variable = rownames(fit$rotation),
    pc1_loading = fit$rotation[, 1], pc2_loading = fit$rotation[, 2],
    pc1_variable_correlation = correlations[, 1], pc2_variable_correlation = correlations[, 2])
  fraction <- fit$sdev^2 / sum(fit$sdev^2)
  variance_rows[[view]] <- tibble(view, axis = seq_along(fraction), percent_variance = 100 * fraction)
}
write_csv(bind_rows(loading_rows), file.path(results_dir, "chemistry_pca_loadings.csv"))
variance_summary <- bind_rows(variance_rows)
write_csv(variance_summary, file.path(results_dir, "chemistry_pca_variance.csv"))
for (view in names(scores_list)) {
  site[[paste0(view, "_pc1")]] <- scores_list[[view]][, 1]
  site[[paste0(view, "_pc2")]] <- scores_list[[view]][, 2]
}
write_csv(site, file.path(results_dir, "chemistry_site_scores.csv"))
write_csv(tibble(comparison = "Raw versus rank integrated PC1",
  spearman_rho = cor(site$integrated_pc1, site$integrated_rank_pc1, method = "spearman")),
  file.path(results_dir, "chemistry_pca_rank_sensitivity.csv"))

# 4. Explain microbial composition using two preselected chemistry axes ----------
# Axes are learned without reference to microbial outcomes; individual ASV tests
# and selection among later chemistry axes are outside this analysis.
metadata <- read_csv(here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  col_types = cols(site_code = col_character())) %>%
  filter(habitat == "sediment", included_10k) %>% select(site_code, sample_id, plate)
paired <- site %>% inner_join(metadata, by = "site_code")
stopifnot(nrow(paired) == 33, !anyDuplicated(paired$sample_id))
counts <- read_csv(here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  show_col_types = FALSE)
microbial_counts <- t(as.matrix(counts[, paired$sample_id]))
rownames(microbial_counts) <- paired$site_code
colnames(microbial_counts) <- counts$asv_id
stopifnot(all(rowSums(microbial_counts) == 10000))
microbial_counts <- microbial_counts[, colSums(microbial_counts) > 0, drop = FALSE]
microbial_hellinger <- decostand(microbial_counts, "hellinger")
# All microbial PCA coordinates preserve the response geometry exactly while
# avoiding permutations over thousands of redundant response columns.
microbial_pca <- prcomp(microbial_hellinger, center = TRUE, scale. = FALSE)
microbial_coordinates <- microbial_pca$x
stopifnot(max(abs(as.vector(dist(microbial_coordinates)) -
  as.vector(dist(microbial_hellinger)))) < 1e-10)
set.seed(permutation_seed)
permutations <- shuffleSet(nrow(paired), nset = n_permutations)
model_design <- tibble(view = names(scores_list), retained_axes = retained_axes,
  n_sites = nrow(paired), permutation_seed = permutation_seed, n_permutations = n_permutations)
write_csv(model_design, file.path(results_dir, "microbial_model_design.csv"))
global_rows <- list()
axis_rows <- list()
validation_rows <- list()
for (view in model_design$view) {
  paired$chemistry_pc1 <- paired[[paste0(view, "_pc1")]]
  paired$chemistry_pc2 <- paired[[paste0(view, "_pc2")]]
  full_fit <- rda(microbial_hellinger ~ chemistry_pc1 + chemistry_pc2 +
    Condition(log_area + detected_primary_features), data = paired)
  fit <- rda(microbial_coordinates ~ chemistry_pc1 + chemistry_pc2 +
    Condition(log_area + detected_primary_features), data = paired)
  inertia_difference <- max(abs(c(fit$tot.chi - full_fit$tot.chi,
    fit$CCA$tot.chi - full_fit$CCA$tot.chi, fit$pCCA$tot.chi - full_fit$pCCA$tot.chi)))
  adjusted_difference <- abs(RsquareAdj(fit)$adj.r.squared - RsquareAdj(full_fit)$adj.r.squared)
  stopifnot(inertia_difference < 1e-10, adjusted_difference < 1e-10)
  validation_rows[[view]] <- tibble(view, max_inertia_difference = inertia_difference,
    adjusted_r_squared_difference = adjusted_difference, full_response_columns = ncol(microbial_hellinger),
    coordinate_response_columns = ncol(microbial_coordinates))
  started <- proc.time()[["elapsed"]]
  test <- anova(fit, permutations = permutations)
  adjusted <- RsquareAdj(fit)
  global_rows[[view]] <- tibble(view, n_sites = nrow(paired), chemistry_df = test$Df[1],
    fraction_total_variance = fit$CCA$tot.chi / fit$tot.chi,
    adjusted_unique_fraction = adjusted$adj.r.squared,
    fraction_after_conditioning = fit$CCA$tot.chi / (fit$tot.chi - fit$pCCA$tot.chi),
    pseudo_f = test$F[1], p_value = test$`Pr(>F)`[1])
  marginal <- anova(fit, by = "margin", permutations = permutations)
  message(view, " permutation tests completed in ",
    round(proc.time()[["elapsed"]] - started, 2), " seconds.")
  axis_rows[[view]] <- tibble(view, term = rownames(marginal), df = marginal$Df,
    variance = marginal$Variance, pseudo_f = marginal$F, p_value = marginal$`Pr(>F)`)
}
write_csv(bind_rows(validation_rows), file.path(results_dir, "microbial_coordinate_equivalence_checks.csv"))
write_csv(bind_rows(global_rows) %>% mutate(p_fdr = p.adjust(p_value, "BH")),
  file.path(results_dir, "microbial_chemistry_partial_rda.csv"))
write_csv(bind_rows(axis_rows) %>% mutate(p_fdr = p.adjust(p_value, "BH")),
  file.path(results_dir, "microbial_chemistry_axis_tests.csv"))
# Keep only the named fixed scores in the export, rather than loop working columns.
paired <- paired %>% select(-chemistry_pc1, -chemistry_pc2)
write_csv(paired, file.path(results_dir, "paired_microbial_chemistry_site_scores.csv"))

# 5. Site influence on the FT-ICR property association ---------------------------
# Scores stay fixed because they were learned without microbial outcomes.
# These fits assess effect-size stability; they are not independent replicate tests.
paired$chemistry_pc1 <- paired$fticr_only_pc1
paired$chemistry_pc2 <- paired$fticr_only_pc2
influence_rows <- list()
for (i in seq_len(nrow(paired))) {
  fit <- tryCatch(rda(microbial_coordinates[-i, , drop = FALSE] ~ chemistry_pc1 + chemistry_pc2 +
    Condition(log_area + detected_primary_features), data = paired[-i, ]),
    error = function(e) e)
  if (inherits(fit, "error")) {
    influence_rows[[i]] <- tibble(omitted_site = paired$site_code[i], completed = FALSE,
      fraction_total_variance = NA_real_, adjusted_unique_fraction = NA_real_,
      error_message = conditionMessage(fit))
  } else {
    influence_rows[[i]] <- tibble(omitted_site = paired$site_code[i], completed = TRUE,
      fraction_total_variance = fit$CCA$tot.chi / fit$tot.chi,
      adjusted_unique_fraction = RsquareAdj(fit)$adj.r.squared, error_message = NA_character_)
  }
}
write_csv(bind_rows(influence_rows), file.path(results_dir, "fticr_microbe_leave_one_site_out.csv"))
paired <- paired %>% select(-chemistry_pc1, -chemistry_pc2)

# 6. Save chemistry ordination, loadings, and an independent microbial display ----
fticr_percent <- variance_summary %>% filter(view == "fticr_only", axis <= retained_axes)
fticr_plot <- ggplot(site, aes(fticr_only_pc1, fticr_only_pc2, color = log_area)) +
  geom_point(size = 2.8) + scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  coord_equal() + labs(x = sprintf("FT-ICR property PC1 (%.1f%%)", fticr_percent$percent_variance[1]),
    y = sprintf("FT-ICR property PC2 (%.1f%%)", fticr_percent$percent_variance[2])) + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_property_pca.pdf"),
  fticr_plot, width = 7, height = 5, bg = "white")
axis_percent <- variance_summary %>% filter(view == "integrated", axis <= retained_axes)
chemistry_plot <- ggplot(site, aes(integrated_pc1, integrated_pc2, color = log_area)) +
  geom_point(size = 2.8) + scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  coord_equal() + labs(x = sprintf("Integrated chemistry PC1 (%.1f%%)", axis_percent$percent_variance[1]),
    y = sprintf("Integrated chemistry PC2 (%.1f%%)", axis_percent$percent_variance[2])) + plot_theme
ggsave(file.path(figures_dir, "2016_integrated_chemistry_pca.pdf"),
  chemistry_plot, width = 7, height = 5, bg = "white")
loading_table <- bind_rows(loading_rows) %>% filter(view == "integrated") %>%
  select(variable, pc1_variable_correlation, pc2_variable_correlation) %>%
  pivot_longer(-variable, names_to = "axis", values_to = "correlation") %>%
  left_join(variable_design, by = "variable") %>%
  mutate(label = factor(label, levels = rev(variable_design$label)),
    axis = factor(axis, levels = c("pc1_variable_correlation", "pc2_variable_correlation"),
      labels = c("PC1", "PC2")))
loadings_plot <- ggplot(loading_table, aes(axis, label, fill = correlation)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", correlation),
    color = abs(correlation) > 0.55), size = 3.4, show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1),
    name = "Variable-axis\ncorrelation") + labs(x = NULL, y = NULL) + plot_theme
ggsave(file.path(figures_dir, "2016_integrated_chemistry_loadings.pdf"),
  loadings_plot, width = 7, height = 6, bg = "white")
microbial_percent <- 100 * microbial_pca$sdev^2 / sum(microbial_pca$sdev^2)
microbial_plot_data <- paired %>% mutate(microbial_pc1 = microbial_pca$x[, 1],
  microbial_pc2 = microbial_pca$x[, 2])
write_csv(microbial_plot_data, file.path(results_dir, "microbial_pca_chemistry_overlay.csv"))
microbial_plot <- ggplot(microbial_plot_data,
  aes(microbial_pc1, microbial_pc2, color = integrated_pc1)) +
  geom_point(size = 2.8) + scale_color_viridis_c(name = "Integrated\nchemistry PC1") +
  coord_equal() + labs(x = sprintf("Microbial PC1 (%.1f%%)", microbial_percent[1]),
    y = sprintf("Microbial PC2 (%.1f%%)", microbial_percent[2])) + plot_theme
ggsave(file.path(figures_dir, "2016_microbial_pca_integrated_chemistry.pdf"),
  microbial_plot, width = 7, height = 5, bg = "white")
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
message("Integrated chemistry PCA completed: 42 chemistry sites, 33 microbial sites, three model views.")
