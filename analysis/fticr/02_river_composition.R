# Explore molecular incidence and chemical signatures along the river network.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
library(adespatial)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "river_composition")
figures_dir <- here("figures")
n_permutations <- 9999
n_subsamples <- 1000
subsample_seed <- 20261003
permutation_seed <- 20261004
ordination_dimensions <- 2
class_order <- c("Lipid", "Protein", "AminoSugar", "Carb", "Lignin", "Tannin", "ConHC", "UnsatHC", "Other")
class_labels <- c("Lipid-like", "Protein-like", "Amino-sugar-like", "Carbohydrate-like",
  "Lignin-like", "Tannin-like", "Condensed-hydrocarbon-like", "Unsaturated-hydrocarbon-like", "Other")
class_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#882255", "#AA4499")
names(class_colors) <- class_labels
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank())

# 2. Align incidence, formula metadata, and environment -------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
pa_table <- read_csv(here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  col_types = cols(site_code = col_character()))
pa <- as.matrix(pa_table[, -1])
rownames(pa) <- pa_table$site_code
properties <- read_csv(here("results", "fticr_2016", "tables", "primary_peak_properties.csv"),
  show_col_types = FALSE)
properties <- properties[match(colnames(pa), properties$peak_id), ]
site <- read_csv(here("results", "fticr_2016", "tables", "jaccard_pcoa_scores.csv"),
  col_types = cols(site_code = col_character()))
site <- site[match(rownames(pa), site$site_code), ]
stopifnot(identical(site$site_code, rownames(pa)),
  identical(properties$peak_id, colnames(pa)), all(pa %in% c(0, 1)),
  !anyDuplicated(site$site_code), all(rowSums(pa) > 0),
  all(properties$molecular_class_source %in% class_order))
site$log_area <- log10(site$site_drainage_area_ha)
site$nitrogen_fraction <- as.vector(pa %*% (properties$nitrogen_count > 0)) / rowSums(pa)

# Retain source classes without changing their unknown original classification rules.
class_counts <- matrix(0, nrow(pa), length(class_order),
  dimnames = list(rownames(pa), class_order))
for (j in seq_along(class_order)) {
  class_counts[, j] <- rowSums(pa[, properties$molecular_class_source == class_order[j], drop = FALSE])
}
class_fractions <- class_counts / rowSums(pa)
stopifnot(all(abs(rowSums(class_fractions) - 1) < 1e-12))
site$lignin_like_fraction <- class_fractions[, "Lignin"]
site$tannin_like_fraction <- class_fractions[, "Tannin"]
site$protein_like_fraction <- class_fractions[, "Protein"]
site$lipid_like_fraction <- class_fractions[, "Lipid"]
site$high_h_c_fraction <- as.vector(pa %*% (properties$hydrogen_carbon_ratio >= 1.5)) / rowSums(pa)
write_csv(site, file.path(results_dir, "site_formula_composition.csv"))
class_table <- bind_cols(tibble(site_code = rownames(pa)), as_tibble(class_fractions)) %>%
  pivot_longer(-site_code, names_to = "source_class", values_to = "fraction_detected") %>%
  left_join(site %>% select(site_code, site_drainage_area_ha, site_stream_order), by = "site_code") %>%
  mutate(class_label = factor(source_class, levels = class_order, labels = class_labels))
write_csv(class_table, file.path(results_dir, "site_class_fractions.csv"))

# 3. River gradients and molecular ordination interpretation --------------------
# Small exploratory summaries, with FDR correction within each predictor family.
trait_names <- c("mean_detected_h_c", "mean_detected_o_c", "nitrogen_fraction",
  "lignin_like_fraction", "tannin_like_fraction", "protein_like_fraction", "lipid_like_fraction", "high_h_c_fraction")
gradient_names <- c("log_area", "site_stream_order", "pcoa1", "pcoa2",
  "surface_suva254", "surface_fluorescence_index")
correlation_rows <- list()
row_index <- 0
for (gradient in gradient_names) {
  for (trait in trait_names) {
    keep <- is.finite(site[[gradient]]) & is.finite(site[[trait]])
    test <- cor.test(site[[gradient]][keep], site[[trait]][keep], method = "spearman", exact = FALSE)
    row_index <- row_index + 1
    correlation_rows[[row_index]] <- tibble(predictor = gradient, trait = trait,
      n = sum(keep), spearman_rho = unname(test$estimate), p_value = test$p.value)
  }
}
correlations <- bind_rows(correlation_rows) %>%
  group_by(predictor) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_csv(correlations, file.path(results_dir, "formula_gradient_correlations.csv"))
ordination_diagnostics <- tibble(
  comparison = c("PCoA1 versus detection count", "PCoA1 versus log drainage area"),
  spearman_rho = c(cor(site$pcoa1, rowSums(pa), method = "spearman"),
    cor(site$pcoa1, site$log_area, method = "spearman")))
write_csv(ordination_diagnostics, file.path(results_dir, "ordination_detection_diagnostics.csv"))

# Mean H:C is a chemical signature; a source fraction is not estimated here.
simple_fit <- lm(mean_detected_h_c ~ log_area, data = site)
adjusted_fit <- lm(mean_detected_h_c ~ log_area + detected_primary_features, data = site)
models <- list(simple = simple_fit, detection_adjusted = adjusted_fit)
coefficient_rows <- list()
for (model_name in names(models)) {
  coefficients <- as.data.frame(coef(summary(models[[model_name]])))
  interval <- confint(models[[model_name]])
  coefficient_rows[[model_name]] <- tibble(model = model_name, term = rownames(coefficients),
    estimate = coefficients[, 1], standard_error = coefficients[, 2],
    p_value = coefficients[, 4], lower_95 = interval[, 1], upper_95 = interval[, 2])
}
write_csv(bind_rows(coefficient_rows), file.path(results_dir, "h_c_drainage_models.csv"))
order_summary <- site %>%
  group_by(site_stream_order) %>%
  summarize(n = n(), median_area_ha = median(site_drainage_area_ha),
    mean_h_c = mean(mean_detected_h_c), mean_o_c = mean(mean_detected_o_c),
    lignin_like_fraction = mean(lignin_like_fraction),
    protein_like_fraction = mean(protein_like_fraction),
    detection_count = mean(detected_primary_features), .groups = "drop")
write_csv(order_summary, file.path(results_dir, "stream_order_composition_summary.csv"))

# Split Jaccard dissimilarity into replacement and nestedness-resultant components.
# Neither component distinguishes biological differences from nondetection.
partition <- beta.div.comp(pa, coef = "BJ", quant = FALSE)
jaccard <- vegdist(pa, method = "jaccard", binary = TRUE)
stopifnot(max(abs(as.vector(partition$D) - as.vector(jaccard))) < 1e-12,
  max(abs(as.vector(partition$D) - as.vector(partition$repl) - as.vector(partition$rich))) < 1e-12)
pair_indices <- which(lower.tri(as.matrix(jaccard)), arr.ind = TRUE)
pair_table <- tibble(site_1 = rownames(pa)[pair_indices[, 1]],
  site_2 = rownames(pa)[pair_indices[, 2]], jaccard = as.vector(jaccard),
  replacement = as.vector(partition$repl), nestedness_resultant = as.vector(partition$rich))
write_csv(pair_table, file.path(results_dir, "molecular_pairwise_components.csv"))
write_csv(tibble(component = c("replacement", "nestedness_resultant"),
  fraction_of_summed_jaccard = c(sum(partition$repl), sum(partition$rich)) / sum(jaccard)),
  file.path(results_dir, "molecular_beta_partition_summary.csv"))

# Site-level permutations assess drainage associations; pairs are not replicates.
set.seed(permutation_seed)
permutations_44 <- shuffleSet(nrow(pa), nset = n_permutations)
distance_views <- list(jaccard = jaccard, replacement = partition$repl)
drainage_tests <- list()
for (view_name in names(distance_views)) {
  test <- adonis2(distance_views[[view_name]] ~ log_area,
    data = site, permutations = permutations_44, add = "lingoes")
  drainage_tests[[view_name]] <- tibble(view = view_name, sites = nrow(pa),
    r_squared = test$R2[1], pseudo_f = test$F[1], p_value = test$`Pr(>F)`[1])
}
write_csv(bind_rows(drainage_tests), file.path(results_dir, "molecular_drainage_permutation_tests.csv"))

# 4. Molecular-microbial comparisons at the 33 paired sites ----------------------
microbial_metadata <- read_csv(here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  col_types = cols(site_code = col_character())) %>%
  filter(habitat == "sediment", included_10k)
paired_site <- site %>%
  inner_join(microbial_metadata %>% select(site_code, sample_id, plate), by = "site_code")
stopifnot(nrow(paired_site) == 33, !anyDuplicated(paired_site$sample_id))
counts_table <- read_csv(here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  show_col_types = FALSE)
microbial_counts <- t(as.matrix(counts_table[, paired_site$sample_id]))
colnames(microbial_counts) <- counts_table$asv_id
rownames(microbial_counts) <- paired_site$site_code
stopifnot(all(rowSums(microbial_counts) == 10000))
microbial_counts <- microbial_counts[, colSums(microbial_counts) > 0, drop = FALSE]
microbial_hellinger <- decostand(microbial_counts, "hellinger")
microbial_distance <- dist(microbial_hellinger)
paired_rows <- match(paired_site$site_code, rownames(pa))
paired_pa <- pa[paired_rows, , drop = FALSE]
paired_jaccard <- vegdist(paired_pa, "jaccard", binary = TRUE)
paired_replacement <- as.dist(as.matrix(partition$repl)[paired_rows, paired_rows])
stopifnot(identical(attr(paired_jaccard, "Labels"), attr(microbial_distance, "Labels")))
set.seed(permutation_seed + 1)
permutations_33 <- shuffleSet(nrow(paired_site), nset = n_permutations)
set.seed(permutation_seed + 2)
order_permutations <- shuffleSet(nrow(paired_site), nset = n_permutations,
  control = how(blocks = factor(paired_site$site_stream_order)))
mantel_tests <- list()
paired_distances <- list(jaccard = paired_jaccard, replacement = paired_replacement)
row_index <- 0
for (view_name in names(paired_distances)) {
  for (scheme in c("unrestricted", "within_stream_order")) {
    permutation_matrix <- permutations_33
    if (scheme == "within_stream_order") {
      permutation_matrix <- order_permutations
    }
    test <- mantel(paired_distances[[view_name]], microbial_distance,
      method = "spearman", permutations = permutation_matrix)
    row_index <- row_index + 1
    mantel_tests[[row_index]] <- tibble(molecular_view = view_name, permutation_scheme = scheme,
      n_sites = nrow(paired_site), mantel_rho = unname(test$statistic), p_value = test$signif)
  }
}
mantel_summary <- bind_rows(mantel_tests) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH"))
write_csv(mantel_summary, file.path(results_dir, "paired_composition_mantel.csv"))

# Refit both two-dimensional ordinations on the identical paired subset.
molecular_pcoa <- wcmdscale(paired_jaccard, k = ordination_dimensions, eig = TRUE, add = "lingoes")
microbial_pca <- prcomp(microbial_hellinger, center = TRUE, scale. = FALSE)
microbial_scores <- microbial_pca$x[, seq_len(ordination_dimensions), drop = FALSE]
procrustes_fit <- procrustes(molecular_pcoa$points, microbial_scores,
  symmetric = TRUE, scale = FALSE)
procrustes_test <- protest(molecular_pcoa$points, microbial_scores,
  permutations = permutations_33)
write_csv(tibble(n_sites = nrow(paired_site), dimensions = ordination_dimensions,
  procrustes_correlation = procrustes_test$t0, p_value = procrustes_test$signif,
  molecular_two_axis_percent = 100 * sum(molecular_pcoa$eig[1:2]) / sum(molecular_pcoa$eig[molecular_pcoa$eig > 0]),
  microbial_two_axis_percent = 100 * sum(microbial_pca$sdev[1:2]^2) / sum(microbial_pca$sdev^2)),
  file.path(results_dir, "paired_procrustes_summary.csv"))
paired_coordinates <- paired_site %>%
  mutate(molecular_axis1 = procrustes_fit$X[, 1], molecular_axis2 = procrustes_fit$X[, 2],
    microbial_axis1 = procrustes_fit$Yrot[, 1], microbial_axis2 = procrustes_fit$Yrot[, 2],
    procrustes_residual = residuals(procrustes_fit))
write_csv(paired_coordinates, file.path(results_dir, "paired_site_coordinates.csv"))

# One chemical signature is assessed after accounting for drainage and detection count. This exploratory model does not establish microbial function.
microbial_fit <- adonis2(microbial_hellinger ~ mean_detected_h_c + log_area +
  detected_primary_features, data = paired_site,
  method = "euclidean", by = "margin", permutations = permutations_33)
write_csv(tibble(term = rownames(microbial_fit), df = microbial_fit$Df,
  sum_of_squares = microbial_fit$SumOfSqs, r_squared = microbial_fit$R2,
  pseudo_f = microbial_fit$F, p_value = microbial_fit$`Pr(>F)`),
  file.path(results_dir, "paired_microbial_h_c_model.csv"))

# 5. Equal-feature subsampling as a sensitivity analysis ------------------------
# This thins observed features uniformly; it does not simulate ions or recover
# undetected compounds. It is not used to estimate standardized molecular richness.
target_features <- min(rowSums(pa))
run_design <- tibble(run_id = seq_len(n_subsamples), seed = subsample_seed + seq_len(n_subsamples))
write_csv(run_design, file.path(results_dir, "subsampling_design.csv"))
run_results <- run_design %>%
  mutate(target_features = target_features, completed = FALSE,
    drainage_r_squared = NA_real_, paired_mantel_rho = NA_real_, error_message = NA_character_)
mean_distance <- matrix(0, nrow(pa), nrow(pa), dimnames = list(rownames(pa), rownames(pa)))
subsampled_h_c <- matrix(NA_real_, n_subsamples, nrow(pa))
run_seconds <- numeric(n_subsamples)
for (run_id in run_design$run_id) {
  started <- proc.time()[["elapsed"]]
  set.seed(run_design$seed[run_id])
  run <- tryCatch({
    subsampled <- matrix(0, nrow(pa), ncol(pa), dimnames = dimnames(pa))
    for (i in seq_len(nrow(pa))) {
      selected <- sample(which(pa[i, ] > 0), size = target_features, replace = FALSE)
      subsampled[i, selected] <- 1
    }
    stopifnot(all(rowSums(subsampled) == target_features))
    run_distance <- vegdist(subsampled, "jaccard", binary = TRUE)
    run_test <- adonis2(run_distance ~ log_area, data = site, permutations = 0)
    run_paired_distance <- as.dist(as.matrix(run_distance)[paired_rows, paired_rows])
    list(distance = as.matrix(run_distance), r_squared = run_test$R2[1],
      paired_rho = cor(as.vector(run_paired_distance), as.vector(microbial_distance), method = "spearman"),
      h_c = as.vector(subsampled %*% properties$hydrogen_carbon_ratio) / target_features)
  }, error = function(e) e)
  if (inherits(run, "error")) {
    run_results$error_message[run_id] <- conditionMessage(run)
  } else {
    run_results$completed[run_id] <- TRUE
    run_results$drainage_r_squared[run_id] <- run$r_squared
    run_results$paired_mantel_rho[run_id] <- run$paired_rho
    mean_distance <- mean_distance + run$distance
    subsampled_h_c[run_id, ] <- run$h_c
  }
  run_seconds[run_id] <- proc.time()[["elapsed"]] - started
  if (run_id == 1) {
    message("First full subsampling run: ", round(run_seconds[1], 3), " seconds.")
  }
}
run_results$elapsed_seconds <- run_seconds
write_csv(run_results, file.path(results_dir, "subsampling_runs.csv"))
stopifnot(any(run_results$completed))
mean_distance <- mean_distance / sum(run_results$completed)
mean_subsampled_distance <- as.dist(mean_distance)
write_csv(bind_cols(tibble(site_code = rownames(mean_distance)), as_tibble(mean_distance)),
  file.path(results_dir, "mean_subsampled_jaccard.csv"))
subsampling_summary <- run_results %>%
  filter(completed) %>%
  summarize(target_features = first(target_features), completed_runs = n(),
    drainage_r_squared_median = median(drainage_r_squared),
    drainage_r_squared_lower = quantile(drainage_r_squared, 0.025),
    drainage_r_squared_upper = quantile(drainage_r_squared, 0.975),
    paired_mantel_rho_median = median(paired_mantel_rho),
    paired_mantel_rho_lower = quantile(paired_mantel_rho, 0.025),
    paired_mantel_rho_upper = quantile(paired_mantel_rho, 0.975))
write_csv(subsampling_summary, file.path(results_dir, "subsampling_summary.csv"))
mean_distance_test <- adonis2(mean_subsampled_distance ~ log_area,
  data = site, permutations = permutations_44, add = "lingoes")
mean_paired_distance <- as.dist(mean_distance[paired_rows, paired_rows])
mean_mantel <- mantel(mean_paired_distance, microbial_distance,
  method = "spearman", permutations = permutations_33)
write_csv(tibble(completed_runs = sum(run_results$completed), target_features = target_features,
  mean_distance_drainage_r_squared = mean_distance_test$R2[1],
  mean_distance_drainage_p = mean_distance_test$`Pr(>F)`[1],
  mean_distance_paired_mantel_rho = unname(mean_mantel$statistic),
  mean_distance_paired_mantel_p = mean_mantel$signif),
  file.path(results_dir, "mean_subsampled_distance_tests.csv"))
subsampled_traits <- site %>% select(site_code, mean_detected_h_c)
subsampled_traits$subsample_mean_h_c <- colMeans(subsampled_h_c, na.rm = TRUE)
subsampled_traits$subsample_h_c_lower <- apply(subsampled_h_c, 2, quantile, probs = 0.025, na.rm = TRUE)
subsampled_traits$subsample_h_c_upper <- apply(subsampled_h_c, 2, quantile, probs = 0.975, na.rm = TRUE)
write_csv(subsampled_traits, file.path(results_dir, "subsampled_h_c_summary.csv"))

# 6. Save figures ---------------------------------------------------------------
site_order <- site$site_code[order(site$site_drainage_area_ha)]
class_table$site_code <- factor(class_table$site_code, levels = site_order)
class_plot <- ggplot(class_table, aes(site_code, fraction_detected, fill = class_label)) +
  geom_col(width = 0.95) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Sites ordered by increasing drainage area", y = "Fraction of detected features", fill = NULL) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3))
ggsave(file.path(figures_dir, "2016_fticr_site_class_composition.pdf"),
  class_plot, width = 11, height = 5.2, bg = "white")

# Correlation vectors cross-reference the same 44-site ordination, without changing it.
vector_fit <- envfit(as.data.frame(site[c("pcoa1", "pcoa2")]),
  as.data.frame(site[c("mean_detected_h_c", "mean_detected_o_c")]), permutations = 0)
vectors <- as.data.frame(scores(vector_fit, display = "vectors"))
names(vectors) <- c("axis1", "axis2")
vectors$label <- c("Mean H:C", "Mean O:C")
vectors$axis1 <- vectors$axis1 * 0.35
vectors$axis2 <- vectors$axis2 * 0.35
vectors$label_vjust <- ifelse(vectors$axis2 > 0, -0.6, 1.6)
overlay_plot <- ggplot(site, aes(pcoa1, pcoa2)) +
  geom_point(aes(color = log_area, size = detected_primary_features), alpha = 0.85) +
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = axis1, yend = axis2),
    inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.12, "inches"))) +
  geom_text(data = vectors, aes(axis1, axis2, label = label, vjust = label_vjust),
    inherit.aes = FALSE, size = 3.5) +
  expand_limits(y = c(min(vectors$axis2) - 0.045, max(vectors$axis2) + 0.045)) +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  scale_size_continuous(name = "Detected\nfeatures", range = c(1.5, 4.5)) +
  coord_equal() +
  labs(x = "Jaccard PCoA1", y = "Jaccard PCoA2") + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_pcoa_formula_overlay.pdf"),
  overlay_plot, width = 7, height = 5, bg = "white")

h_c_plot <- ggplot(site, aes(site_drainage_area_ha, mean_detected_h_c)) +
  geom_point(aes(color = factor(site_stream_order)), size = 2.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  scale_x_log10() + scale_color_viridis_d(name = "Stream order") +
  labs(x = "Drainage area (ha)", y = "Mean H:C of detected features") + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_h_c_drainage.pdf"),
  h_c_plot, width = 6, height = 4.8, bg = "white")

protein_plot <- ggplot(site, aes(site_drainage_area_ha, protein_like_fraction)) +
  geom_point(aes(color = factor(site_stream_order)), size = 2.8) +
  scale_x_log10() + scale_color_viridis_d(name = "Stream order") +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Drainage area (ha)", y = "Protein-like fraction of detected features") + plot_theme
ggsave(file.path(figures_dir, "2016_fticr_protein_like_drainage.pdf"),
  protein_plot, width = 6, height = 4.8, bg = "white")

procrustes_plot <- ggplot(paired_coordinates) +
  geom_segment(aes(x = molecular_axis1, y = molecular_axis2,
    xend = microbial_axis1, yend = microbial_axis2), color = "gray65") +
  geom_point(aes(molecular_axis1, molecular_axis2, color = log_area, shape = "Molecular"), size = 2.5) +
  geom_point(aes(microbial_axis1, microbial_axis2, color = log_area, shape = "Microbial"), size = 2.5) +
  scale_shape_manual(name = NULL, values = c(Molecular = 16, Microbial = 2)) +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") + coord_equal() +
  labs(x = "Procrustes axis 1", y = "Procrustes axis 2") + plot_theme
ggsave(file.path(figures_dir, "2016_paired_molecular_microbes_procrustes.pdf"),
  procrustes_plot, width = 6, height = 4.8, bg = "white")

ordination_views <- list(replacement = partition$repl, equal_feature = mean_subsampled_distance)
ordination_rows <- list()
ordination_summary_rows <- list()
for (view_name in names(ordination_views)) {
  ordination <- wcmdscale(ordination_views[[view_name]], k = 2, eig = TRUE, add = "lingoes")
  uncorrected <- wcmdscale(ordination_views[[view_name]], k = 2, eig = TRUE)
  plot_data <- site %>% mutate(axis1 = ordination$points[, 1], axis2 = ordination$points[, 2])
  ordination_rows[[view_name]] <- plot_data %>% mutate(view = view_name)
  percent <- 100 * ordination$eig / sum(ordination$eig[ordination$eig > 0])
  ordination_summary_rows[[view_name]] <- tibble(view = view_name,
    n_sites = nrow(pa), lingoes_constant = ordination$ac,
    uncorrected_negative_eigenvalues = sum(uncorrected$eig < -1e-10),
    axis1_percent = percent[1], axis2_percent = percent[2])
  sensitivity_plot <- ggplot(plot_data, aes(axis1, axis2, color = log_area)) +
    geom_point(size = 2.8) + scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
    coord_equal() + labs(x = sprintf("PCoA1 (%.1f%%)", percent[1]),
      y = sprintf("PCoA2 (%.1f%%)", percent[2])) + plot_theme
  ggsave(file.path(figures_dir, paste0("2016_fticr_", view_name, "_pcoa.pdf")),
    sensitivity_plot, width = 6, height = 4.8, bg = "white")
}
write_csv(bind_rows(ordination_rows), file.path(results_dir, "sensitivity_pcoa_scores.csv"))
write_csv(bind_rows(ordination_summary_rows), file.path(results_dir, "sensitivity_pcoa_summary.csv"))
session_lines <- capture.output(sessionInfo())
writeLines(trimws(session_lines, which = "right"), file.path(results_dir, "session_info.txt"))
message("River composition analysis completed: ", nrow(pa), " molecular sites, ",
  nrow(paired_site), " paired sites, ", sum(run_results$completed), " subsampling runs.")
