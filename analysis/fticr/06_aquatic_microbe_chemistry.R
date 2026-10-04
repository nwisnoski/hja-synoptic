# Match water-associated microbes to water optics and contextual sediment chemistry.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "aquatic_microbe_chemistry")
derived_dir <- here("data", "derived", "fticr_2016")
figures_dir <- here("figures")
n_permutations <- 9999
permutation_seed <- 20261011
presence_threshold <- 0
optical_variables <- c("suva254_z", "fluorescence_index_z", "log_t_c_z")
rank_variables <- c("rank_suva254", "rank_fluorescence_index", "rank_log_t_c")
formula_variables <- c("mean_detected_h_c", "mean_detected_o_c", "nitrogen_fraction",
  "protein_like_fraction", "tannin_like_fraction")
habitats <- c("planktonic", "hyporheic")
habitat_labels <- c(planktonic = "Planktonic", hyporheic = "Hyporheic")
habitat_colors <- c(planktonic = "#0072B2", hyporheic = "#D55E00")
group_levels <- c("Headwater", "Intermediate", "Mainstem")
model_design <- tribble(
  ~model_id, ~habitat, ~chemistry,
  1, "planktonic", "matched_water_optics",
  2, "hyporheic", "matched_water_optics",
  3, "planktonic", "sediment_fticr_properties",
  4, "hyporheic", "sediment_fticr_properties")
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "right")

# 2. Retain exact microbial site identities and matched water compartments -------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
inventory <- read_csv(here("data", "derived", "analysis_inputs", "audit", "source_file_inventory.csv"),
  show_col_types = FALSE)
source_checks <- inventory %>% mutate(current_md5 = unname(tools::md5sum(here(source_path))),
  checksum_matches = current_md5 == md5)
stopifnot(all(source_checks$checksum_matches))
write_csv(source_checks, file.path(results_dir, "source_checks.csv"))
metadata <- read_csv(here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  col_types = cols(site_code = col_character())) %>%
  filter(included_10k, habitat %in% habitats) %>%
  select(sample_id, site_code, habitat, plate, site_drainage_area_ha, site_stream_order,
    site_segment_number, site_utm_x_m, site_utm_y_m)
stopifnot(nrow(metadata) == 51, !anyDuplicated(metadata$sample_id),
  !anyDuplicated(paste(metadata$site_code, metadata$habitat)))
metadata$compartment <- ifelse(metadata$habitat == "planktonic", "surface", "hyporheic")
metadata$log_area <- log10(metadata$site_drainage_area_ha)
metadata$network_group <- factor(ifelse(metadata$site_stream_order <= 2, "Headwater",
  ifelse(metadata$site_stream_order == 5, "Mainstem", "Intermediate")), levels = group_levels)
optical <- read_csv(here("data", "derived", "fticr_2016", "dom_optical_site_profiles.csv"),
  col_types = cols(site_code = col_character()))
# Rank all available sites within each water compartment before joining microbes.
# This assesses index/outlier sensitivity without selecting a microbial outcome.
optical <- optical %>% group_by(compartment) %>%
  mutate(rank_suva254 = rank(suva254, na.last = "keep"),
    rank_fluorescence_index = rank(fluorescence_index, na.last = "keep"),
    rank_log_t_c = rank(log_t_c, na.last = "keep")) %>% ungroup()
matched <- metadata %>% left_join(optical %>% select(site_code, compartment,
  complete_primary_optics, optical_pc1, optical_pc2, all_of(optical_variables),
  all_of(rank_variables)), by = c("site_code", "compartment"))
stopifnot(nrow(matched) == 51, all(matched$complete_primary_optics),
  all(is.finite(as.matrix(matched[, c(optical_variables, rank_variables)]))))

# 3. Extend incidence summaries to all 60 profiles; project fixed chemistry PCs --
# The preparation filter is unchanged. Signal defines detection only.
intensity <- readRDS(here("data", "derived", "analysis_inputs", "fticr",
  "fticr_site_by_peak_primary.rds"))
properties <- read_csv(here("results", "fticr_2016", "tables", "primary_peak_properties.csv"),
  show_col_types = FALSE)
properties <- properties[match(colnames(intensity), properties$peak_id), ]
stopifnot(identical(properties$peak_id, colnames(intensity)),
  all(is.finite(intensity)), all(intensity >= 0), !anyDuplicated(rownames(intensity)))
pa <- 1L * (intensity > presence_threshold)
stopifnot(nrow(pa) == 60, ncol(pa) == 4760, all(rowSums(pa) > 0))
write_csv(tibble(site_code = rownames(pa), as_tibble(pa, .name_repair = "minimal")),
  file.path(derived_dir, "site_60_presence_absence.csv"))
formula <- tibble(site_code = rownames(pa), detected_primary_features = rowSums(pa),
  mean_detected_h_c = as.vector(pa %*% properties$hydrogen_carbon_ratio) / rowSums(pa),
  mean_detected_o_c = as.vector(pa %*% properties$oxygen_carbon_ratio) / rowSums(pa),
  nitrogen_fraction = as.vector(pa %*% (properties$nitrogen_count > 0)) / rowSums(pa),
  protein_like_fraction = as.vector(pa %*% (properties$molecular_class_source == "Protein")) / rowSums(pa),
  tannin_like_fraction = as.vector(pa %*% (properties$molecular_class_source == "Tannin")) / rowSums(pa))
loadings <- read_csv(here("results", "fticr_2016", "integrated_chemistry", "chemistry_pca_loadings.csv"),
  show_col_types = FALSE) %>% filter(view == "fticr_only")
loadings <- loadings[match(formula_variables, loadings$variable), ]
scaling <- read_csv(here("results", "fticr_2016", "integrated_chemistry", "integrated_variable_scaling.csv"),
  show_col_types = FALSE)
scaling <- scaling[match(formula_variables, scaling$variable), ]
formula_z <- sweep(as.matrix(formula[, formula_variables]), 2, scaling$center, "-")
formula_z <- sweep(formula_z, 2, scaling$standard_deviation, "/")
projection <- formula_z %*% as.matrix(loadings[, c("pc1_loading", "pc2_loading")])
formula$fticr_pc1 <- projection[, 1]
formula$fticr_pc2 <- projection[, 2]
# Preserve the previously learned 42-site reference, rather than refit to outcomes.
old_scores <- read_csv(here("results", "fticr_2016", "integrated_chemistry", "chemistry_site_scores.csv"),
  col_types = cols(site_code = col_character()))
reproduced <- formula[match(old_scores$site_code, formula$site_code), ]
maximum_projection_error <- max(abs(c(reproduced$fticr_pc1 - old_scores$fticr_only_pc1,
  reproduced$fticr_pc2 - old_scores$fticr_only_pc2)))
old_formula <- read_csv(here("results", "fticr_2016", "river_composition", "site_formula_composition.csv"),
  col_types = cols(site_code = col_character()))
reproduced_traits <- formula[match(old_formula$site_code, formula$site_code), formula_variables]
maximum_trait_error <- max(abs(as.matrix(reproduced_traits) - as.matrix(old_formula[, formula_variables])))
stopifnot(maximum_projection_error < 1e-10, maximum_trait_error < 1e-10)
write_csv(tibble(reference_sites = nrow(old_scores), checked_formula_sites = nrow(old_formula),
  maximum_projection_error, maximum_trait_error), file.path(results_dir, "chemistry_reference_checks.csv"))
write_csv(formula, file.path(derived_dir, "site_60_formula_properties.csv"))
matched <- matched %>% left_join(formula, by = "site_code") %>%
  mutate(has_contextual_fticr = !is.na(detected_primary_features))
write_csv(matched, file.path(results_dir, "aquatic_site_chemistry_matches.csv"))
write_csv(matched %>% group_by(habitat) %>% summarize(microbial_samples = n(),
  matched_water_optics = sum(complete_primary_optics), sediment_fticr_matches = sum(has_contextual_fticr),
  missing_fticr_sites = paste(site_code[!has_contextual_fticr], collapse = ";"), .groups = "drop"),
  file.path(results_dir, "habitat_overlap_audit.csv"))

# 4. Habitat-specific Hellinger composition and descriptive PCA ------------------
counts <- read_csv(here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  show_col_types = FALSE)
stopifnot(!anyDuplicated(counts$asv_id), all(metadata$sample_id %in% names(counts)))
pca_rows <- list()
pca_summary <- list()
for (habitat in habitats) {
  sites <- matched %>% filter(.data$habitat == .env$habitat)
  community <- t(as.matrix(counts[, sites$sample_id]))
  colnames(community) <- counts$asv_id
  rownames(community) <- sites$site_code
  community <- community[, colSums(community) > 0, drop = FALSE]
  stopifnot(all(rowSums(community) == 10000), !anyDuplicated(rownames(community)))
  fit <- prcomp(decostand(community, "hellinger"), center = TRUE, scale. = FALSE)
  percent <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  sites$microbial_pc1 <- fit$x[, 1]
  sites$microbial_pc2 <- fit$x[, 2]
  pca_rows[[habitat]] <- sites
  pca_summary[[habitat]] <- tibble(habitat, sites = nrow(sites), asvs = ncol(community),
    pc1_percent = percent[1], pc2_percent = percent[2])
  p <- ggplot(sites, aes(microbial_pc1, microbial_pc2, color = optical_pc1, shape = network_group)) +
    geom_point(size = 2.8) + scale_color_viridis_c(option = "C",
      limits = range(matched$optical_pc1), name = "Matched water\noptical PC1") +
    scale_shape_manual(values = c(16, 17, 15), drop = FALSE, name = "Network group") +
    guides(color = guide_colorbar(order = 1), shape = guide_legend(order = 2)) +
    labs(x = paste0("Microbial PC1 (", round(percent[1], 1), "%)"),
      y = paste0("Microbial PC2 (", round(percent[2], 1), "%)")) + plot_theme
  ggsave(file.path(figures_dir, paste0("2016_", habitat, "_microbes_water_optics.pdf")),
    p, width = 7.3, height = 4.8, bg = "white")
}
write_csv(bind_rows(pca_rows), file.path(results_dir, "aquatic_microbial_pca_scores.csv"))
write_csv(bind_rows(pca_summary), file.path(results_dir, "aquatic_microbial_pca_summary.csv"))

# 5. Four predeclared composition models and full-distance concordance -----------
# Analyze habitats separately so repeated sites across habitats are not replicates.
# Water optics use all three variables, rather than discard optical PC3.
model_rows <- list()
mantel_rows <- list()
rank_rows <- list()
validation_rows <- list()
influence_rows <- list()
for (i in seq_len(nrow(model_design))) {
  design <- model_design[i, ]
  sites <- matched %>% filter(habitat == design$habitat)
  if (design$chemistry == "sediment_fticr_properties") {
    sites <- filter(sites, has_contextual_fticr)
  }
  community <- t(as.matrix(counts[, sites$sample_id]))
  rownames(community) <- sites$site_code
  community <- community[, colSums(community) > 0, drop = FALSE]
  hellinger <- decostand(community, "hellinger")
  coordinates <- prcomp(hellinger, center = TRUE, scale. = FALSE)$x
  geometry_error <- max(abs(as.vector(dist(coordinates)) - as.vector(dist(hellinger))))
  stopifnot(geometry_error < 1e-10, all(rowSums(community) == 10000))
  sites$chemistry_pc1 <- sites$fticr_pc1
  sites$chemistry_pc2 <- sites$fticr_pc2
  model_formula <- hellinger ~ suva254_z + fluorescence_index_z + log_t_c_z + Condition(log_area)
  if (design$chemistry == "sediment_fticr_properties") {
    model_formula <- hellinger ~ chemistry_pc1 + chemistry_pc2 +
      Condition(log_area + detected_primary_features)
  }
  full_fit <- rda(model_formula, data = sites)
  model_formula[[2]] <- quote(coordinates)
  fit <- rda(model_formula, data = sites)
  inertia_error <- max(abs(c(fit$tot.chi - full_fit$tot.chi,
    fit$CCA$tot.chi - full_fit$CCA$tot.chi, fit$pCCA$tot.chi - full_fit$pCCA$tot.chi)))
  adjusted_error <- abs(RsquareAdj(fit)$adj.r.squared - RsquareAdj(full_fit)$adj.r.squared)
  stopifnot(inertia_error < 1e-10, adjusted_error < 1e-10)
  validation_rows[[i]] <- tibble(model_id = design$model_id, geometry_error, inertia_error,
    adjusted_error, response_asvs = ncol(hellinger), response_coordinates = ncol(coordinates))
  set.seed(permutation_seed + design$model_id)
  permutations <- shuffleSet(nrow(sites), nset = n_permutations)
  set.seed(permutation_seed + 100 + design$model_id)
  order_permutations <- shuffleSet(nrow(sites), nset = n_permutations,
    control = how(blocks = factor(sites$site_stream_order)))
  permutation_result <- anova(fit, permutations = permutations)
  order_result <- anova(fit, permutations = order_permutations)
  model_rows[[i]] <- tibble(habitat = design$habitat, chemistry = design$chemistry,
    model_id = design$model_id, sites = nrow(sites), seed = permutation_seed + design$model_id,
    constraints = fit$CCA$rank, condition_rank = fit$pCCA$rank,
    residual_df = df.residual(fit), fraction_total_variation = fit$CCA$tot.chi / fit$tot.chi,
    adjusted_unique_fraction = RsquareAdj(fit)$adj.r.squared,
    pseudo_f = permutation_result$F[1], p_value = permutation_result$`Pr(>F)`[1],
    p_within_stream_order = order_result$`Pr(>F)`[1],
    used_order_permutations = nrow(order_permutations))
  chemistry_distance <- dist(as.matrix(sites[, optical_variables]))
  if (design$chemistry == "sediment_fticr_properties") {
    chemistry_distance <- vegdist(pa[sites$site_code, ], method = "jaccard", binary = TRUE)
  }
  concordance <- mantel(dist(hellinger), chemistry_distance, method = "spearman",
    permutations = permutations, parallel = 1)
  order_concordance <- mantel(dist(hellinger), chemistry_distance, method = "spearman",
    permutations = order_permutations, parallel = 1)
  mantel_rows[[i]] <- tibble(habitat = design$habitat,
    signature = ifelse(design$chemistry == "matched_water_optics", "Matched water optical distance",
      "Contextual sediment binary Jaccard"), sites = nrow(sites),
    spearman_rho = concordance$statistic, p_value = concordance$signif,
    p_within_stream_order = order_concordance$signif)
  if (design$chemistry == "matched_water_optics") {
    rank_fit <- rda(coordinates ~ rank_suva254 + rank_fluorescence_index + rank_log_t_c +
      Condition(log_area), data = sites)
    rank_result <- anova(rank_fit, permutations = permutations)
    rank_rows[[design$habitat]] <- tibble(habitat = design$habitat, sites = nrow(sites),
      adjusted_unique_fraction = RsquareAdj(rank_fit)$adj.r.squared,
      p_value = rank_result$`Pr(>F)`[1])
  }
  # Leave-one-site-out fractions assess influence, not predictive performance.
  for (j in seq_len(nrow(sites))) {
    reduced_coordinates <- coordinates[-j, , drop = FALSE]
    leave_out_formula <- model_formula
    leave_out_formula[[2]] <- quote(reduced_coordinates)
    leave_out_fit <- rda(leave_out_formula, data = sites[-j, ])
    influence_rows[[length(influence_rows) + 1]] <- tibble(model_id = design$model_id,
      omitted_site = sites$site_code[j], adjusted_unique_fraction = RsquareAdj(leave_out_fit)$adj.r.squared)
  }
}
models <- bind_rows(model_rows) %>% mutate(p_bh_four_models = p.adjust(p_value, "BH"),
  p_order_bh_four_models = p.adjust(p_within_stream_order, "BH"))
write_csv(models, file.path(results_dir, "aquatic_chemistry_partial_rda.csv"))
write_csv(bind_rows(mantel_rows) %>% mutate(p_bh_four_tests = p.adjust(p_value, "BH"),
  p_order_bh_four_tests = p.adjust(p_within_stream_order, "BH")),
  file.path(results_dir, "aquatic_chemistry_mantel.csv"))
write_csv(bind_rows(rank_rows) %>% mutate(p_bh_two_tests = p.adjust(p_value, "BH")),
  file.path(results_dir, "water_optics_rank_sensitivity.csv"))
write_csv(bind_rows(validation_rows), file.path(results_dir, "microbial_coordinate_equivalence_checks.csv"))
write_csv(bind_rows(influence_rows), file.path(results_dir, "leave_one_site_out_fractions.csv"))
write_csv(bind_rows(influence_rows) %>% group_by(model_id) %>%
  summarize(minimum = min(adjusted_unique_fraction), maximum = max(adjusted_unique_fraction),
    median = median(adjusted_unique_fraction), .groups = "drop") %>% left_join(model_design, by = "model_id"),
  file.path(results_dir, "leave_one_site_out_summary.csv"))

# 6. Summarize associations without truncating negative adjusted estimates -------
p <- ggplot(models, aes(factor(chemistry, levels = c("matched_water_optics", "sediment_fticr_properties")),
  100 * adjusted_unique_fraction, color = habitat)) +
  geom_hline(yintercept = 0, color = "#777777", linetype = "dashed") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = habitat_colors, labels = habitat_labels, name = NULL) +
  scale_x_discrete(labels = c("Matched water optics", "Sediment FT-ICR properties")) +
  coord_flip() + labs(x = NULL, y = "Adjusted unique microbial variation (%)") +
  plot_theme + theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "2016_aquatic_microbes_chemistry_associations.pdf"),
  p, width = 7.5, height = 4.2, bg = "white")
capture.output(sessionInfo(), file = file.path(results_dir, "session_info.txt"))
print(models)
print(read_csv(file.path(results_dir, "aquatic_chemistry_mantel.csv"), show_col_types = FALSE))
