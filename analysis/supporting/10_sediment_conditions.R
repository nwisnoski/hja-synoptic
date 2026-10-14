# Compare sediment composition and properties using all eligible samples.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "fticr_2016", "sediment_characteristics")
derived_dir <- here("data", "derived", "sediment_characteristics_2016")
figures_dir <- here("figures", "supplement", "sediment_conditions")
n_permutations <- 9999
permutation_seed <- 20261010
enzyme_axes <- 2
target_reads <- 10000
paired_permutation_seed <- 20261080
paired_design <- tribble(
  ~model, ~conditions,
  "reference", "log_area + detected_features",
  "sediment_adjusted", "log_area + detected_features + organic + enzyme_pc1 + enzyme_pc2",
  "hydraulic_subset", "log_area + detected_features + organic + enzyme_pc1 + enzyme_pc2 + hydraulic")
variable_design <- tribble(
  ~variable, ~block, ~transformation, ~label,
  "sediment_organic_content_percent", "organic", "log10", "Organic content",
  "sediment_hydraulic_conductivity_geomean_m_s", "hydraulic", "log10", "Hydraulic conductivity",
  "sediment_nag_umol_g_hr", "enzyme", "asinh_scaled", "NAG activity",
  "sediment_lap_umol_g_hr", "enzyme", "asinh_scaled", "LAP activity",
  "sediment_glu_umol_g_hr", "enzyme", "asinh_scaled", "GLU activity",
  "sediment_ap_umol_g_hr", "enzyme", "asinh_scaled", "AP activity",
  "sediment_d50_um", "texture", "log10", "Median grain size")
comparison_design <- tribble(
  ~comparison, ~tested, ~extra_condition, ~role, ~label,
  "organic", "organic", "", "primary", "Organic content",
  "enzyme", "enzyme_pc1 + enzyme_pc2", "", "primary", "Enzyme profile",
  "hydraulic", "hydraulic", "", "primary", "Hydraulic conductivity",
  "texture", "texture", "", "grain_subset", "Median grain size",
  "organic_unique", "organic", "enzyme_pc1 + enzyme_pc2", "secondary", "Organic content | enzymes",
  "enzyme_unique", "enzyme_pc1 + enzyme_pc2", "organic", "secondary", "Enzymes | organic content")
response_design <- tribble(
  ~response, ~pool, ~label,
  "microbial_hellinger", "microbial", "Microbial composition",
  "molecular_jaccard", "molecular", "Molecular composition",
  "molecular_equal_features", "molecular", "Equal-feature molecular composition")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
response_colors <- c(microbial_hellinger = "#0072B2", molecular_jaccard = "#D55E00",
  molecular_equal_features = "#009E73")

# 2. Audit every sample; missing predictors affect only their own analyses --------
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "river_composition", "mean_subsampled_jaccard.csv"),
  here("results", "fticr_2016", "river_composition", "subsampling_summary.csv"),
  here("data", "derived", "analysis_inputs", "environment", "sediment_44_environment.csv"),
  here("results", "fticr_2016", "integrated_chemistry", "chemistry_site_scores.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
metadata <- read_csv(input_paths[1], col_types = cols(site_code = col_character())) %>%
  filter(habitat == "sediment", included_10k)
molecular <- read_csv(input_paths[3], col_types = cols(site_code = col_character()))
environment <- read_csv(input_paths[6], col_types = cols(site_code = col_character()))
stopifnot(nrow(metadata) == 34, nrow(molecular) == 44, nrow(environment) == 44,
  !anyDuplicated(metadata$site_code), !anyDuplicated(molecular$site_code),
  !anyDuplicated(environment$site_code), setequal(environment$site_code, molecular$site_code))
site_columns <- c("site_code", "site_segment_number", "site_stream_order",
  "site_drainage_area_ha", variable_design$variable)
paired <- metadata %>% filter(site_code %in% molecular$site_code) %>% arrange(site_code)
reference <- environment[match(paired$site_code, environment$site_code), ]
stopifnot(nrow(paired) == 33,
  isTRUE(all.equal(as.data.frame(paired[site_columns]), as.data.frame(reference[site_columns]),
    check.attributes = FALSE)))
site <- bind_rows(metadata %>% select(all_of(site_columns)),
  environment %>% filter(!site_code %in% metadata$site_code) %>% select(all_of(site_columns))) %>%
  arrange(site_code) %>% left_join(metadata %>% select(site_code, sample_id), by = "site_code")
site$microbial_available <- !is.na(site$sample_id)
site$molecular_available <- site$site_code %in% molecular$site_code
pa_all <- as.matrix(molecular[, -1])
rownames(pa_all) <- molecular$site_code
site$detected_features <- rowSums(pa_all)[match(site$site_code, rownames(pa_all))]
stopifnot(nrow(site) == 45, sum(site$microbial_available) == 34,
  sum(site$molecular_available) == 44, all(pa_all %in% c(0, 1)))
write_csv(site, file.path(derived_dir, "all_sediment_site_inputs.csv"))
coverage <- site %>% select(site_code, microbial_available, molecular_available,
  all_of(variable_design$variable)) %>%
  pivot_longer(all_of(variable_design$variable), names_to = "variable", values_to = "value") %>%
  left_join(variable_design, by = "variable") %>% mutate(available = is.finite(value))
write_csv(coverage, file.path(results_dir, "variable_site_coverage.csv"))
write_csv(coverage %>% pivot_longer(c(microbial_available, molecular_available),
  names_to = "pool", values_to = "in_pool") %>% filter(in_pool) %>%
  group_by(pool, variable, label, block) %>% summarize(sites = sum(available),
    missing = sum(!available), zeros = sum(value == 0, na.rm = TRUE), .groups = "drop"),
  file.path(results_dir, "variable_coverage_summary.csv"))

# 3. Build response distances on their full eligible pools ------------------------
counts <- read_csv(input_paths[2], show_col_types = FALSE)
community <- t(as.matrix(counts[, metadata$sample_id]))
colnames(community) <- counts$asv_id
rownames(community) <- metadata$site_code
community <- community[, colSums(community) > 0, drop = FALSE]
stopifnot(all(rowSums(community) == target_reads))
thinned <- read_csv(input_paths[4], col_types = cols(site_code = col_character()))
thinning_summary <- read_csv(input_paths[5], show_col_types = FALSE)
stopifnot(thinning_summary$target_features == 391, thinning_summary$completed_runs == 1000)
equal_distance <- as.matrix(thinned[, thinned$site_code])
rownames(equal_distance) <- thinned$site_code
distances <- list(microbial_hellinger = as.matrix(dist(decostand(community, "hellinger"))),
  molecular_jaccard = as.matrix(vegdist(pa_all, "jaccard", binary = TRUE)),
  molecular_equal_features = equal_distance)

# 4. Transform predictors within each full response pool --------------------------
# asinh(x / median positive x) retains zeros without an arbitrary additive offset.
models <- list()
transformation_rows <- list()
loading_rows <- list()
variance_rows <- list()
model_rows <- list()
for (pool in unique(response_design$pool)) {
  s <- site[site[[paste0(pool, "_available")]], ]
  transformed <- as.matrix(s[variable_design$variable])
  for (j in seq_len(nrow(variable_design))) {
    x <- transformed[, j]
    positive_scale <- NA_real_
    if (variable_design$transformation[j] == "log10") {
      stopifnot(all(x[is.finite(x)] > 0))
      transformed[, j] <- log10(x)
    } else {
      stopifnot(all(x[is.finite(x)] >= 0))
      positive_scale <- median(x[is.finite(x) & x > 0])
      transformed[, j] <- asinh(x / positive_scale)
    }
    transformation_rows[[length(transformation_rows) + 1]] <-
      variable_design[j, ] %>% mutate(pool, asinh_positive_scale = positive_scale)
  }
  enzyme_variables <- variable_design$variable[variable_design$block == "enzyme"]
  enzyme_pca <- prcomp(transformed[, enzyme_variables], center = TRUE, scale. = TRUE)
  for (axis in seq_len(ncol(enzyme_pca$rotation))) {
    largest <- which.max(abs(enzyme_pca$rotation[, axis]))
    if (enzyme_pca$rotation[largest, axis] < 0) {
      enzyme_pca$rotation[, axis] <- -enzyme_pca$rotation[, axis]
      enzyme_pca$x[, axis] <- -enzyme_pca$x[, axis]
    }
  }
  loading_rows[[pool]] <- tibble(pool, variable = enzyme_variables,
    center = enzyme_pca$center, sd = enzyme_pca$scale, as_tibble(enzyme_pca$rotation))
  variance_rows[[pool]] <- tibble(pool, axis = seq_along(enzyme_pca$sdev),
    proportion = enzyme_pca$sdev^2 / sum(enzyme_pca$sdev^2))
  model_data <- s %>% select(site_code, sample_id, site_segment_number,
    site_stream_order, detected_features) %>%
    mutate(log_area = log10(s$site_drainage_area_ha), log_detection = log10(detected_features),
      organic = transformed[, "sediment_organic_content_percent"],
      hydraulic = transformed[, "sediment_hydraulic_conductivity_geomean_m_s"],
      texture = transformed[, "sediment_d50_um"],
      enzyme_pc1 = enzyme_pca$x[, 1], enzyme_pc2 = enzyme_pca$x[, 2])
  stopifnot(enzyme_axes == 2)
  models[[pool]] <- model_data
  model_rows[[pool]] <- model_data %>% mutate(pool)
}
write_csv(bind_rows(transformation_rows), file.path(results_dir, "transformation_design.csv"))
write_csv(bind_rows(loading_rows), file.path(results_dir, "enzyme_pca_loadings.csv"))
write_csv(bind_rows(variance_rows), file.path(results_dir, "enzyme_pca_variance.csv"))
write_csv(bind_rows(model_rows), file.path(derived_dir, "response_pool_model_inputs.csv"))

# 5. Fit small models using only the predictors required for each question --------
# Null: the tested block adds no association after the specified conditions.
# Permute reduced-model residuals. Within-order permutations and segment deletion
# assess sensitivity; neither fully accounts for geographic autocorrelation.
test_rows <- list()
removal_rows <- list()
geometry_rows <- list()
inclusion_rows <- list()
design_rows <- list()
for (r in seq_len(nrow(response_design))) {
  response <- response_design$response[r]
  pool <- response_design$pool[r]
  model_data <- models[[pool]]
  base_condition <- "log_area"
  if (pool == "molecular") base_condition <- c(base_condition, "detected_features")
  for (comparison in seq_len(nrow(comparison_design))) {
    name <- comparison_design$comparison[comparison]
    conditioned <- base_condition
    if (comparison_design$extra_condition[comparison] != "") {
      conditioned <- c(conditioned, comparison_design$extra_condition[comparison])
    }
    formula <- as.formula(paste("y ~", comparison_design$tested[comparison],
      "+ Condition(", paste(conditioned, collapse = " + "), ")"))
    variables <- setdiff(all.vars(formula), "y")
    included <- complete.cases(model_data[variables])
    d <- model_data[included, ]
    inclusion_rows[[length(inclusion_rows) + 1]] <- model_data %>%
      transmute(response, comparison = name, site_code, sample_id, included,
        reason = if_else(included, "Included", "Missing predictor required for this model"))
    raw <- distances[[response]][d$site_code, d$site_code]
    stopifnot(all(is.finite(raw)), max(abs(raw - t(raw))) < 1e-12, all(diag(raw) == 0))
    # All axes preserve the response geometry; correction is refitted per subset.
    ordination <- wcmdscale(as.dist(raw), eig = TRUE, add = "lingoes")
    y <- ordination$points
    corrected <- sqrt(raw^2 + 2 * ordination$ac)
    diag(corrected) <- 0
    error <- max(abs(as.matrix(dist(y)) - corrected))
    stopifnot(error < 1e-7)
    geometry_rows[[length(geometry_rows) + 1]] <- tibble(response, comparison = name,
      sites = nrow(d), axes = ncol(y), lingoes_constant = ordination$ac,
      reconstruction_error = error)
    seed <- permutation_seed + comparison * 10
    set.seed(seed)
    free_permutations <- shuffleSet(nrow(d), nset = n_permutations)
    set.seed(seed + 1)
    order_permutations <- shuffleSet(nrow(d), nset = n_permutations,
      control = how(blocks = factor(d$site_stream_order)))
    fit <- rda(formula, data = d)
    expected_rank <- length(all.vars(as.formula(paste("~", comparison_design$tested[comparison]))))
    stopifnot(fit$CCA$rank == expected_rank, fit$CA$rank > 0)
    free_test <- anova(fit, permutations = free_permutations, model = "reduced")
    order_test <- anova(fit, permutations = order_permutations, model = "reduced")
    fraction <- RsquareAdj(fit)
    role <- comparison_design$role[comparison]
    if (response == "molecular_equal_features") role <- paste0(role, "_sensitivity")
    test_rows[[length(test_rows) + 1]] <- tibble(response, comparison = name, role,
      sites = nrow(d), segments = n_distinct(d$site_segment_number),
      tested_df = fit$CCA$rank, residual_df = free_test$Df[2],
      conditional_raw_fraction = fraction$r.squared,
      conditional_adjusted_fraction = fraction$adj.r.squared,
      pseudo_f = free_test$F[1], p = free_test$`Pr(>F)`[1],
      within_order_p = order_test$`Pr(>F)`[1], status = "ok")
    design_rows[[length(design_rows) + 1]] <- tibble(response, comparison = name,
      formula = paste(deparse(formula), collapse = " "), free_seed = seed,
      within_order_seed = seed + 1, requested_permutations = n_permutations,
      free_permutations = nrow(free_permutations), order_permutations = nrow(order_permutations))
    for (segment in unique(d$site_segment_number)) {
      keep <- d$site_segment_number != segment
      reduced_data <- d[keep, ]
      # Hold the parent geometry and enzyme reference fixed for influence checks.
      y <- ordination$points[keep, , drop = FALSE]
      removal <- tryCatch(rda(formula, data = reduced_data), error = identity)
      status <- "ok"
      fraction_removed <- NA_real_
      if (inherits(removal, "error")) {
        status <- conditionMessage(removal)
      } else if (removal$CCA$rank != expected_rank) {
        status <- "Rank deficient after segment removal"
      } else {
        fraction_removed <- RsquareAdj(removal)$adj.r.squared
      }
      removal_rows[[length(removal_rows) + 1]] <- tibble(response, comparison = name,
        removed_segment = segment, retained_sites = sum(keep),
        conditional_adjusted_fraction = fraction_removed, status)
    }
  }
}
tests <- bind_rows(test_rows) %>% left_join(response_design, by = "response") %>%
  group_by(role) %>% mutate(p_bh = p.adjust(p, "BH"),
    within_order_p_bh = p.adjust(within_order_p, "BH")) %>% ungroup()
write_csv(tests, file.path(results_dir, "composition_characteristic_tests.csv"))
write_csv(bind_rows(inclusion_rows), file.path(results_dir, "model_site_inclusion.csv"))
write_csv(bind_rows(design_rows), file.path(results_dir, "model_permutation_design.csv"))
write_csv(bind_rows(geometry_rows), file.path(results_dir, "response_geometry_checks.csv"))
removals <- bind_rows(removal_rows)
write_csv(removals, file.path(results_dir, "leave_one_segment_out.csv"))
write_csv(removals %>% group_by(response, comparison) %>%
  summarize(min_adjusted = min(conditional_adjusted_fraction, na.rm = TRUE),
    max_adjusted = max(conditional_adjusted_fraction, na.rm = TRUE),
    failed = sum(status != "ok"), .groups = "drop"),
  file.path(results_dir, "leave_one_segment_out_summary.csv"))

# 6. Ask whether the fixed molecular-property association survives adjustment ----
# Direct cross-dataset tests require both profiles; unrelated microbial models
# above retain site 43. Keep the published 42-site molecular PCA reference fixed.
chemistry <- read_csv(input_paths[7], col_types = cols(site_code = col_character())) %>%
  select(site_code, fticr_only_pc1, fticr_only_pc2)
paired_data <- models$microbial %>% filter(site_code %in% molecular$site_code) %>%
  left_join(chemistry, by = "site_code")
stopifnot(nrow(paired_data) == 33, all(complete.cases(paired_data[c("fticr_only_pc1", "fticr_only_pc2")])) )
paired_rows <- list()
paired_inclusion_rows <- list()
for (j in seq_len(nrow(paired_design))) {
  formula <- as.formula(paste("y ~ fticr_only_pc1 + fticr_only_pc2 + Condition(",
    paired_design$conditions[j], ")"))
  variables <- setdiff(all.vars(formula), "y")
  included <- complete.cases(paired_data[variables])
  d <- paired_data[included, ]
  y <- wcmdscale(as.dist(distances$microbial_hellinger[d$site_code, d$site_code]),
    eig = TRUE, add = "lingoes")$points
  fit <- rda(formula, data = d)
  stopifnot(fit$CCA$rank == 2, fit$CA$rank > 0)
  set.seed(paired_permutation_seed)
  free <- shuffleSet(nrow(d), nset = n_permutations)
  set.seed(paired_permutation_seed + 1)
  order <- shuffleSet(nrow(d), nset = n_permutations,
    control = how(blocks = factor(d$site_stream_order)))
  test <- anova(fit, permutations = free, model = "reduced")
  order_test <- anova(fit, permutations = order, model = "reduced")
  paired_rows[[j]] <- tibble(model = paired_design$model[j], sites = nrow(d),
    conditional_adjusted_fraction = RsquareAdj(fit)$adj.r.squared,
    pseudo_f = test$F[1], p = test$`Pr(>F)`[1], within_order_p = order_test$`Pr(>F)`[1],
    formula = paste(deparse(formula), collapse = " "), seed = paired_permutation_seed,
    within_order_seed = paired_permutation_seed + 1)
  paired_inclusion_rows[[j]] <- paired_data %>% transmute(model = paired_design$model[j],
    site_code, sample_id, included)
}
write_csv(bind_rows(paired_rows) %>% mutate(p_bh = p.adjust(p, "BH"),
  within_order_p_bh = p.adjust(within_order_p, "BH")),
  file.path(results_dir, "paired_molecular_property_tests.csv"))
write_csv(bind_rows(paired_inclusion_rows), file.path(results_dir, "paired_model_site_inclusion.csv"))

# 7. Save figures and verify that the inputs remained intact ---------------------
profile_coverage <- site %>% select(site_code, microbial_available, molecular_available) %>%
  pivot_longer(-site_code, names_to = "profile", values_to = "available") %>%
  mutate(label = if_else(profile == "microbial_available", "Microbial library", "Molecular profile"))
coverage_plot <- bind_rows(coverage %>% select(site_code, label, available), profile_coverage) %>%
  mutate(label = factor(label, levels = rev(c("Microbial library", "Molecular profile", variable_design$label)))) %>%
  ggplot(aes(site_code, label, fill = available)) + geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(`FALSE` = "#EEEEEE", `TRUE` = "#0072B2"),
    labels = c(`FALSE` = "Missing", `TRUE` = "Available"), name = NULL) +
  labs(x = "Sediment site (34 microbial and 44 molecular sites; 45 unique)", y = NULL) + plot_theme +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top")
ggsave(file.path(figures_dir, "2016_sediment_characteristics_coverage.pdf"), coverage_plot,
  width = 10, height = 4.2, bg = "white")
effect_plot <- tests %>% filter(grepl("^primary", role)) %>%
  left_join(comparison_design %>% select(comparison, characteristic = label), by = "comparison") %>%
  mutate(characteristic = factor(characteristic,
    levels = c("Organic content", "Hydraulic conductivity", "Enzyme profile")),
    label = factor(label, levels = response_design$label)) %>%
  ggplot(aes(100 * conditional_adjusted_fraction, characteristic, color = response)) +
  geom_vline(xintercept = 0, color = "gray60", linewidth = 0.35) + geom_point(size = 2.7) +
  geom_text(aes(label = paste0("n = ", sites)), nudge_y = 0.20, size = 3, show.legend = FALSE) +
  facet_wrap(~ label, nrow = 1) + scale_color_manual(values = response_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  labs(x = "Adjusted variation beyond conditioning variables (%)", y = NULL) + plot_theme +
  theme(panel.grid.major.y = element_blank())
ggsave(file.path(figures_dir, "2016_sediment_characteristics_unique_fractions.pdf"), effect_plot,
  width = 10, height = 3.4, bg = "white")
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "input_checksums.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
message("Sediment characteristics comparison complete: ", results_dir)
