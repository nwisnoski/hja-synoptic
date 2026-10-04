# Compare molecular and sediment microbial dispersion on identical sites.
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
figures_dir <- here("figures")
target_depth <- 10000
target_features <- 391
n_feature_draws <- 1000
n_permutations <- 9999
n_balanced_draws <- 1000
permutation_seed <- 20261015
resampling_seed <- 20261016
plot_seed <- 20261017
group_levels <- c("Headwater", "Intermediate", "Mainstem")
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
contrast_names <- c("Headwater-Intermediate", "Headwater-Mainstem", "Intermediate-Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
test_design <- tribble(
  ~metric, ~role, ~label,
  "molecular_equal_features", "primary", "Molecular signatures (equal-feature Jaccard)",
  "microbial_bray", "primary", "Sediment microbes (Bray-Curtis)",
  "molecular_jaccard", "sensitivity", "Molecular signatures (complete Jaccard)",
  "microbial_jaccard", "sensitivity", "Sediment microbes (binary Jaccard)") %>%
  mutate(seed = permutation_seed)
draw_design <- tibble(run = seq_len(n_balanced_draws), seed = resampling_seed + run)
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "none")

# 2. Audit the intersection without translating or renaming sample identities -----
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
  list.files(here("results", "diversity_2016", "network_dispersion"), full.names = TRUE),
  here("figures", "2016_fticr_dispersion_equal_features.pdf"),
  here("figures", "2016_sediment_microbes_dispersion_bray.pdf"))
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
group_counts <- site %>% group_by(network_group) %>%
  summarize(sites = n(), segments = n_distinct(site_segment_number), .groups = "drop")
write_csv(group_counts, file.path(results_dir, "group_sample_counts.csv"))

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
  microbial_bray = as.matrix(vegdist(community, method = "bray")),
  molecular_jaccard = as.matrix(vegdist(pa, method = "jaccard", binary = TRUE)),
  microbial_jaccard = as.matrix(vegdist(community, method = "jaccard", binary = TRUE)))
write_csv(tibble(sites = nrow(site), microbial_asvs = ncol(community),
  molecular_features = ncol(pa), target_depth, target_features, n_feature_draws),
  file.path(results_dir, "community_summary.csv"))

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
write_csv(test_design, file.path(results_dir, "test_design.csv"))

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
write_csv(selection_audit, file.path(results_dir, "balanced_selection_audit.csv"))
write_csv(draw_design, file.path(results_dir, "balanced_draw_design.csv"))
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
    geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.85) +
    scale_color_manual(values = group_colors) + scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = paste0(test_design$label[test], "\nDistance to group spatial median")) + plot_theme
  panels[[metric]] <- p
  ggsave(file.path(figures_dir, paste0("2016_matched_sediment_dispersion_", metric, ".pdf")),
    p, width = 6.5, height = 4.5, bg = "white")
}
ggsave(file.path(figures_dir, "2016_matched_sediment_dispersion_comparison.pdf"),
  wrap_plots(panels, ncol = 2), width = 13, height = 4.5, bg = "white")
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
