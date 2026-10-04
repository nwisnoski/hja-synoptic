# Aquatic microbial beta diversity across headwater, intermediate, and mainstem sites.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
results_dir <- here("results", "diversity_2016", "network_dispersion")
figures_dir <- here("figures")
target_depth <- 10000
habitats <- c("sediment", "hyporheic", "planktonic")
habitat_labels <- c(sediment = "Sediment", hyporheic = "Hyporheic", planktonic = "Planktonic")
metrics <- c("bray", "jaccard")
n_permutations <- 9999
n_balanced_draws <- 1000
permutation_seed <- 20261012
resampling_seed <- 20261013
plot_seed <- 20261014
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
group_levels <- c("Headwater", "Intermediate", "Mainstem")
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
contrast_names <- c("Headwater-Intermediate", "Headwater-Mainstem", "Intermediate-Mainstem")
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "none")
test_design <- expand_grid(habitat = habitats, metric = metrics) %>%
  mutate(test_id = row_number(), seed = permutation_seed + test_id,
    role = ifelse(metric == "bray", "primary", "incidence sensitivity"))
draw_design <- expand_grid(habitat = habitats, run = seq_len(n_balanced_draws)) %>%
  mutate(seed = resampling_seed + row_number())

# 2. Align the fixed community table and audit sample/network identities ----------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "StreamCenterline.csv"))
input_checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()))
count_table <- read_csv(input_paths[1], show_col_types = FALSE)
stopifnot(!anyDuplicated(metadata$sample_id), !anyDuplicated(count_table$asv_id))
counts <- t(as.matrix(count_table[, -1]))
colnames(counts) <- count_table$asv_id
storage.mode(counts) <- "numeric"
metadata <- metadata[match(rownames(counts), metadata$sample_id), ]
stopifnot(!anyNA(metadata$sample_id), identical(metadata$sample_id, rownames(counts)),
  all(metadata$included_10k), all(is.finite(counts)), all(counts >= 0),
  all(counts == floor(counts)), all(rowSums(counts) == target_depth))
write_csv(metadata %>% transmute(sample_id, site_code, habitat,
  included = habitat %in% habitats,
  reason = ifelse(habitat %in% habitats, "aquatic fixed 10K sample", "regional soil comparison pool")),
  file.path(results_dir, "sample_inclusion_audit.csv"))
metadata <- metadata %>% filter(habitat %in% habitats) %>%
  select(sample_id, site_code, habitat, site_stream_order, site_segment_number,
    site_utm_x_m, site_utm_y_m, site_drainage_area_ha, has_fticr_profile) %>%
  mutate(network_group = factor(case_when(
    site_stream_order %in% headwater_orders ~ "Headwater",
    site_stream_order %in% intermediate_orders ~ "Intermediate",
    site_stream_order == mainstem_order ~ "Mainstem"), levels = group_levels))
stopifnot(!anyNA(metadata), !anyDuplicated(paste(metadata$habitat, metadata$site_code)))

# Derive basin-outlet distance from the recorded segment, not valley-local D.
geometry <- read_csv(input_paths[3], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
stopifnot(nrow(problems(geometry)) == 0,
  all(is.finite(as.matrix(geometry[, c("x", "y", "outlet_distance", "order", "segment")]))))
metadata$centerline_match_error_m <- NA_real_
metadata$basin_outlet_stream_distance_m <- NA_real_
metadata$centerline_stream_order <- NA_real_
for (i in seq_len(nrow(metadata))) {
  candidates <- which(geometry$segment == metadata$site_segment_number[i])
  stopifnot(length(candidates) > 0)
  distance_squared <- (geometry$x[candidates] - metadata$site_utm_x_m[i])^2 +
    (geometry$y[candidates] - metadata$site_utm_y_m[i])^2
  nearest <- candidates[which.min(distance_squared)]
  metadata$centerline_match_error_m[i] <- sqrt(min(distance_squared))
  metadata$basin_outlet_stream_distance_m[i] <- geometry$outlet_distance[nearest]
  metadata$centerline_stream_order[i] <- geometry$order[nearest]
}
stopifnot(all(metadata$centerline_match_error_m < 1e-5),
  all(metadata$site_stream_order == metadata$centerline_stream_order))
write_csv(metadata, file.path(results_dir, "sample_network_position_audit.csv"))
group_counts <- metadata %>% group_by(habitat, network_group) %>%
  summarize(samples = n(), segments = n_distinct(site_segment_number), .groups = "drop")
stopifnot(nrow(group_counts) == length(habitats) * length(group_levels), all(group_counts$samples >= 3))
write_csv(group_counts, file.path(results_dir, "group_sample_counts.csv"))
write_csv(test_design, file.path(results_dir, "test_design.csv"))
write_csv(draw_design, file.path(results_dir, "balanced_draw_design.csv"))

# 3. Dispersion around spatial medians, with common geometry within each habitat --
# Bray-Curtis describes abundance composition; binary Jaccard assesses incidence.
# Estimate a Lingoes constant on each full habitat matrix and retain it in subsets.
dispersion_rows <- list()
omnibus_rows <- list()
contrast_rows <- list()
correction_rows <- list()
pair_rows <- list()
corrected_distances <- list()
for (test in seq_len(nrow(test_design))) {
  habitat <- test_design$habitat[test]
  metric <- test_design$metric[test]
  sample_data <- metadata %>% filter(.data$habitat == .env$habitat)
  community <- counts[sample_data$sample_id, , drop = FALSE]
  community <- community[, colSums(community) > 0, drop = FALSE]
  raw_distance <- vegdist(community, method = metric, binary = metric == "jaccard")
  raw <- as.matrix(raw_distance)
  ordination <- wcmdscale(raw_distance, eig = TRUE, add = "lingoes")
  correction <- ordination$ac
  if (is.null(correction)) {
    correction <- 0
  }
  corrected <- sqrt(raw^2 + 2 * correction)
  diag(corrected) <- 0
  corrected_distances[[paste(habitat, metric, sep = "_")]] <- corrected
  correction_rows[[test]] <- tibble(habitat, metric, samples = nrow(community),
    asvs = ncol(community), lingoes_constant = correction,
    max_distance_change = max(abs(corrected - raw)))
  fit <- betadisper(as.dist(corrected), sample_data$network_group,
    type = "median", bias.adjust = TRUE)
  set.seed(test_design$seed[test])
  permutations <- shuffleSet(nrow(sample_data), nset = n_permutations)
  tested <- permutest(fit, permutations = permutations, pairwise = TRUE, parallel = 1)
  means <- tapply(fit$distances, sample_data$network_group, mean)
  dispersion_rows[[test]] <- sample_data %>% mutate(metric, distance = fit$distances)
  omnibus_rows[[test]] <- tibble(habitat, metric, role = test_design$role[test],
    seed = test_design$seed[test], samples = nrow(sample_data),
    headwater_mean = means["Headwater"], intermediate_mean = means["Intermediate"],
    mainstem_mean = means["Mainstem"], headwater_mainstem_ratio = means["Headwater"] / means["Mainstem"],
    statistic = tested$tab[1, "F"], p_value = tested$tab[1, "Pr(>F)"])
  for (contrast in contrast_names) {
    groups <- strsplit(contrast, "-", fixed = TRUE)[[1]]
    contrast_rows[[length(contrast_rows) + 1]] <- tibble(habitat, metric, contrast,
      difference = means[groups[1]] - means[groups[2]], ratio = means[groups[1]] / means[groups[2]],
      p_value = tested$pairwise$permuted[contrast])
  }
  indices <- which(upper.tri(raw), arr.ind = TRUE)
  pair_rows[[test]] <- tibble(habitat, metric,
    sample_1 = sample_data$sample_id[indices[, 1]], sample_2 = sample_data$sample_id[indices[, 2]],
    site_1 = sample_data$site_code[indices[, 1]], site_2 = sample_data$site_code[indices[, 2]],
    group_1 = sample_data$network_group[indices[, 1]], group_2 = sample_data$network_group[indices[, 2]],
    same_group = group_1 == group_2,
    same_segment = sample_data$site_segment_number[indices[, 1]] == sample_data$site_segment_number[indices[, 2]],
    geographic_distance_km = sqrt(
      (sample_data$site_utm_x_m[indices[, 1]] - sample_data$site_utm_x_m[indices[, 2]])^2 +
      (sample_data$site_utm_y_m[indices[, 1]] - sample_data$site_utm_y_m[indices[, 2]])^2) / 1000,
    raw_dissimilarity = raw[indices], corrected_dissimilarity = corrected[indices])
}
dispersion <- bind_rows(dispersion_rows)
omnibus <- bind_rows(omnibus_rows) %>% group_by(metric) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
contrasts <- bind_rows(contrast_rows) %>% group_by(metric) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
pairs <- bind_rows(pair_rows)
stopifnot(all(is.finite(dispersion$distance)), !anyNA(omnibus$p_value), !anyNA(contrasts$p_value))
write_csv(dispersion, file.path(results_dir, "site_dispersion.csv"))
write_csv(omnibus, file.path(results_dir, "three_group_tests.csv"))
write_csv(contrasts, file.path(results_dir, "pairwise_group_contrasts.csv"))
write_csv(bind_rows(correction_rows), file.path(results_dir, "distance_correction_audit.csv"))
write_csv(pairs, file.path(results_dir, "microbial_network_pairs.csv"))
write_csv(pairs %>% filter(same_group) %>% group_by(habitat, metric, group_1) %>%
  summarize(pairs = n(), mean_raw_dissimilarity = mean(raw_dissimilarity),
    mean_geographic_distance_km = mean(geographic_distance_km), .groups = "drop"),
  file.path(results_dir, "within_group_pair_summary.csv"))

# 4. Balance sample counts and retain one sample per distinct mapped segment ------
# Keep every mainstem site; its count defines five/five/three samples per group.
# Draws describe design sensitivity, not confidence intervals or independent tests.
balanced_rows <- list()
selection_rows <- list()
start_time <- proc.time()[["elapsed"]]
for (draw in seq_len(nrow(draw_design))) {
  habitat <- draw_design$habitat[draw]
  sample_data <- metadata %>% filter(.data$habitat == .env$habitat)
  mainstem_indices <- which(sample_data$network_group == "Mainstem")
  balanced_n <- length(mainstem_indices)
  stopifnot(n_distinct(sample_data$site_segment_number[mainstem_indices]) == balanced_n)
  selected <- integer(0)
  set.seed(draw_design$seed[draw])
  for (group in group_levels) {
    candidates <- which(sample_data$network_group == group)
    segments <- unique(sample_data$site_segment_number[candidates])
    stopifnot(length(segments) >= balanced_n)
    if (group == "Mainstem") {
      selected <- c(selected, mainstem_indices)
    } else {
      selected_segments <- segments[sample.int(length(segments), balanced_n)]
      for (segment in selected_segments) {
        segment_indices <- candidates[sample_data$site_segment_number[candidates] == segment]
        selected <- c(selected, segment_indices[sample.int(length(segment_indices), 1)])
      }
    }
  }
  labels <- sample_data$network_group[selected]
  stopifnot(all(table(labels) == balanced_n))
  selection_rows[[draw]] <- sample_data[selected, ] %>%
    transmute(habitat, run = draw_design$run[draw], seed = draw_design$seed[draw],
      sample_id, site_code, network_group, segment = site_segment_number)
  for (metric in metrics) {
    corrected <- corrected_distances[[paste(habitat, metric, sep = "_")]]
    diagnostic <- character(0)
    fit <- tryCatch(withCallingHandlers(
      betadisper(as.dist(corrected[selected, selected]), labels, type = "median", bias.adjust = TRUE),
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
    balanced_rows[[length(balanced_rows) + 1]] <- tibble(habitat, metric,
      run = draw_design$run[draw], seed = draw_design$seed[draw], samples_per_group = balanced_n,
      succeeded, diagnostic = paste(diagnostic, collapse = " | "),
      headwater_mean = means["Headwater"], intermediate_mean = means["Intermediate"],
      mainstem_mean = means["Mainstem"],
      headwater_mainstem_ratio = means["Headwater"] / means["Mainstem"],
      intermediate_mainstem_ratio = means["Intermediate"] / means["Mainstem"],
      headwater_intermediate_ratio = means["Headwater"] / means["Intermediate"])
  }
  if (draw == 1) {
    message("Representative two-metric balanced draw: ", round(proc.time()[["elapsed"]] - start_time, 3), " s")
  }
}
balanced <- bind_rows(balanced_rows)
selections <- bind_rows(selection_rows)
stopifnot(nrow(balanced) == nrow(draw_design) * length(metrics))
selection_audit <- selections %>% group_by(habitat, run, network_group) %>%
  summarize(samples = n(), segments = n_distinct(segment), .groups = "drop")
stopifnot(all(selection_audit$samples == selection_audit$segments))
write_csv(selections, file.path(results_dir, "balanced_site_selection.csv"))
write_csv(selection_audit, file.path(results_dir, "balanced_selection_audit.csv"))
write_csv(balanced, file.path(results_dir, "balanced_segment_draws.csv"))
balanced_summary <- balanced %>% group_by(habitat, metric) %>%
  summarize(draws = n(), successful_draws = sum(succeeded), failed_draws = sum(!succeeded),
    draws_with_diagnostics = sum(nchar(diagnostic) > 0),
    median_headwater_mainstem_ratio = median(headwater_mainstem_ratio[succeeded]),
    min_headwater_mainstem_ratio = min(headwater_mainstem_ratio[succeeded]),
    max_headwater_mainstem_ratio = max(headwater_mainstem_ratio[succeeded]),
    fraction_headwater_greater_mainstem = mean(headwater_mainstem_ratio[succeeded] > 1),
    fraction_intermediate_greater_mainstem = mean(intermediate_mainstem_ratio[succeeded] > 1),
    fraction_headwater_greater_intermediate = mean(headwater_intermediate_ratio[succeeded] > 1),
    .groups = "drop")
write_csv(balanced_summary, file.path(results_dir, "balanced_segment_summary.csv"))

# 5. Individual vector panels and a combined primary figure -----------------------
primary_panels <- list()
for (test in seq_len(nrow(test_design))) {
  habitat <- test_design$habitat[test]
  metric <- test_design$metric[test]
  panel <- dispersion %>% filter(.data$habitat == .env$habitat, .data$metric == .env$metric)
  sizes <- table(panel$network_group)
  x_labels <- paste0(group_levels, "\n(n = ", as.integer(sizes[group_levels]), ")")
  set.seed(plot_seed + test)
  p <- ggplot(panel, aes(network_group, distance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.5) +
    geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.85) +
    scale_color_manual(values = group_colors) +
    scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = paste0(habitat_labels[habitat], "\nDistance to group spatial median")) +
    plot_theme
  ggsave(file.path(figures_dir, paste0("2016_", habitat, "_microbes_dispersion_", metric, ".pdf")),
    p, width = 6, height = 4.5, bg = "white")
  if (metric == "bray") {
    primary_panels[[habitat]] <- p
  }
}
combined <- wrap_plots(primary_panels, ncol = 3) + plot_layout(axes = "keep")
ggsave(file.path(figures_dir, "2016_aquatic_microbes_dispersion_bray.pdf"),
  combined, width = 16, height = 4.6, bg = "white")

# 6. Preserve reproducibility and distinguish deferred phylogenetic work ----------
input_checks$final_md5 <- unname(tools::md5sum(input_paths))
input_checks$unchanged <- input_checks$initial_md5 == input_checks$final_md5
stopifnot(all(input_checks$unchanged))
write_csv(input_checks, file.path(results_dir, "input_checksums.csv"))
write_csv(tibble(method = "beta-NTI", calculated = FALSE,
  local_master_tree_present = file.exists(here("results", "phylogeny_2016", "asv_tree_screened_fasttree.nwk")),
  status = "Deferred: verify master tree, phylogenetic signal, null pool and runtime before calculation"),
  file.path(results_dir, "phylogenetic_analysis_status.csv"))
write_csv(tibble(n_permutations, n_balanced_draws, target_depth,
  balanced_runtime_seconds = proc.time()[["elapsed"]] - start_time),
  file.path(results_dir, "run_settings.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(results_dir, "session_info.txt"))
print(omnibus)
print(balanced_summary)
