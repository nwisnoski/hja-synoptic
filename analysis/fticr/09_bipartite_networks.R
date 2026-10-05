# Infer exploratory sediment microbe-molecule groups and compare explicit nulls.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(igraph)

# 1. Settings and plotting theme ------------------------------------------------
test_run <- FALSE
target_depth <- 10000
minimum_present_and_absent <- 10L
edge_threshold <- 0.6
threshold_sensitivity <- c(0.5, 0.7)
n_null <- 999L
null_seed <- 20261040L
observed_seed <- 20261041L
resampling_seed <- 20261042L
start_groups <- c(0L, 0L, 4L, 8L, 16L, 32L, 64L, 128L)
max_sweeps <- 100L
modularity_tolerance <- 1e-10
quasiswap_thin <- 100L
minimum_group_nodes <- 3L
results_dir <- here("results", "fticr_2016", "bipartite_networks")
figure_prefix <- "2016_sediment_bipartite"
figures_dir <- here("figures")
if (test_run) {
  n_null <- 9L
  results_dir <- here("results", "fticr_2016", "bipartite_network_pilot")
  figure_prefix <- "pilot_2016_sediment_bipartite"
}
null_design <- expand_grid(
  scheme = c("profile_unrestricted", "profile_within_order", "molecular_fixed_margins", "graph_fixed_degrees"),
  run = seq_len(n_null)) %>%
  mutate(seed = null_seed + match(scheme, unique(scheme)) * 10000L + run,
    module_seed = seed + 100000L)
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
source(here("analysis", "fticr", "bipartite_helpers.R"))

# 2. Read and audit the paired inputs -------------------------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_metadata_screened.csv"),
  here("results", "fticr_2016", "tables", "primary_peak_properties.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
stopifnot(!anyNA(checks$initial_md5))
metadata <- read_csv(input_paths[1], col_types = cols(site_code = col_character()),
  show_col_types = FALSE) %>% filter(included_10k, habitat == "sediment")
molecular <- read_csv(input_paths[3], col_types = cols(site_code = col_character()),
  show_col_types = FALSE)
stopifnot(!anyDuplicated(metadata$site_code), !anyDuplicated(metadata$sample_id),
  !anyDuplicated(molecular$site_code))
inclusion <- full_join(metadata %>% transmute(site_code, sample_id, sediment_10k = TRUE),
  molecular %>% transmute(site_code, molecular_profile = TRUE), by = "site_code") %>%
  mutate(sediment_10k = replace_na(sediment_10k, FALSE),
    molecular_profile = replace_na(molecular_profile, FALSE),
    included = sediment_10k & molecular_profile)
write_csv(inclusion, file.path(results_dir, "site_inclusion_audit.csv"))
site <- metadata %>% filter(site_code %in% molecular$site_code) %>% arrange(site_code) %>%
  select(site_code, sample_id, site_segment_number, site_stream_order, site_drainage_area_ha,
    site_utm_x_m, site_utm_y_m)
stopifnot(!anyNA(site), all(site$site_drainage_area_ha > 0))
count_table <- read_csv(input_paths[2], show_col_types = FALSE)
stopifnot(!anyDuplicated(count_table$asv_id), all(site$sample_id %in% names(count_table)))
counts <- t(as.matrix(count_table[, site$sample_id]))
colnames(counts) <- count_table$asv_id
rownames(counts) <- site$site_code
molecular_pa <- as.matrix(molecular[match(site$site_code, molecular$site_code), -1])
rownames(molecular_pa) <- site$site_code
stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)),
  all(rowSums(counts) == target_depth), all(molecular_pa %in% c(0, 1)))
site$detected_molecular_features <- rowSums(molecular_pa)
site$detected_asvs <- rowSums(counts > 0)
write_csv(site, file.path(results_dir, "matched_site_manifest.csv"))
taxonomy <- read_csv(input_paths[4], show_col_types = FALSE) %>% select(-sequence)
properties <- read_csv(input_paths[5], show_col_types = FALSE)
stopifnot(!anyDuplicated(taxonomy$asv_id), !anyDuplicated(properties$peak_id))

# 3. Filter prevalence and collapse observationally identical incidence patterns -
asv_prevalence <- colSums(counts > 0)
molecular_prevalence <- colSums(molecular_pa)
feature_audit <- bind_rows(
  tibble(node_type = "ASV", feature_id = colnames(counts), present = asv_prevalence),
  tibble(node_type = "Molecular", feature_id = colnames(molecular_pa), present = molecular_prevalence)) %>%
  mutate(absent = nrow(site) - present,
    included = present >= minimum_present_and_absent & absent >= minimum_present_and_absent)
write_csv(feature_audit, file.path(results_dir, "feature_inclusion_audit.csv"))
pattern_maps <- list()
matrices <- list()
for (node_type in c("ASV", "Molecular")) {
  if (node_type == "ASV") {
    incidence <- 1L * (counts > 0)
    prefix <- "ASVP"
  } else {
    incidence <- molecular_pa
    prefix <- "MOLP"
  }
  selected <- feature_audit$feature_id[feature_audit$node_type == node_type & feature_audit$included]
  incidence <- incidence[, selected, drop = FALSE]
  signatures <- apply(incidence, 2, paste0, collapse = "")
  representatives <- which(!duplicated(signatures))
  pattern <- match(signatures, signatures[representatives])
  ids <- sprintf("%s_%04d", prefix, seq_along(representatives))
  pattern_maps[[node_type]] <- tibble(node_type, feature_id = selected,
    pattern_id = ids[pattern], representative_id = selected[representatives][pattern])
  matrices[[node_type]] <- incidence[, representatives, drop = FALSE]
  colnames(matrices[[node_type]]) <- ids
}
pattern_map <- bind_rows(pattern_maps)
write_csv(pattern_map, file.path(results_dir, "feature_pattern_map.csv"))
x <- matrices$ASV
y <- matrices$Molecular
write_csv(tibble(site_code = site$site_code, as_tibble(x)), file.path(results_dir, "asv_pattern_incidence.csv"))
write_csv(tibble(site_code = site$site_code, as_tibble(y)), file.path(results_dir, "molecular_pattern_incidence.csv"))
stopifnot(all(apply(x, 2, sd) > 0), all(apply(y, 2, sd) > 0))
associations <- cor(x, y)
set.seed(observed_seed)
observed <- bipartite_statistics(associations, edge_threshold, start_groups,
  max_sweeps, modularity_tolerance)
write_csv(as_tibble(observed$statistics), file.path(results_dir, "observed_network.csv"))
write_csv(observed$fit$starts, file.path(results_dir, "observed_optimizer_starts.csv"))
web <- 1L * (associations >= edge_threshold)
active_rows <- which(rowSums(web) > 0)
active_columns <- which(colSums(web) > 0)
active_web <- web[active_rows, active_columns, drop = FALSE]
selected_edges <- which(abs(associations) >= edge_threshold, arr.ind = TRUE)
edges <- tibble(asv_pattern = rownames(associations)[selected_edges[, 1]],
  molecular_pattern = colnames(associations)[selected_edges[, 2]],
  asv_index = selected_edges[, 1], molecular_index = selected_edges[, 2],
  phi = associations[selected_edges], sign = if_else(phi > 0, "Positive", "Negative"))
memberships <- bind_rows(
  tibble(node_type = "ASV", pattern_id = colnames(x), module = observed$fit$row, degree = rowSums(web)),
  tibble(node_type = "Molecular", pattern_id = colnames(y), module = observed$fit$column, degree = colSums(web)))
write_csv(memberships, file.path(results_dir, "pattern_membership.csv"))

# 4. Rebuild networks under data nulls and fixed-degree bipartite topology nulls --
# Quasiswap is nonsequential: each draw preserves both matrix margins.
molecular_null <- nullmodel(y, method = "quasiswap")
graph_null <- nullmodel(active_web, method = "quasiswap")
null_rows <- vector("list", nrow(null_design))
benchmark_start <- proc.time()[["elapsed"]]
for (run in seq_len(nrow(null_design))) {
  design <- null_design[run, ]
  set.seed(design$seed)
  start_time <- proc.time()[["elapsed"]]
  diagnostic <- NA_character_
  margins_preserved <- NA
  changed_fraction <- NA_real_
  fit <- tryCatch({
    if (design$scheme == "graph_fixed_degrees") {
      random_web <- simulate(graph_null, nsim = 1, thin = quasiswap_thin)[, , 1]
      margins_preserved <- all(rowSums(random_web) == rowSums(active_web)) &
        all(colSums(random_web) == colSums(active_web)) & all(random_web %in% c(0, 1))
      stopifnot(margins_preserved)
      changed_fraction <- mean(random_web != active_web)
      set.seed(design$module_seed)
      module_fit <- bipartite_modules(random_web, start_groups, max_sweeps, modularity_tolerance)
      list(statistics = data.frame(positive_edges = sum(random_web), negative_edges = NA_integer_,
        active_asv_patterns = nrow(random_web), active_molecular_patterns = ncol(random_web),
        modularity = module_fit$modularity, modules = module_fit$modules, converged = module_fit$converged))
    } else {
      if (design$scheme == "molecular_fixed_margins") {
        random_y <- simulate(molecular_null, nsim = 1, thin = quasiswap_thin)[, , 1]
        margins_preserved <- all(rowSums(random_y) == rowSums(y)) &
          all(colSums(random_y) == colSums(y)) & all(random_y %in% c(0, 1))
        stopifnot(margins_preserved)
      } else {
        permutation <- seq_len(nrow(site))
        if (design$scheme == "profile_unrestricted") {
          permutation <- sample.int(nrow(site))
        } else {
          for (order in sort(unique(site$site_stream_order))) {
            positions <- which(site$site_stream_order == order)
            permutation[positions] <- positions[sample.int(length(positions))]
          }
        }
        random_y <- y[permutation, , drop = FALSE]
      }
      changed_fraction <- mean(random_y != y)
      random_associations <- cor(x, random_y)
      set.seed(design$module_seed)
      bipartite_statistics(random_associations, edge_threshold, start_groups,
        max_sweeps, modularity_tolerance)
    }
  }, error = function(e) {
    diagnostic <<- conditionMessage(e)
    NULL
  })
  if (is.null(fit)) {
    statistics <- tibble(positive_edges = NA_integer_, negative_edges = NA_integer_,
      active_asv_patterns = NA_integer_, active_molecular_patterns = NA_integer_,
      modularity = NA_real_, modules = NA_integer_, converged = FALSE)
  } else {
    statistics <- as_tibble(fit$statistics)
  }
  null_rows[[run]] <- bind_cols(design, statistics) %>%
    mutate(margins_preserved, changed_fraction, diagnostic,
      successful = is.na(diagnostic) & converged,
      elapsed_seconds = proc.time()[["elapsed"]] - start_time)
  if (run %% 50L == 0L || run == nrow(null_design)) {
    write_csv(bind_rows(null_rows), file.path(results_dir, "null_runs.csv"))
    message("Null run ", run, "/", nrow(null_design), "; elapsed ",
      round(proc.time()[["elapsed"]] - benchmark_start, 1), " s")
  }
}
nulls <- bind_rows(null_rows)
write_csv(null_design, file.path(results_dir, "null_design.csv"))
test_rows <- list()
for (scheme in unique(nulls$scheme)) {
  metrics <- c("positive_edges", "negative_edges", "modularity")
  if (scheme == "graph_fixed_degrees") {
    metrics <- "modularity"
  }
  for (metric in metrics) {
    eligible <- nulls$scheme == scheme & nulls$successful & is.finite(nulls[[metric]])
    values <- nulls[[metric]][eligible]
    observed_value <- observed$statistics[[metric]]
    complete <- length(values) == n_null
    test_rows[[length(test_rows) + 1L]] <- tibble(scheme, metric,
      observed = observed_value, completed = length(values), expected = n_null,
      null_mean = mean(values), null_sd = sd(values),
      null_q025 = unname(quantile(values, 0.025)), null_q975 = unname(quantile(values, 0.975)),
      p_value = if (complete) (1 + sum(values >= observed_value - modularity_tolerance)) / (n_null + 1) else NA_real_,
      complete)
  }
}
tests <- bind_rows(test_rows) %>%
  mutate(family = if_else(metric == "modularity", "Module structure", "Association counts")) %>%
  group_by(family) %>% mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
write_csv(tests, file.path(results_dir, "null_tests.csv"))

# 5. Refit groups after removing entire stream segments --------------------------
segments <- sort(unique(site$site_segment_number))
resampling_design <- tibble(run = seq_along(segments), omitted_segment = segments,
  seed = resampling_seed + seq_along(segments))
positive_edges <- edges %>% filter(sign == "Positive")
within_module_pairs <- which(outer(observed$fit$row, observed$fit$column, "==") &
  outer(!is.na(observed$fit$row), !is.na(observed$fit$column)), arr.ind = TRUE)
pair_stability <- matrix(FALSE, nrow(within_module_pairs), length(segments))
edge_stability <- matrix(FALSE, nrow(edges), length(segments))
edge_coassignment <- matrix(FALSE, nrow(positive_edges), length(segments))
resampling_rows <- list()
membership_rows <- list()
for (run in seq_along(segments)) {
  retained <- site$site_segment_number != segments[run]
  refit_associations <- cor(x[retained, , drop = FALSE], y[retained, , drop = FALSE])
  set.seed(resampling_design$seed[run])
  refit <- bipartite_statistics(refit_associations, edge_threshold, start_groups,
    max_sweeps, modularity_tolerance)
  edge_phi <- refit_associations[cbind(edges$asv_index, edges$molecular_index)]
  edge_stability[, run] <- if_else(edges$sign == "Positive",
    edge_phi >= edge_threshold, edge_phi <= -edge_threshold)
  coassignment <- refit$fit$row[positive_edges$asv_index] ==
    refit$fit$column[positive_edges$molecular_index]
  coassignment[is.na(coassignment)] <- FALSE
  edge_coassignment[, run] <- coassignment
  same_pair <- refit$fit$row[within_module_pairs[, 1]] == refit$fit$column[within_module_pairs[, 2]]
  same_pair[is.na(same_pair)] <- FALSE
  pair_stability[, run] <- same_pair
  resampling_rows[[run]] <- bind_cols(resampling_design[run, ], as_tibble(refit$statistics)) %>%
    mutate(retained_sites = sum(retained))
  membership_rows[[run]] <- bind_rows(
    tibble(node_type = "ASV", pattern_id = colnames(x), module = refit$fit$row),
    tibble(node_type = "Molecular", pattern_id = colnames(y), module = refit$fit$column)) %>%
    mutate(omitted_segment = segments[run])
}
edges$segment_threshold_support <- rowMeans(edge_stability)
edges$segment_coassignment <- NA_real_
edges$segment_coassignment[edges$sign == "Positive"] <- rowMeans(edge_coassignment)
write_csv(bind_rows(resampling_rows), file.path(results_dir, "leave_segment_out_runs.csv"))
write_csv(bind_rows(membership_rows), file.path(results_dir, "leave_segment_out_memberships.csv"))
pairs <- tibble(asv_pattern = colnames(x)[within_module_pairs[, 1]],
  molecular_pattern = colnames(y)[within_module_pairs[, 2]],
  module = observed$fit$row[within_module_pairs[, 1]],
  segment_coassignment = rowMeans(pair_stability))
write_csv(pairs, file.path(results_dir, "within_module_pair_stability.csv"))

# 6. Examine detection/gradient adjustment and two edge-threshold sensitivities --
# Partial correlations are descriptive; there is no conditional permutation P.
adjustment <- model.matrix(~ log10(site_drainage_area_ha) + detected_molecular_features,
  data = site)
stopifnot(qr(adjustment)$rank == ncol(adjustment))
x_residual <- qr.resid(qr(adjustment), x)
y_residual <- qr.resid(qr(adjustment), y)
adjusted_associations <- cor(x_residual, y_residual)
edges$partial_phi <- adjusted_associations[cbind(edges$asv_index, edges$molecular_index)]
write_csv(edges %>% select(-asv_index, -molecular_index), file.path(results_dir, "candidate_edges.csv"))
sensitivity_rows <- list()
for (threshold in c(edge_threshold, threshold_sensitivity)) {
  set.seed(observed_seed)
  fitted <- bipartite_statistics(associations, threshold, start_groups, max_sweeps, modularity_tolerance)
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- as_tibble(fitted$statistics) %>%
    mutate(view = "Raw incidence", threshold)
}
set.seed(observed_seed)
adjusted <- bipartite_statistics(adjusted_associations, edge_threshold, start_groups,
  max_sweeps, modularity_tolerance)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- as_tibble(adjusted$statistics) %>%
  mutate(view = "Detection and drainage adjusted", threshold = edge_threshold)
write_csv(bind_rows(sensitivity_rows), file.path(results_dir, "network_sensitivities.csv"))
write_csv(bind_rows(
  tibble(node_type = "ASV", pattern_id = colnames(x), module = adjusted$fit$row),
  tibble(node_type = "Molecular", pattern_id = colnames(y), module = adjusted$fit$column)),
  file.path(results_dir, "adjusted_pattern_membership.csv"))

# 7. Map modules to the original features and summarize their composition --------
feature_membership <- pattern_map %>% left_join(memberships, by = c("node_type", "pattern_id"))
write_csv(feature_membership, file.path(results_dir, "feature_membership.csv"))
asv_members <- feature_membership %>% filter(node_type == "ASV") %>%
  left_join(taxonomy, by = c("feature_id" = "asv_id"))
molecular_members <- feature_membership %>% filter(node_type == "Molecular") %>%
  left_join(properties, by = c("feature_id" = "peak_id"))
stopifnot(!anyNA(asv_members$Kingdom), !anyNA(molecular_members$carbon_count))
write_csv(asv_members, file.path(results_dir, "asv_module_taxonomy.csv"))
write_csv(molecular_members, file.path(results_dir, "molecular_module_properties.csv"))
module_summary <- memberships %>% filter(!is.na(module)) %>% group_by(module) %>%
  summarize(asv_patterns = sum(node_type == "ASV"), molecular_patterns = sum(node_type == "Molecular"),
    edges = sum(degree[node_type == "ASV"]), .groups = "drop")
positive_edges <- edges %>% filter(sign == "Positive") %>%
  mutate(asv_module = observed$fit$row[asv_index], molecular_module = observed$fit$column[molecular_index])
module_summary <- module_summary %>%
  left_join(positive_edges %>% filter(asv_module == molecular_module) %>% group_by(module = asv_module) %>%
    summarize(within_module_edges = n(), median_phi = median(phi),
      mean_edge_support = mean(segment_threshold_support),
      mean_partial_phi = mean(partial_phi), .groups = "drop"), by = "module") %>%
  left_join(pairs %>% group_by(module) %>% summarize(mean_pair_coassignment = mean(segment_coassignment),
    fraction_pairs_support_80 = mean(segment_coassignment >= 0.8), .groups = "drop"), by = "module") %>%
  left_join(asv_members %>% filter(!is.na(module)) %>% count(module, name = "original_asvs"), by = "module") %>%
  left_join(molecular_members %>% filter(!is.na(module)) %>% group_by(module) %>%
    summarize(original_molecular_features = n(), mean_h_c = mean(hydrogen_carbon_ratio),
      mean_o_c = mean(oxygen_carbon_ratio), fraction_n = mean(nitrogen_count > 0),
      fraction_s = mean(sulfur_count > 0), fraction_p = mean(phosphorus_count > 0), .groups = "drop"), by = "module") %>%
  mutate(substantive_group = asv_patterns >= minimum_group_nodes & molecular_patterns >= minimum_group_nodes) %>%
  arrange(desc(within_module_edges))
write_csv(module_summary, file.path(results_dir, "module_summary.csv"))
write_csv(asv_members %>% filter(!is.na(module)) %>% count(module, Phylum, name = "asvs"),
  file.path(results_dir, "module_phylum_counts.csv"))
profile_rows <- list()
for (module_id in module_summary$module) {
  module_asvs <- which(observed$fit$row == module_id)
  module_molecules <- which(observed$fit$column == module_id)
  asv_fraction <- rowMeans(x[, module_asvs, drop = FALSE])
  molecular_fraction <- rowMeans(y[, module_molecules, drop = FALSE])
  profile_rows[[length(profile_rows) + 1L]] <- site %>% transmute(site_code, site_stream_order,
    site_segment_number, site_drainage_area_ha, detected_molecular_features, module = .env$module_id,
    asv_detection_fraction = .env$asv_fraction,
    molecular_detection_fraction = .env$molecular_fraction)
}
write_csv(bind_rows(profile_rows), file.path(results_dir, "module_site_profiles.csv"))

# 8. Save compact null and group figures -----------------------------------------
labels <- c(profile_unrestricted = "Whole profiles", profile_within_order = "Profiles within stream order",
  molecular_fixed_margins = "Molecular incidence with fixed margins", graph_fixed_degrees = "Bipartite graph with fixed degrees")
edge_plot <- ggplot(nulls %>% filter(scheme != "graph_fixed_degrees", successful), aes(positive_edges)) +
  geom_histogram(bins = 30, fill = "#0072B2", color = "white") +
  geom_vline(xintercept = observed$statistics$positive_edges, color = "#D55E00", linewidth = 0.7) +
  facet_wrap(~scheme, scales = "free", ncol = 1, labeller = as_labeller(labels)) +
  labs(x = "Positive cross-type links (phi >= 0.6)", y = "Null draws") + plot_theme +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))
ggsave(file.path(figures_dir, paste0(figure_prefix, "_association_nulls.pdf")), edge_plot,
  width = 6.2, height = 7.2, bg = "white")
module_plot <- ggplot(nulls %>% filter(scheme == "graph_fixed_degrees", successful), aes(modularity)) +
  geom_histogram(bins = 30, fill = "#0072B2", color = "white") +
  geom_vline(xintercept = observed$statistics$modularity, color = "#D55E00", linewidth = 0.7) +
  labs(x = "Optimized bipartite modularity", y = "Fixed-degree null draws") + plot_theme
ggsave(file.path(figures_dir, paste0(figure_prefix, "_module_null.pdf")), module_plot,
  width = 6, height = 4.6, bg = "white")
group_plot_data <- module_summary %>% filter(substantive_group) %>%
  mutate(label = paste0(module, " (", asv_patterns, " / ", molecular_patterns, ")"),
    label = factor(label, levels = rev(label)))
summary_plot <- ggplot(group_plot_data,
  aes(mean_pair_coassignment, label, size = within_module_edges, color = mean_partial_phi)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", midpoint = 0,
    limits = c(-1, 1), name = "Mean partial phi") +
  scale_size_continuous(name = "Within-module links", range = c(3, 9)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Mean co-membership after stream-segment removal",
    y = "Module (ASV / molecular patterns)") + plot_theme
ggsave(file.path(figures_dir, paste0(figure_prefix, "_groups.pdf")), summary_plot,
  width = 7.8, height = 7.2, bg = "white")
# The matrix has tens of thousands of cells; raster avoids a cumbersome PDF.
row_order <- active_rows[order(observed$fit$row[active_rows], -rowSums(web)[active_rows])]
column_order <- active_columns[order(observed$fit$column[active_columns], -colSums(web)[active_columns])]
heat <- expand_grid(row = seq_along(row_order), column = seq_along(column_order)) %>%
  mutate(phi = associations[cbind(row_order[row], column_order[column])],
    selected = phi >= edge_threshold)
heat_plot <- ggplot(heat, aes(column, row, fill = phi)) + geom_raster() +
  scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", midpoint = 0,
    limits = c(-1, 1), name = "Phi") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  labs(x = "Molecular patterns ordered by module", y = "ASV patterns ordered by module") +
  plot_theme + theme(panel.grid = element_blank())
ggsave(file.path(figures_dir, paste0(figure_prefix, "_association_matrix.png")), heat_plot,
  width = 7.2, height = 6.4, dpi = 600, bg = "white")

# 9. Verify input integrity and retain settings and software ----------------------
settings <- tibble(target_depth, minimum_present_and_absent, edge_threshold, n_null,
  observed_seed, null_seed, resampling_seed, max_sweeps, modularity_tolerance,
  quasiswap_thin, sites = nrow(site), segments = length(segments),
  asv_patterns = ncol(x), molecular_patterns = ncol(y), candidate_pairs = length(associations),
  optimizer_starts = paste(start_groups, collapse = ","), test_run)
write_csv(settings, file.path(results_dir, "settings.csv"))
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
write_csv(checks, file.path(results_dir, "input_integrity.csv"))
stopifnot(all(checks$unchanged), all(nulls$successful),
  all(bind_rows(resampling_rows)$converged), observed$statistics$converged)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
print(tests)
print(module_summary %>% filter(substantive_group))
