# Figure 2: first detection and matched-habitat iCAMP process importance.
# Adapt Ruiz-Gonzalez et al. (2015), doi:10.1111/ele.12499, to parallel
# aquatic compartments. Earliest detection is not a demonstrated source.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(permute)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
target_depth <- 10000
detection_threshold <- 1
balanced_draws <- 1000
balanced_samples_per_cell <- 3
point_size <- 3
point_alpha <- 0.5
panel_aspect <- 0.8
figure_width <- 9
single_row_height <- 4.8
composite_width <- 11.5
composite_height <- 10
base_seed <- 20261005
rank_seed <- 20261010L
rank_permutations <- 9999L
icamp_axis_max <- 60
icamp_panel_width <- 5.8
icamp_pair_point_size <- 1.6
icamp_pair_point_alpha <- 0.35
icamp_jitter_width <- 0.16
icamp_point_seed <- 20261014L
icamp_processes <- tibble::tribble(
  ~panel, ~process, ~label,
  "icamp_dispersal_limitation", "Dispersal.Limitation", "Dispersal limitation",
  "icamp_homogeneous_selection", "Homogeneous.Selection", "Homogeneous selection"
)
taxonomic_lineage <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rank_levels <- c("ASV", "Genus", "Family", "Order", "Class", "Phylum")
rank_scopes <- c("All reads", "Genus-classified ASVs")
unresolved_labels <- "^(uncultured|unclassified|unknown|unidentified|metagenome|Ambiguous_taxa)$"
rank_nulls <- tibble(null = c("All paired sites", "Within stream order"),
  seed = rank_seed + seq_len(2))
comparison_labels <- c(same_site = "Same site", other_site_mean = "Other sites (mean)")
habitats <- c("sediment", "hyporheic")
habitat_labels <- c(sediment = "Sediment", hyporheic = "Porewater")
source_pools <- c("Sediment", "Hyporheic")
network_groups <- c("Headwater", "Intermediate", "Mainstem")
category_labels <- c("Soil", "Sediment", "Hyporheic", "Sediment + Hyporheic")
category_display_labels <- c("Soil" = "Soil", "Sediment" = "Sediment",
  "Hyporheic" = "Porewater", "Sediment + Hyporheic" = "Sediment + Porewater")
category_colors <- c("Soil" = "#D55E00", "Sediment" = "#009E73",
  "Hyporheic" = "#0072B2", "Sediment + Hyporheic" = "#CC79A7")
results_dir <- here("results", "diversity_2016", "first_detection_flowpaths")
manuscript_dir <- here("figures", "manuscript")
supplement_dir <- here("figures", "supplement", "first_detection")
habitat_results_dir <- here("results", "diversity_2016", "soil_porewater_sediment")
icamp_results_dir <- here("results", "icamp_2016", "reviewed_pool_castor")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    strip.background = element_blank(), legend.position = "bottom",
    legend.title = element_text(size = 10), legend.text = element_text(size = 10))
theme_set(plot_theme)

# 2. Align fixed-depth samples and retain every identity in an audit -------------
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(here("data", "derived", "soil_porewater_sediment_2016", "asv_counts_10k.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sample_metadata.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sediment_porewater_pair_audit.csv"),
  file.path(habitat_results_dir, "soil_sharing_by_sample.csv"),
  file.path(habitat_results_dir, "paired_soil_sharing_test.csv"),
  here("results", "dada2_2016", "asv_taxonomy.csv"),
  file.path(icamp_results_dir, "matched_stream_samples.csv"),
  file.path(icamp_results_dir, "matched_within_group_pair_contrasts.csv"),
  file.path(icamp_results_dir, "matched_within_group_summary.csv"),
  file.path(icamp_results_dir, "pairwise_process_importance.csv"),
  file.path(icamp_results_dir, "input_checksums.csv"),
  file.path(icamp_results_dir, "icamp_settings.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
stopifnot(!anyNA(checks$initial_md5))
count_table <- read_csv(input_paths[1], show_col_types = FALSE)
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()), show_col_types = FALSE)
counts <- as.matrix(count_table[, -1])
rownames(counts) <- count_table$asv_id
metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
stopifnot(identical(metadata$sample_id, colnames(counts)), !anyDuplicated(metadata$sample_id),
  !anyDuplicated(rownames(counts)), all(metadata$included_manuscript),
  !any(metadata$excluded_suspect_library),
  all(colSums(counts) == target_depth), all(rowSums(counts) > 0),
  all(counts >= 0 & counts == floor(counts)), detection_threshold == 1)
metadata <- metadata %>% mutate(network_group = case_when(
  habitat == "soil" ~ "Soil", site_stream_order %in% 1:2 ~ "Headwater",
  site_stream_order %in% 3:4 ~ "Intermediate", site_stream_order == 5 ~ "Mainstem"),
  source_pool = case_when(habitat == "soil" ~ "Soil", habitat == "sediment" ~ "Sediment",
    habitat == "hyporheic" ~ "Hyporheic"),
  cell = paste(habitat, network_group, sep = "_"))
stopifnot(nrow(metadata) == 76, sum(metadata$habitat == "soil") == 15,
  sum(metadata$habitat == "sediment") == 34, sum(metadata$habitat == "hyporheic") == 27,
  !anyNA(metadata$network_group), all(metadata$habitat %in% c("soil", habitats)))
cell_design <- metadata %>% group_by(cell, habitat, source_pool, network_group) %>%
  summarize(samples = n(), segments = n_distinct(site_segment_number, na.rm = TRUE), .groups = "drop")
stopifnot(nrow(cell_design) == 7, min(cell_design$samples) >= balanced_samples_per_cell)
presence <- counts >= detection_threshold
cell_indices <- split(seq_len(nrow(metadata)), metadata$cell)

# 3. Repeat the detection inventory with all samples and with balanced sampling --
# Soils have precedence. For remaining taxa, use the earliest network group
# containing a detection. Sediment and hyporheic porewater are parallel pools,
# with no assumed flow order between them; earliest-stage ties remain shared.
# Equal-sample draws retain three profiles per habitat-stage cell and three
# soils, giving 21 libraries per draw. Equal sample counts do not establish
# equal detection completeness. Receiving habitats remain separate. Each draw
# rebuilds the inventory and evaluates its selected communities; this is a
# sampling sensitivity, not a confidence interval or independent replication.
draw_summary_rows <- vector("list", balanced_draws + 1)
selection_rows <- vector("list", balanced_draws)
draw_diagnostics <- vector("list", balanced_draws + 1)
for (draw_id in 0:balanced_draws) {
  draw_start <- proc.time()[["elapsed"]]
  if (draw_id == 0) {
    selected <- seq_len(nrow(metadata))
  } else {
    set.seed(base_seed + draw_id)
    selected <- sort(unlist(lapply(cell_indices, sample, size = balanced_samples_per_cell), use.names = FALSE))
    selection_rows[[draw_id]] <- metadata[selected, ] %>% select(sample_id, habitat, source_pool, network_group, cell) %>%
      mutate(draw_id = draw_id, seed = base_seed + draw_id, .before = 1)
  }
  first_category <- integer(nrow(counts))
  first_stage <- rep(NA_integer_, nrow(counts))
  soil_columns <- selected[metadata$habitat[selected] == "soil"]
  soil_detected <- rowSums(presence[, soil_columns, drop = FALSE]) > 0
  first_category[soil_detected] <- 1L
  first_stage[soil_detected] <- 0L
  for (stage_id in seq_along(network_groups)) {
    stage_presence <- matrix(FALSE, nrow(counts), length(source_pools))
    for (pool_id in seq_along(source_pools)) {
      columns <- selected[metadata$source_pool[selected] == source_pools[pool_id] &
        metadata$network_group[selected] == network_groups[stage_id]]
      stage_presence[, pool_id] <- rowSums(presence[, columns, drop = FALSE]) > 0
    }
    new_detection <- first_category == 0 & rowSums(stage_presence) > 0
    single_pool <- new_detection & rowSums(stage_presence) == 1
    first_category[new_detection] <- 4L
    first_category[single_pool] <- max.col(stage_presence[single_pool, , drop = FALSE], ties.method = "first") + 1L
    first_stage[new_detection] <- stage_id
  }
  aquatic_columns <- selected[metadata$habitat[selected] != "soil"]
  observed_asvs <- colSums(presence[, aquatic_columns, drop = FALSE])
  profiles <- vector("list", length(category_labels))
  for (category_id in seq_along(category_labels)) {
    rows <- first_category == category_id
    profiles[[category_id]] <- metadata[aquatic_columns, ] %>%
      select(sample_id, site_code, habitat, network_group, cell) %>%
      mutate(category = category_labels[category_id],
        asvs_detected = colSums(presence[rows, aquatic_columns, drop = FALSE]),
        observed_asvs = observed_asvs, richness_percent = 100 * asvs_detected / observed_asvs,
        reads = colSums(counts[rows, aquatic_columns, drop = FALSE]),
        read_percent = 100 * reads / target_depth)
  }
  profiles <- bind_rows(profiles)
  totals <- profiles %>% group_by(sample_id) %>% summarize(
    richness_percent = sum(richness_percent), read_percent = sum(read_percent), .groups = "drop")
  stopifnot(max(abs(totals$richness_percent - 100)) < 1e-8,
    max(abs(totals$read_percent - 100)) < 1e-8,
    all(first_category[rowSums(presence[, selected, drop = FALSE]) > 0] > 0))
  draw_summary_rows[[draw_id + 1]] <- profiles %>% group_by(habitat, network_group, category) %>%
    summarize(samples = n(), richness_percent = mean(richness_percent),
      read_percent = mean(read_percent), .groups = "drop") %>% mutate(draw_id = draw_id, .before = 1)
  if (draw_id == 0) {
    write_csv(profiles, file.path(results_dir, "sample_first_detection_profiles.csv"))
    assignments <- tibble(asv_id = rownames(counts),
      first_detection_pool = category_labels[first_category],
      earliest_network_group = c("Soil", network_groups)[first_stage + 1L])
    stopifnot(!anyNA(assignments$first_detection_pool))
    write_csv(assignments, file.path(results_dir, "asv_first_detection_assignments.csv"))
    # Track first-stage additions separately from the habitat-origin display.
    stage_profiles <- vector("list", 4)
    for (stage_id in 0:3) {
      rows <- first_stage == stage_id
      stage_profiles[[stage_id + 1]] <- metadata[aquatic_columns, ] %>% select(sample_id, habitat, network_group) %>%
        mutate(first_network_group = c("Soil", network_groups)[stage_id + 1L],
          richness_percent = 100 * colSums(presence[rows, aquatic_columns, drop = FALSE]) / observed_asvs)
    }
    write_csv(bind_rows(stage_profiles), file.path(results_dir, "sample_first_stage_profiles.csv"))
  }
  draw_diagnostics[[draw_id + 1]] <- tibble(draw_id = draw_id,
    seed = ifelse(draw_id == 0, NA_integer_, base_seed + draw_id), samples = length(selected),
    assigned_asvs = sum(first_category > 0), seconds = proc.time()[["elapsed"]] - draw_start,
    diagnostic = "")
  if (draw_id %in% c(0, 1, balanced_draws) || draw_id %% 100 == 0) {
    message("Completed detection inventory ", draw_id, "/", balanced_draws,
      "; elapsed for inventory: ", round(draw_diagnostics[[draw_id + 1]]$seconds, 2), " s")
  }
}
summaries <- bind_rows(draw_summary_rows)
full_summary <- summaries %>% filter(draw_id == 0)
balanced_summary <- summaries %>% filter(draw_id > 0) %>% group_by(habitat, network_group, category) %>%
  summarize(draws = n(), mean_richness_percent = mean(richness_percent),
    richness_p025 = quantile(richness_percent, 0.025), richness_p975 = quantile(richness_percent, 0.975),
    mean_read_percent = mean(read_percent), read_p025 = quantile(read_percent, 0.025),
    read_p975 = quantile(read_percent, 0.975), .groups = "drop")
write_csv(full_summary, file.path(results_dir, "full_inventory_network_summary.csv"))
write_csv(balanced_summary, file.path(results_dir, "balanced_inventory_summary.csv"))
write_csv(summaries %>% filter(draw_id > 0), file.path(results_dir, "balanced_inventory_draws.csv"))
write_csv(bind_rows(selection_rows), file.path(results_dir, "balanced_sample_selections.csv"))
write_csv(bind_rows(draw_diagnostics), file.path(results_dir, "inventory_diagnostics.csv"))

# 4. Compare paired habitats across nested taxonomic resolutions ----------------
pair_audit <- read_csv(input_paths[3], col_types = cols(site_code = col_character())) %>%
  filter(paired_manuscript) %>% arrange(site_code)
stopifnot(nrow(pair_audit) == 18, !anyDuplicated(pair_audit$site_code),
  all(c(pair_audit$sediment_id, pair_audit$porewater_id) %in% colnames(counts)),
  all(metadata$site_code[match(pair_audit$sediment_id, metadata$sample_id)] == pair_audit$site_code),
  all(metadata$site_code[match(pair_audit$porewater_id, metadata$sample_id)] == pair_audit$site_code),
  all(metadata$habitat[match(pair_audit$sediment_id, metadata$sample_id)] == "sediment"),
  all(metadata$habitat[match(pair_audit$porewater_id, metadata$sample_id)] == "hyporheic"))
paired_counts <- counts[, c(pair_audit$sediment_id, pair_audit$porewater_id), drop = FALSE]
active <- rowSums(paired_counts) > 0
paired_counts <- paired_counts[active, , drop = FALSE]
taxonomy <- read_csv(input_paths[6], show_col_types = FALSE)
stopifnot(!anyDuplicated(taxonomy$asv_id))
taxonomy <- taxonomy[match(rownames(paired_counts), taxonomy$asv_id), ]
stopifnot(identical(taxonomy$asv_id, rownames(paired_counts)))
for (rank in taxonomic_lineage) {
  taxonomy[[rank]] <- trimws(taxonomy[[rank]])
  unresolved <- is.na(taxonomy[[rank]]) | taxonomy[[rank]] == "" |
    grepl(unresolved_labels, taxonomy[[rank]], ignore.case = TRUE)
  taxonomy[[rank]][unresolved] <- NA_character_
}
# Reuse identical intact-profile permutations at every rank and in both scopes.
orders <- metadata$site_stream_order[match(pair_audit$porewater_id, metadata$sample_id)]
stopifnot(!anyNA(orders),
  all(orders == metadata$site_stream_order[match(pair_audit$sediment_id, metadata$sample_id)]))
rank_shuffles <- list()
for (i in seq_len(nrow(rank_nulls))) {
  set.seed(rank_nulls$seed[i])
  control <- how()
  if (rank_nulls$null[i] == "Within stream order") control <- how(blocks = factor(orders))
  rank_shuffles[[rank_nulls$null[i]]] <- shuffleSet(nrow(pair_audit),
    nset = rank_permutations, control = control)
}
# The sensitivity fixes the same genus-classified ASVs across all ranks and
# renormalizes each sample. Rank-dependent filtering would change composition.
genus_classified <- complete.cases(taxonomy[, taxonomic_lineage])
rank_sample_rows <- list()
rank_test_rows <- list()
coverage_rows <- list()
for (scope in rank_scopes) {
  retained <- if (scope == "All reads") rep(TRUE, nrow(paired_counts)) else genus_classified
  current_counts <- paired_counts[retained, , drop = FALSE]
  current_taxonomy <- taxonomy[retained, ]
  stopifnot(all(colSums(current_counts) > 0))
  previous_distance <- NULL
  for (rank in rank_levels) {
    if (rank == "ASV") {
      assigned <- rep(TRUE, nrow(current_counts))
      aggregated <- current_counts
    } else {
      prefix <- taxonomic_lineage[seq_len(match(rank, taxonomic_lineage))]
      assigned <- complete.cases(current_taxonomy[, prefix])
      keys <- apply(current_taxonomy[, prefix], 1, paste, collapse = "|")
      # Missing ranks stay ASV-specific, preserving reads without creating a
      # shared "unclassified" taxon. Complete lineage prefixes avoid homonyms.
      keys[!assigned] <- paste0("unresolved_ASV|", current_taxonomy$asv_id[!assigned])
      aggregated <- rowsum(current_counts, keys, reorder = FALSE)
    }
    stopifnot(identical(unname(colSums(aggregated)), unname(colSums(current_counts))))
    distance <- as.matrix(dist(decostand(t(aggregated), method = "hellinger")))
    stopifnot(all(is.finite(distance)), all(distance >= 0), all(diag(distance) == 0),
      max(abs(distance - t(distance))) < 1e-12, max(distance) <= sqrt(2) + 1e-12)
    # Nested aggregation must not increase Hellinger distance. Its decline
    # alone is expected mathematically and does not demonstrate convergence.
    if (!is.null(previous_distance)) stopifnot(all(distance <= previous_distance + 1e-12))
    previous_distance <- distance
    cross <- distance[pair_audit$sediment_id, pair_audit$porewater_id, drop = FALSE]
    observed <- mean(diag(cross))
    rank_sample_rows[[length(rank_sample_rows) + 1L]] <- pair_audit %>%
      select(site_code, sediment_id, porewater_id) %>%
      mutate(scope, rank, metric = "Hellinger", same_site = diag(cross),
        other_site_mean = (colSums(cross) - diag(cross)) / (nrow(pair_audit) - 1L))
    for (i in seq_len(nrow(rank_nulls))) {
      shuffles <- rank_shuffles[[rank_nulls$null[i]]]
      null_values <- apply(shuffles, 1, function(p) mean(cross[cbind(seq_len(nrow(pair_audit)), p)]))
      rank_test_rows[[length(rank_test_rows) + 1L]] <- tibble(scope, rank,
        metric = "Hellinger", null = rank_nulls$null[i], pairs = nrow(pair_audit),
        observed_mean = observed, observed_median = median(diag(cross)),
        offsite_mean = mean(cross[row(cross) != col(cross)]),
        random_assignment_mean = mean(null_values),
        null_p025 = unname(quantile(null_values, 0.025)),
        null_p975 = unname(quantile(null_values, 0.975)),
        p = (1L + sum(null_values <= observed)) / (1L + length(null_values)),
        permutations = length(null_values), seed = rank_nulls$seed[i],
        alternative = "Lower same-site distance")
    }
    if (scope == "All reads") {
      coverage_rows[[length(coverage_rows) + 1L]] <- metadata[
        match(colnames(paired_counts), metadata$sample_id), ] %>%
        select(sample_id, site_code, habitat) %>% mutate(rank,
          classified_reads = colSums(current_counts[assigned, , drop = FALSE]),
          total_reads = colSums(current_counts),
          classified_read_percent = 100 * classified_reads / total_reads)
    }
  }
}
rank_samples <- bind_rows(rank_sample_rows)
rank_tests <- bind_rows(rank_test_rows) %>% group_by(scope, null) %>%
  mutate(p_bh_family6 = p.adjust(p, "BH")) %>% ungroup()
write_csv(rank_samples, file.path(results_dir, "taxonomic_rank_sample_comparisons.csv"))
write_csv(rank_tests, file.path(results_dir, "taxonomic_rank_tests.csv"))
write_csv(bind_rows(coverage_rows), file.path(results_dir, "taxonomic_rank_coverage.csv"))

# 5. Validate completed iCAMP summaries on the same 18 paired sites --------------
# Reuse assignments from the reviewed regional pool. Both habitats contribute
# the same unordered within-position site pairs, including same-segment pairs.
# These rows share sites; means are descriptive, without independent-pair tests.
icamp_samples <- read_csv(input_paths[7], col_types = cols(site_code = col_character()))
icamp_pairs <- read_csv(input_paths[8],
  col_types = cols(site_low = col_character(), site_high = col_character()))
icamp_summary <- read_csv(input_paths[9], show_col_types = FALSE)
icamp_assignments <- read_csv(input_paths[10], show_col_types = FALSE)
icamp_checks <- read_csv(input_paths[11], show_col_types = FALSE)
icamp_settings <- read_csv(input_paths[12], show_col_types = FALSE)
stopifnot(all(unname(tools::md5sum(here(icamp_checks$input))) == icamp_checks$md5),
  icamp_settings$value[icamp_settings$setting == "samples"] == "76",
  icamp_settings$value[icamp_settings$setting == "asvs"] == "33974",
  icamp_settings$value[icamp_settings$setting == "randomizations"] == "1000")
paired_manifest <- metadata %>%
  filter(sample_id %in% c(pair_audit$sediment_id, pair_audit$porewater_id)) %>%
  select(sample_id, site_code, habitat, network_group, site_segment_number, site_stream_order) %>%
  arrange(sample_id)
stopifnot(identical(as.data.frame(icamp_samples %>% arrange(sample_id)),
  as.data.frame(paired_manifest)), nrow(icamp_samples) == 36L,
  !anyDuplicated(icamp_samples$sample_id), !anyNA(icamp_pairs))
icamp_group_sizes <- icamp_samples %>% distinct(site_code, network_group) %>%
  count(network_group, name = "sites") %>%
  mutate(network_group = factor(network_group, levels = network_groups)) %>%
  arrange(network_group) %>% mutate(pairs = choose(sites, 2))
stopifnot(identical(icamp_group_sizes$sites, c(11L, 3L, 4L)),
  nrow(icamp_pairs) == sum(icamp_group_sizes$pairs) * 5L,
  !anyDuplicated(icamp_pairs[, c("site_low", "site_high", "process")]),
  all(icamp_pairs$Sediment >= 0 & icamp_pairs$Sediment <= 1),
  all(icamp_pairs$Porewater >= 0 & icamp_pairs$Porewater <= 1),
  max(abs(icamp_pairs$porewater_minus_sediment -
    (icamp_pairs$Porewater - icamp_pairs$Sediment))) < 1e-12)
# Check every plotted source pair against its canonical pooled assignment.
assignment_keys <- paste(pmin(icamp_assignments$sample1, icamp_assignments$sample2),
  pmax(icamp_assignments$sample1, icamp_assignments$sample2), sep = "|")
stopifnot(!anyDuplicated(assignment_keys), nrow(icamp_assignments) == choose(76, 2))
assignment_values <- as.matrix(icamp_assignments[, unique(icamp_pairs$process)])
for (habitat in habitats) {
  samples <- icamp_samples %>% filter(.data$habitat == .env$habitat)
  first <- samples$sample_id[match(icamp_pairs$site_low, samples$site_code)]
  second <- samples$sample_id[match(icamp_pairs$site_high, samples$site_code)]
  pair_rows <- match(paste(pmin(first, second), pmax(first, second), sep = "|"), assignment_keys)
  expected <- assignment_values[cbind(pair_rows, match(icamp_pairs$process, colnames(assignment_values)))]
  observed <- icamp_pairs[[habitat_labels[habitat]]]
  stopifnot(!anyNA(expected), max(abs(expected - observed)) < 1e-12,
    all(samples$network_group[match(icamp_pairs$site_low, samples$site_code)] == icamp_pairs$network_group),
    all(samples$network_group[match(icamp_pairs$site_high, samples$site_code)] == icamp_pairs$network_group))
}
icamp_expected <- icamp_pairs %>% group_by(network_group, process) %>%
  summarize(pairs = n(), sediment_mean = mean(Sediment), porewater_mean = mean(Porewater),
    mean_difference = mean(porewater_minus_sediment), .groups = "drop") %>%
  arrange(network_group, process)
icamp_summary <- icamp_summary %>% mutate(pairs = as.integer(pairs)) %>% arrange(network_group, process)
stopifnot(identical(icamp_summary[, c("network_group", "process", "pairs")],
  icamp_expected[, c("network_group", "process", "pairs")]),
  max(abs(as.matrix(icamp_summary[, c("sediment_mean", "porewater_mean", "mean_difference")]) -
    as.matrix(icamp_expected[, c("sediment_mean", "porewater_mean", "mean_difference")]))) < 1e-12)
icamp_plot_data <- icamp_summary %>% filter(process %in% icamp_processes$process) %>%
  pivot_longer(c(sediment_mean, porewater_mean), names_to = "habitat", values_to = "importance") %>%
  mutate(habitat = if_else(habitat == "sediment_mean", "sediment", "hyporheic"),
    habitat = factor(habitat, levels = habitats),
    network_group = factor(network_group, levels = network_groups),
    plot_x = as.integer(network_group))
stopifnot(nrow(icamp_plot_data) == 12L, all(100 * icamp_plot_data$importance <= icamp_axis_max))
# Give each unordered site pair one horizontal offset, shared by both habitats
# and both processes. Vertical coordinates remain the observed assignments.
icamp_pair_positions <- icamp_pairs %>% distinct(site_low, site_high, network_group) %>%
  mutate(network_group = factor(network_group, levels = network_groups)) %>%
  arrange(network_group, site_low, site_high)
set.seed(icamp_point_seed)
icamp_pair_positions <- icamp_pair_positions %>%
  mutate(plot_x = as.integer(network_group) + runif(n(), -icamp_jitter_width, icamp_jitter_width))
icamp_pair_points <- icamp_pairs %>% filter(process %in% icamp_processes$process) %>%
  select(site_low, site_high, network_group, process, Sediment, Porewater) %>%
  mutate(network_group = factor(network_group, levels = network_groups)) %>%
  left_join(icamp_pair_positions, by = c("site_low", "site_high", "network_group")) %>%
  pivot_longer(c(Sediment, Porewater), names_to = "habitat", values_to = "importance") %>%
  mutate(habitat = if_else(habitat == "Sediment", "sediment", "hyporheic"),
    habitat = factor(habitat, levels = habitats))
stopifnot(nrow(icamp_pair_points) == sum(icamp_group_sizes$pairs) * 4L,
  !anyNA(icamp_pair_points), all(100 * icamp_pair_points$importance <= icamp_axis_max))

# 6. Show inventories and iCAMP contrasts; retain earlier panels as supplements ---
plot_data <- full_summary %>% mutate(category = factor(category, levels = category_labels),
  habitat = factor(habitat, levels = habitats), network_group = factor(network_group, levels = network_groups))
plots <- list()
for (metric in c("richness_percent", "read_percent")) {
  y_label <- ifelse(metric == "richness_percent", "ASV richness (%)", "ASV reads (%)")
  plots[[metric]] <- ggplot(plot_data, aes(network_group, .data[[metric]], fill = category)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    geom_text(data = distinct(plot_data, habitat, network_group, samples),
      aes(network_group, 104, label = paste0("n = ", samples)), inherit.aes = FALSE, size = 3) +
    facet_wrap(~habitat, nrow = 1, labeller = as_labeller(habitat_labels)) +
    scale_fill_manual(values = category_colors, labels = category_display_labels,
      drop = FALSE, name = "First detection pool",
      guide = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = c("Headwater\n(order 1-2)", "Intermediate\n(order 3-4)", "Mainstem\n(order 5)")) +
    labs(x = NULL, y = y_label) + theme(aspect.ratio = panel_aspect)
  ggsave(file.path(supplement_dir, paste0("2016_first_detection_", metric, ".pdf")),
    plots[[metric]], width = figure_width, height = single_row_height, bg = "white")
}
# Reuse script 02's formal paired soil-sharing test and verify it from counts.
sharing <- read_csv(input_paths[4], col_types = cols(site_code = col_character()))
soil_test <- read_csv(input_paths[5], show_col_types = FALSE)
stopifnot(!anyDuplicated(sharing$sample_id), soil_test$pairs == 18,
  all(c(pair_audit$sediment_id, pair_audit$porewater_id) %in% sharing$sample_id))
soil_inventory <- rowSums(presence[, metadata$habitat == "soil", drop = FALSE]) > 0
expected_sharing <- 100 * colSums(counts[soil_inventory, sharing$sample_id, drop = FALSE]) / target_depth
stopifnot(max(abs(expected_sharing - sharing$soil_read_percent)) < 1e-10)
paired_sharing <- pair_audit %>% transmute(site_code,
  sediment = sharing$soil_read_percent[match(sediment_id, sharing$sample_id)],
  hyporheic = sharing$soil_read_percent[match(porewater_id, sharing$sample_id)])
stopifnot(abs(mean(paired_sharing$sediment - paired_sharing$hyporheic) -
  soil_test$mean_difference_percentage_points) < 1e-10)
paired_sharing <- paired_sharing %>%
  pivot_longer(c(sediment, hyporheic), names_to = "habitat", values_to = "read_percent") %>%
  mutate(habitat = factor(habitat, levels = habitats))
plots$paired_soil_sharing <- ggplot(paired_sharing, aes(habitat, read_percent, group = site_code)) +
  geom_line(color = "gray70", linewidth = 0.4) +
  geom_point(aes(color = habitat), size = point_size, alpha = point_alpha) +
  scale_color_manual(values = c(sediment = "#009E73", hyporheic = "#0072B2")) +
  scale_x_discrete(labels = habitat_labels) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Reads shared with soil inventory (%)") +
  plot_theme + theme(legend.position = "none", aspect.ratio = panel_aspect)
rank_plot_data <- rank_samples %>%
  pivot_longer(c(same_site, other_site_mean), names_to = "comparison", values_to = "distance") %>%
  mutate(rank = factor(rank, levels = rank_levels), scope = factor(scope, levels = rank_scopes),
    comparison = factor(comparison, levels = names(comparison_labels)))
rank_means <- rank_plot_data %>% group_by(scope, rank, comparison) %>%
  summarize(distance = mean(distance), .groups = "drop")
for (view in c("taxonomic_resolution", "taxonomic_resolution_sensitivity")) {
  profiles <- rank_plot_data %>% filter(comparison == "same_site")
  means <- rank_means
  if (view == "taxonomic_resolution") {
    # The two mean curves nearly overlap. Show the same-site mean alone in
    # the supplementary panel; retain both comparisons in the sensitivity
    # plot and tables.
    profiles <- profiles %>% filter(scope == "All reads")
    means <- means %>% filter(scope == "All reads", comparison == "same_site")
  }
  plots[[view]] <- ggplot(profiles, aes(rank, distance, group = site_code)) +
    geom_line(color = "gray80", linewidth = 0.35) +
    geom_line(data = means, aes(group = comparison, linetype = comparison),
      color = "black", linewidth = 0.8) +
    scale_linetype_manual(values = c(same_site = "solid", other_site_mean = "dashed"),
      labels = comparison_labels, name = NULL) +
    scale_y_continuous(limits = c(0, sqrt(2)), breaks = c(0, 0.4, 0.8, 1.2)) +
    labs(x = "Taxonomic resolution", y = "Sediment-porewater Hellinger distance") +
    plot_theme + theme(aspect.ratio = panel_aspect, legend.text = element_text(size = 9),
      legend.key.width = grid::unit(1.2, "lines"))
  if (view == "taxonomic_resolution") {
    plots[[view]] <- plots[[view]] + theme(legend.position = "none")
  } else {
    plots[[view]] <- plots[[view]] + facet_wrap(~scope)
  }
}
for (panel in c("paired_soil_sharing", "taxonomic_resolution", "taxonomic_resolution_sensitivity")) {
  ggsave(file.path(supplement_dir, paste0("2016_", panel, ".pdf")), plots[[panel]],
    width = ifelse(panel == "taxonomic_resolution_sensitivity", 11.5, 5.8),
    height = single_row_height, bg = "white")
}
for (i in seq_len(nrow(icamp_processes))) {
  panel <- icamp_plot_data %>% filter(process == icamp_processes$process[i])
  pair_points <- icamp_pair_points %>% filter(process == icamp_processes$process[i])
  # Lines connect descriptive group means; they are not fitted regressions.
  plots[[icamp_processes$panel[i]]] <- ggplot(panel,
    aes(plot_x, 100 * importance, color = habitat, shape = habitat, group = habitat)) +
    geom_point(data = pair_points, size = icamp_pair_point_size,
      alpha = icamp_pair_point_alpha, show.legend = FALSE) +
    geom_line(linewidth = 0.6) + geom_point(size = point_size) +
    scale_color_manual(values = setNames(category_colors[c("Sediment", "Hyporheic")], habitats),
      labels = habitat_labels, name = "Habitat") +
    scale_shape_manual(values = c(sediment = 16, hyporheic = 17),
      labels = habitat_labels, name = "Habitat") +
    scale_x_continuous(breaks = seq_along(network_groups), limits = c(0.5, 3.5),
      labels = paste0(network_groups, "\n", icamp_group_sizes$sites, " paired sites")) +
    scale_y_continuous(limits = c(0, icamp_axis_max), breaks = seq(0, icamp_axis_max, 20),
      expand = expansion(mult = c(0.02, 0.04))) +
    facet_wrap(~process, labeller = as_labeller(setNames(icamp_processes$label, icamp_processes$process))) +
    labs(x = NULL, y = "Pairwise importance (%)") + plot_theme +
    theme(aspect.ratio = panel_aspect, axis.text.x = element_text(size = 9))
  ggsave(file.path(supplement_dir, paste0("2016_first_detection_", icamp_processes$panel[i], ".pdf")),
    plots[[icamp_processes$panel[i]]], width = icamp_panel_width, height = single_row_height, bg = "white")
}
inventory_panel <- plots$read_percent +
  geom_text(data = tibble(habitat = factor(habitats, levels = habitats), tag = LETTERS[1:2]),
    aes(x = -Inf, y = Inf, label = tag), inherit.aes = FALSE,
    hjust = -0.5, vjust = 1.3, fontface = "bold", size = 4)
icamp_panels <- ((plots$icamp_dispersal_limitation + labs(tag = "C")) |
  (plots$icamp_homogeneous_selection + labs(tag = "D"))) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined <- inventory_panel / icamp_panels
ggsave(file.path(manuscript_dir, "2016_first_detection_flowpaths.pdf"), combined,
  width = composite_width, height = composite_height, bg = "white")
if (interactive()) {
  print(combined)
}

# 7. Verify that the community and completed iCAMP inputs were not changed -------
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "input_checksums.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
print(full_summary %>% filter(category == "Soil"))
print(balanced_summary %>% filter(category == "Soil"))
