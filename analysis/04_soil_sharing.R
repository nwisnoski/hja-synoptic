# Figures 3 and 4: taxa that sediment and porewater share with soil.
# Figure 3 adapts the first-detection inventory of Ruiz-Gonzalez et al. (2015),
# doi:10.1111/ele.12499, to parallel aquatic compartments; earliest detection is
# not a demonstrated source. Figure 4 asks where the soil-core share of sediment
# reads is high, whether porewater carries the same share, and which taxa the
# soil-core ASVs are. ASVs found in at least five of the 15 soils are soil-core,
# defined from the soil samples alone. Run 02_habitat_diversity.R first.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(permute)
library(igraph)
library(spdep)
library(ggplot2)
library(patchwork)
source(here("analysis", "figure_style.R"))

# 1. Settings and plotting theme ------------------------------------------------
target_depth <- 10000L
detection_threshold <- 1
# An ASV is soil-core if it occurs in at least soil_core_min of the 15 soils. The
# class uses soil samples only, so the stream data it is tested on play no part
# in defining it. In Figure 3 the remaining soil ASVs, found in fewer soils, are
# the part of the soil share that grows with the number of soils sampled.
soil_core_min <- 5L
# First-detection inventory (Figure 3).
balanced_draws <- 1000
balanced_samples_per_cell <- 3
point_alpha <- 0.5
detection_seed <- 20261005
rank_seed <- 20261010L
rank_permutations <- 9999L
taxonomic_lineage <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rank_levels <- c("ASV", "Genus", "Family", "Order", "Class", "Phylum")
rank_scopes <- c("All reads", "Genus-classified ASVs")
unresolved_labels <- "^(uncultured|unclassified|unknown|unidentified|metagenome|Ambiguous_taxa)$"
rank_nulls <- tibble(null = c("All paired sites", "Within stream order"),
  seed = rank_seed + seq_len(2))
comparison_labels <- c(same_site = "Same site", other_site_mean = "Other sites (mean)")
habitats <- c("sediment", "hyporheic")
source_pools <- c("Sediment", "Hyporheic")
network_groups <- c("Headwater", "Intermediate", "Mainstem")
category_labels <- c("Soil", "Sediment", "Hyporheic", "Sediment + Hyporheic")
display_levels <- c("Soil core", "Soil non-core", "Sediment", "Hyporheic", "Sediment + Hyporheic")
display_labels <- c("Soil core" = "Soil: core (5+ sites)", "Soil non-core" = "Soil: non-core (<5 sites)",
  "Sediment" = "Sediment", "Hyporheic" = "Porewater", "Sediment + Hyporheic" = "Sediment + Porewater")
display_colors <- c("Soil core" = "#D55E00", "Soil non-core" = "#E69F00", "Sediment" = "#009E73",
  "Hyporheic" = "#0072B2", "Sediment + Hyporheic" = "#CC79A7")
# Soil-core share across the network (Figure 4).
n_permutations <- 9999L
soil_core_seed <- 20261116L
significance_threshold <- 0.05
# The major fork of the network: segment 134 drains the north branch and
# segment 39 the upper Lookout branch; they meet at 7.7 km from the outlet.
north_segment <- 134
lookout_segment <- 39
distance_tolerance_m <- 1
# Number of soil-core lineages shown in the ranked panel: the most abundant in
# soil, so the panel shows the lineages that dominate soils rather than taxa
# that are rare in soil and common elsewhere.
n_ranked_lineages <- 25L
detection_results_dir <- here("results", "diversity_2016", "first_detection_flowpaths")
soil_core_results_dir <- here("results", "diversity_2016", "soil_affinity")
detection_figures_dir <- here("figures", "supplement", "first_detection")
soil_core_figures_dir <- here("figures", "supplement", "soil_affinity")
manuscript_dir <- here("figures", "manuscript")
habitat_results_dir <- here("results", "diversity_2016", "soil_porewater_sediment")
# Habitat colors and shapes match Figure 1. Text sizes, panel sizes, and the
# exporter come from figure_style.R.
habitat_labels <- c(soil = "Soil", sediment = "Sediment", hyporheic = "Porewater")
habitat_colors <- c(soil = "#D55E00", sediment = "#009E73", hyporheic = "#0072B2")
habitat_shapes <- c(soil = 17, sediment = 16, hyporheic = 15)
map_point_size <- 1.5 * figure_point_size
plot_theme <- figure_theme
theme_set(plot_theme)

# 2. Align fixed-depth samples and retain every identity in an audit -------------
for (directory in c(manuscript_dir, detection_figures_dir, soil_core_figures_dir,
    detection_results_dir, soil_core_results_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
input_paths <- c(here("data", "derived", "soil_porewater_sediment_2016", "asv_counts_10k.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sample_metadata.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sediment_porewater_pair_audit.csv"),
  file.path(habitat_results_dir, "soil_sharing_by_sample.csv"),
  file.path(habitat_results_dir, "paired_soil_sharing_test.csv"),
  here("results", "dada2_2016", "asv_taxonomy.csv"),
  file.path(habitat_results_dir, "soil_inventory.csv"),
  here("data", "StreamCenterline.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_metadata_screened.csv"))
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

# The 18 paired sites in file order, as used by the soil-core tests.
pairs <- read_csv(input_paths[3], col_types = cols(site_code = col_character()), show_col_types = FALSE) %>%
  filter(paired_manuscript)
stopifnot(nrow(pairs) == 18)

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
    set.seed(detection_seed + draw_id)
    selected <- sort(unlist(lapply(cell_indices, sample, size = balanced_samples_per_cell), use.names = FALSE))
    selection_rows[[draw_id]] <- metadata[selected, ] %>% select(sample_id, habitat, source_pool, network_group, cell) %>%
      mutate(draw_id = draw_id, seed = detection_seed + draw_id, .before = 1)
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
    write_csv(profiles, file.path(detection_results_dir, "sample_first_detection_profiles.csv"))
    assignments <- tibble(asv_id = rownames(counts),
      first_detection_pool = category_labels[first_category],
      earliest_network_group = c("Soil", network_groups)[first_stage + 1L])
    stopifnot(!anyNA(assignments$first_detection_pool))
    write_csv(assignments, file.path(detection_results_dir, "asv_first_detection_assignments.csv"))
    # Track first-stage additions separately from the habitat-origin display.
    stage_profiles <- vector("list", 4)
    for (stage_id in 0:3) {
      rows <- first_stage == stage_id
      stage_profiles[[stage_id + 1]] <- metadata[aquatic_columns, ] %>% select(sample_id, habitat, network_group) %>%
        mutate(first_network_group = c("Soil", network_groups)[stage_id + 1L],
          richness_percent = 100 * colSums(presence[rows, aquatic_columns, drop = FALSE]) / observed_asvs)
    }
    write_csv(bind_rows(stage_profiles), file.path(detection_results_dir, "sample_first_stage_profiles.csv"))
  }
  draw_diagnostics[[draw_id + 1]] <- tibble(draw_id = draw_id,
    seed = ifelse(draw_id == 0, NA_integer_, detection_seed + draw_id), samples = length(selected),
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
write_csv(full_summary, file.path(detection_results_dir, "full_inventory_network_summary.csv"))
write_csv(balanced_summary, file.path(detection_results_dir, "balanced_inventory_summary.csv"))
write_csv(summaries %>% filter(draw_id > 0), file.path(detection_results_dir, "balanced_inventory_draws.csv"))
write_csv(bind_rows(selection_rows), file.path(detection_results_dir, "balanced_sample_selections.csv"))
write_csv(bind_rows(draw_diagnostics), file.path(detection_results_dir, "inventory_diagnostics.csv"))

# 4. Identify the soil-core ASVs and split the soil first detections ------------
soil_ids <- metadata$sample_id[metadata$habitat == "soil"]
soil_samples_detected <- rowSums(counts[, soil_ids] > 0)
# Script 02 records the same soil occurrence; the two must agree.
saved_soil_inventory <- read_csv(input_paths[7], show_col_types = FALSE)
stopifnot(identical(unname(soil_samples_detected[saved_soil_inventory$asv_id]),
  saved_soil_inventory$soil_samples_detected))
sediment <- metadata %>% filter(habitat == "sediment")
sediment_counts <- counts[, sediment$sample_id]
sediment_prevalence <- rowSums(sediment_counts > 0)
# Every soil-core ASV counts toward the shares in both habitats, including the
# few never detected in sediment, so the habitats are compared on the same set.
in_core <- soil_samples_detected >= soil_core_min
# Mean abundance of each soil-core ASV in soil and in sediment separates ASVs
# that are more abundant in soil from those more abundant in sediment.
core_soil_mean <- rowMeans(counts[in_core, soil_ids])
core_sediment_mean <- rowMeans(sediment_counts[in_core, ])
more_in_sediment <- core_sediment_mean > core_soil_mean
class_summary <- tibble(
  soil_core_asvs = sum(in_core),
  soil_core_asvs_more_abundant_in_soil = sum(!more_in_sediment),
  soil_core_asvs_more_abundant_in_sediment = sum(more_in_sediment),
  soil_core_sediment_reads_from_sediment_favoring = sum(core_sediment_mean[more_in_sediment]) /
    sum(core_sediment_mean),
  soil_core_asvs_in_sediment = sum(in_core & sediment_prevalence > 0),
  sediment_asvs = sum(sediment_prevalence > 0),
  soil_core_sediment_read_share = sum(sediment_counts[in_core, ]) / sum(sediment_counts),
  median_sediment_prevalence = median(sediment_prevalence[in_core & sediment_prevalence > 0]))
write_csv(class_summary, file.path(soil_core_results_dir, "class_summary.csv"))
# Figure 3 splits the soil first detections of the full inventory into soil-core
# and non-core ASVs. Balanced draws keep three soils, too few for the same
# prevalence rule, so they remain a sensitivity of the total soil share.
soil_first <- assignments$first_detection_pool == "Soil"
stopifnot(identical(assignments$asv_id, rownames(counts)), all(soil_first[in_core]),
  sum(in_core) == 1730L)
aquatic_samples <- metadata$habitat != "soil"
split_rows <- list()
for (part in c("Soil core", "Soil non-core")) {
  rows <- in_core
  if (part == "Soil non-core") rows <- soil_first & !in_core
  split_rows[[part]] <- metadata[aquatic_samples, ] %>%
    select(sample_id, habitat, network_group) %>%
    mutate(category = part,
      richness_percent = 100 * colSums(presence[rows, aquatic_samples, drop = FALSE]) /
        colSums(presence[, aquatic_samples, drop = FALSE]),
      read_percent = 100 * colSums(counts[rows, aquatic_samples, drop = FALSE]) / target_depth)
}
soil_split <- bind_rows(split_rows) %>% group_by(habitat, network_group, category) %>%
  summarize(samples = n(), richness_percent = mean(richness_percent),
    read_percent = mean(read_percent), .groups = "drop")
# The two parts add up to the soil category of the full inventory.
soil_totals <- soil_split %>% group_by(habitat, network_group) %>%
  summarize(richness_percent = sum(richness_percent), read_percent = sum(read_percent), .groups = "drop") %>%
  left_join(full_summary %>% filter(category == "Soil"), by = c("habitat", "network_group"),
    suffix = c("_split", "_full"))
stopifnot(max(abs(soil_totals$read_percent_split - soil_totals$read_percent_full)) < 1e-8,
  max(abs(soil_totals$richness_percent_split - soil_totals$richness_percent_full)) < 1e-8)
write_csv(soil_split, file.path(detection_results_dir, "full_inventory_soil_core_split.csv"))

# 5. Compare paired habitats across nested taxonomic resolutions ----------------
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
write_csv(rank_samples, file.path(detection_results_dir, "taxonomic_rank_sample_comparisons.csv"))
write_csv(rank_tests, file.path(detection_results_dir, "taxonomic_rank_tests.csv"))
write_csv(bind_rows(coverage_rows), file.path(detection_results_dir, "taxonomic_rank_coverage.csv"))

# 6. Draw Figure 3; retain earlier panels as supplements -----------------------
plot_data <- full_summary %>% filter(category != "Soil") %>%
  select(habitat, network_group, category, samples, richness_percent, read_percent) %>%
  bind_rows(soil_split) %>%
  mutate(category = factor(category, levels = display_levels),
    habitat = factor(habitat, levels = habitats), network_group = factor(network_group, levels = network_groups))
plots <- list()
for (metric in c("richness_percent", "read_percent")) {
  y_label <- ifelse(metric == "richness_percent", "ASV richness (%)", "ASV reads (%)")
  plots[[metric]] <- ggplot(plot_data, aes(network_group, .data[[metric]], fill = category)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    geom_text(data = distinct(plot_data, habitat, network_group, samples),
      aes(network_group, 106, label = paste0("n = ", samples)), inherit.aes = FALSE,
      size = figure_label_size) +
    facet_wrap(~habitat, nrow = 1, labeller = as_labeller(habitat_labels)) +
    scale_fill_manual(values = display_colors, labels = display_labels,
      drop = FALSE, name = "First detection pool",
      guide = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 114), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = c("Headwater\n(order 1-2)", "Intermediate\n(order 3-4)", "Mainstem\n(order 5)")) +
    labs(x = NULL, y = y_label) + fixed_panels(n_col = 2)
  save_figure(plots[[metric]], file.path(detection_figures_dir, paste0("2016_first_detection_", metric, ".pdf")))
}
# Reuse script 02's formal paired soil-sharing test and verify it from counts.
sharing <- read_csv(input_paths[4], col_types = cols(site_code = col_character()))
soil_test <- read_csv(input_paths[5], show_col_types = FALSE)
stopifnot(!anyDuplicated(sharing$sample_id), soil_test$pairs == 18,
  all(c(pair_audit$sediment_id, pair_audit$porewater_id) %in% sharing$sample_id))
soil_detected <- rowSums(presence[, metadata$habitat == "soil", drop = FALSE]) > 0
expected_sharing <- 100 * colSums(counts[soil_detected, sharing$sample_id, drop = FALSE]) / target_depth
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
  geom_point(aes(color = habitat), size = figure_point_size, alpha = point_alpha) +
  scale_color_manual(values = c(sediment = "#009E73", hyporheic = "#0072B2")) +
  scale_x_discrete(labels = habitat_labels) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Reads shared with\nsoil inventory (%)") +
  plot_theme + theme(legend.position = "none") + fixed_panels()
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
    labs(x = "Taxonomic resolution", y = "Sediment-porewater\nHellinger distance") +
    plot_theme + theme(legend.key.width = grid::unit(1.2, "lines"))
  if (view == "taxonomic_resolution") {
    plots[[view]] <- plots[[view]] + theme(legend.position = "none") +
      fixed_panels(figure_large_panel_width, figure_large_panel_height)
  } else {
    plots[[view]] <- plots[[view]] + facet_wrap(~scope, nrow = 1) +
      fixed_panels(figure_large_panel_width, figure_large_panel_height, n_col = 2)
  }
}
for (panel in c("paired_soil_sharing", "taxonomic_resolution", "taxonomic_resolution_sensitivity")) {
  save_figure(plots[[panel]], file.path(detection_figures_dir, paste0("2016_", panel, ".pdf")))
}
# Figure 3 shows read shares; richness shares remain a supplementary panel.
inventory_panel <- plots$read_percent +
  guides(fill = guide_legend(ncol = 1, title.position = "top")) +
  theme(legend.position = "right", legend.justification = "left")
save_figure(inventory_panel, file.path(manuscript_dir, "fig3_first_detection.pdf"))
if (interactive()) {
  print(inventory_panel)
}

# 7. Place every sample on the stream network and assign its branch ---------------
# The graph follows 03_molecular_gradients.R: consecutive centerline points within a segment
# are joined by their change in outlet distance, and segments join at shared
# confluence points.
geometry <- read_csv(input_paths[8], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
geometry$node <- seq_len(nrow(geometry))
same_segment <- geometry$segment[-1] == geometry$segment[-nrow(geometry)]
channel_edges <- tibble(from = geometry$node[-nrow(geometry)][same_segment],
  to = geometry$node[-1][same_segment])
channel_edges$weight <- abs(geometry$outlet_distance[channel_edges$from] -
  geometry$outlet_distance[channel_edges$to])
shared_points <- split(geometry$node, paste(geometry$x, geometry$y))
shared_points <- shared_points[lengths(shared_points) > 1]
confluence_rows <- list()
for (i in seq_along(shared_points)) {
  nodes <- shared_points[[i]]
  confluence_rows[[i]] <- tibble(from = nodes[1], to = nodes[-1], weight = 0)
}
stream_graph <- graph_from_data_frame(bind_rows(channel_edges, bind_rows(confluence_rows)),
  directed = FALSE, vertices = tibble(name = geometry$node))
# Each sample is placed at the nearest centerline point of its own segment.
receiving <- metadata %>% filter(habitat != "soil")
receiving$centerline_node <- NA_integer_
receiving$centerline_match_error_m <- NA_real_
for (i in seq_len(nrow(receiving))) {
  candidates <- which(geometry$segment == receiving$site_segment_number[i])
  offset <- sqrt((geometry$x[candidates] - receiving$site_utm_x_m[i])^2 +
    (geometry$y[candidates] - receiving$site_utm_y_m[i])^2)
  receiving$centerline_node[i] <- candidates[which.min(offset)]
  receiving$centerline_match_error_m[i] <- min(offset)
}
stopifnot(all(receiving$centerline_match_error_m < 1))
receiving$outlet_distance_m <- geometry$outlet_distance[receiving$centerline_node]
# A sample lies on a branch when water from it passes a point just above the
# fork on that branch: its channel distance to that point equals the difference
# of their outlet distances. The anchor is the third-lowest point of the branch's
# lowest segment, so it lies above the shared confluence point.
upstream_of <- function(segment_number) {
  points <- which(geometry$segment == segment_number)
  anchor <- points[order(geometry$outlet_distance[points])[3]]
  channel <- as.numeric(distances(stream_graph, v = receiving$centerline_node, to = anchor,
    weights = E(stream_graph)$weight))
  difference <- receiving$outlet_distance_m - geometry$outlet_distance[anchor]
  difference >= 0 & abs(channel - difference) < distance_tolerance_m
}
on_north <- upstream_of(north_segment)
on_lookout <- upstream_of(lookout_segment)
stopifnot(!any(on_north & on_lookout))
receiving$branch <- case_when(on_north ~ "North branch", on_lookout ~ "Lookout branch",
  TRUE ~ "Below the fork")

# 8. Measure each sample's soil-core share ----------------------------------------
# Shares are fractions of each sample's 10,000 reads.
receiving$soil_core_share <- colSums(counts[in_core, receiving$sample_id]) / target_depth
sample_shares <- receiving %>% select(sample_id, site_code, habitat, branch, site_stream_order,
  site_drainage_area_ha, site_utm_x_m, site_utm_y_m, outlet_distance_m, soil_core_share)
write_csv(sample_shares, file.path(soil_core_results_dir, "sample_shares.csv"))

# 9. Test where the soil-core share is high ----------------------------------------
sediment_sites <- sample_shares %>% filter(habitat == "sediment")
porewater_sites <- sample_shares %>% filter(habitat == "hyporheic")
paired <- pairs %>% transmute(site_code,
  sediment = sample_shares$soil_core_share[match(sediment_id, sample_shares$sample_id)],
  porewater = sample_shares$soil_core_share[match(porewater_id, sample_shares$sample_id)])
stopifnot(!anyNA(paired))
# Paired difference: sediment minus porewater at the same site, tested by
# exchanging the two habitats within sites (random sign flips).
set.seed(soil_core_seed + 1L)
observed_difference <- mean(paired$sediment - paired$porewater)
sign_null <- replicate(n_permutations,
  mean((paired$sediment - paired$porewater) * sample(c(-1, 1), nrow(paired), replace = TRUE)))
# Spatial autocorrelation along the network: Moran's I with inverse watercourse
# distance weights among the 34 sediment sites.
sediment_index <- match(sediment_sites$sample_id, receiving$sample_id)
watercourse_m <- distances(stream_graph, v = receiving$centerline_node[sediment_index],
  to = receiving$centerline_node[sediment_index], weights = E(stream_graph)$weight)
inverse_distance <- 1 / pmax(watercourse_m / 1000, 0.05)
diag(inverse_distance) <- 0
set.seed(soil_core_seed + 2L)
moran <- moran.mc(sediment_sites$soil_core_share, mat2listw(inverse_distance, style = "W"),
  nsim = n_permutations)
spearman_row <- function(test, x, y) {
  fit <- cor.test(x, y, method = "spearman", exact = FALSE)
  tibble(test, n = length(x), estimate = unname(fit$estimate), p_value = fit$p.value)
}
branch_fit <- kruskal.test(soil_core_share ~ branch, data = sediment_sites)
site_tests <- bind_rows(
  tibble(test = "Soil-core share in sediment minus porewater at paired sites", n = nrow(paired),
    estimate = observed_difference,
    p_value = (1 + sum(abs(sign_null) >= abs(observed_difference))) / (n_permutations + 1)),
  spearman_row("Soil-core share in sediment vs porewater at paired sites", paired$sediment, paired$porewater),
  tibble(test = "Moran's I of sediment soil-core share (inverse watercourse distance)",
    n = nrow(sediment_sites), estimate = unname(moran$statistic), p_value = moran$p.value),
  spearman_row("Sediment soil-core share vs drainage area",
    sediment_sites$soil_core_share, sediment_sites$site_drainage_area_ha),
  spearman_row("Porewater soil-core share vs drainage area",
    porewater_sites$soil_core_share, porewater_sites$site_drainage_area_ha),
  tibble(test = "Sediment soil-core share among branches (Kruskal-Wallis chi-squared)",
    n = nrow(sediment_sites), estimate = unname(branch_fit$statistic), p_value = branch_fit$p.value)) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH"))
write_csv(site_tests, file.path(soil_core_results_dir, "site_tests.csv"))

# 10. Rank the soil-core taxa by their abundance in soil -------------------------------
# ASVs are grouped by family; an ASV without a family assignment is grouped
# under its lowest assigned rank as "Unclassified <rank>". SILVA's placeholder
# families are named by their order. Abundances are mean percentages of
# each habitat's reads, summed over the lineage's soil-core ASVs.
placeholder_families <- c("Unknown Family", "Incertae Sedis")
screened_taxonomy <- read_csv(input_paths[9], show_col_types = FALSE) %>% select(-sequence)
stopifnot(all(rownames(counts) %in% screened_taxonomy$asv_id))
core_taxonomy <- screened_taxonomy[match(rownames(counts)[in_core], screened_taxonomy$asv_id), ] %>%
  mutate(lineage = case_when(
    Family %in% placeholder_families ~ Order,
    !is.na(Family) ~ Family,
    !is.na(Order) ~ paste("Unclassified", Order),
    !is.na(Class) ~ paste("Unclassified", Class),
    TRUE ~ paste("Unclassified", Phylum)))
for (current_habitat in names(habitat_labels)) {
  core_taxonomy[[paste0(current_habitat, "_percent")]] <-
    100 * rowMeans(counts[core_taxonomy$asv_id, metadata$habitat == current_habitat]) / target_depth
}
core_lineages <- core_taxonomy %>% group_by(lineage, Phylum) %>%
  summarize(soil_core_asvs = n(), soil_percent = sum(soil_percent),
    sediment_percent = sum(sediment_percent), porewater_percent = sum(hyporheic_percent),
    .groups = "drop") %>%
  arrange(desc(soil_percent)) %>%
  mutate(soil_rank = row_number(), shown = soil_rank <= n_ranked_lineages,
    soil_to_sediment = soil_percent / sediment_percent,
    share_of_soil_core_soil_reads = soil_percent / sum(soil_percent))
stopifnot(!anyDuplicated(core_lineages$lineage))
write_csv(core_lineages, file.path(soil_core_results_dir, "soil_core_lineages.csv"))

# 11. Draw Figure 4 ---------------------------------------------------------------
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% 10 == 1 | row_number() == n()) %>% ungroup()
map_extent_x <- c(map_lines$x, sediment_sites$site_utm_x_m)
map_extent_y <- c(map_lines$y, sediment_sites$site_utm_y_m)
# The legend gives the number of soil-core ASVs. Read shares run to 50%, above
# every observed value.
soil_core_label <- paste0("Soil-core ASVs (n = ", format(sum(in_core), big.mark = ","), ")\nReads (%)")
share_limits <- c(0, 50)
stopifnot(100 * max(sample_shares$soil_core_share) <= share_limits[2])
panels <- list()
# A: soil-core share of sediment reads.
panels$map <- ggplot(map_lines, aes(x, y, group = segment)) +
  geom_path(color = "#C0C0C0", linewidth = 0.35) +
  geom_point(data = sediment_sites %>% arrange(soil_core_share),
    aes(site_utm_x_m, site_utm_y_m, fill = 100 * soil_core_share),
    inherit.aes = FALSE, shape = 21, size = map_point_size, color = "black", stroke = 0.4) +
  scale_fill_viridis_c(option = "D", limits = share_limits, breaks = seq(0, 50, by = 10),
    name = soil_core_label) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = grid::unit(3, "cm"),
    barheight = grid::unit(0.25, "cm"))) +
  coord_equal() + utm_km_scales + plot_theme +
  fixed_equal_panels(map_extent_x, map_extent_y)
map_panel_height <- figure_map_width * diff(range(map_extent_y)) / diff(range(map_extent_x))
# B: soil-core share against drainage area in both habitats, paired sites
# joined. No line is drawn: neither drainage-area slope passes the threshold.
drainage_tests <- site_tests %>% filter(grepl("vs drainage area", test))
stopifnot(all(drainage_tests$p_bh > significance_threshold))
pair_segments <- pairs %>% transmute(
  site_drainage_area_ha = sample_shares$site_drainage_area_ha[match(sediment_id, sample_shares$sample_id)],
  sediment = 100 * paired$sediment, porewater = 100 * paired$porewater)
drainage_panel_height <- 2.0
panels$drainage <- ggplot(sample_shares %>% mutate(habitat = factor(habitat, levels = names(habitat_labels))),
    aes(site_drainage_area_ha, 100 * soil_core_share, color = habitat, shape = habitat)) +
  geom_segment(data = pair_segments, aes(x = site_drainage_area_ha, xend = site_drainage_area_ha,
    y = porewater, yend = sediment), inherit.aes = FALSE, color = "#BBBBBB", linewidth = 0.3) +
  geom_point(size = figure_point_size, alpha = 0.85) +
  scale_x_log10(limits = c(5, 10^4), breaks = 10^(1:4), labels = scales::label_log()) +
  scale_y_continuous(limits = share_limits) +
  scale_color_manual(values = habitat_colors, labels = habitat_labels, name = NULL) +
  scale_shape_manual(values = habitat_shapes, labels = habitat_labels, name = NULL) +
  labs(x = "Drainage area (ha)", y = "Soil-core reads (%)") +
  plot_theme + fixed_panels(figure_map_width, drainage_panel_height)
# C: the soil-core lineages most abundant in soil, with their abundance in
# sediment and porewater. Rows run from soil-dominant lineages at the top to
# lineages more abundant in sediment than in soil at the bottom, ordered by the
# ratio of soil to sediment abundance; the dashed line marks where the two are
# equal. Abundances span orders of magnitude, so the axis is logarithmic; a
# lineage absent from a habitat is not drawn there. The panel spans both rows
# on the left; its height is set to their combined height.
shown_lineages <- core_lineages %>% filter(shown) %>% arrange(soil_to_sediment)
ranked <- shown_lineages %>%
  mutate(lineage = factor(lineage, levels = lineage)) %>%
  pivot_longer(c(soil_percent, sediment_percent, porewater_percent),
    names_to = "habitat", values_to = "percent") %>%
  mutate(habitat = factor(recode(habitat, soil_percent = "soil", sediment_percent = "sediment",
    porewater_percent = "hyporheic"), levels = names(habitat_labels))) %>%
  filter(percent > 0)
equal_abundance_row <- sum(shown_lineages$soil_to_sediment < 1) + 0.5
lineage_panel_height <- map_panel_height + drainage_panel_height + 1.15
# Every other row is shaded so each lineage can be followed across the panel.
# On the log axis, xmin = 0 maps to the left edge.
row_bands <- tibble(row = seq(1, n_ranked_lineages, by = 2))
panels$lineages <- ggplot(ranked, aes(percent, lineage)) +
  geom_rect(data = row_bands, aes(xmin = 0, xmax = Inf, ymin = row - 0.5, ymax = row + 0.5),
    inherit.aes = FALSE, fill = "#EFEFEF") +
  geom_hline(yintercept = equal_abundance_row, linetype = "dashed", color = "#666666", linewidth = 0.3) +
  geom_line(aes(group = lineage), color = "#BBBBBB", linewidth = 0.3) +
  geom_point(aes(color = habitat, shape = habitat), size = figure_point_size) +
  scale_x_log10(labels = scales::label_log()) +
  scale_color_manual(values = habitat_colors, labels = habitat_labels, name = NULL) +
  scale_shape_manual(values = habitat_shapes, labels = habitat_labels, name = NULL) +
  labs(x = "Mean reads from soil-core ASVs (%)", y = NULL) +
  plot_theme + fixed_panels(figure_panel_width, lineage_panel_height)
panel_files <- c(map = "soil_core_map", drainage = "soil_core_drainage", lineages = "soil_core_lineages")
for (panel in names(panel_files)) {
  save_figure(panels[[panel]], file.path(soil_core_figures_dir, paste0("2016_sediment_", panel_files[[panel]], ".pdf")))
}
figure <- wrap_plots(panels[c("map", "drainage", "lineages")], design = "AC\nBC") +
  plot_annotation(tag_levels = "A")
save_figure(figure, file.path(manuscript_dir, "fig4_soil_affinity.pdf"))

# 12. Verify inputs were not changed and record the settings -----------------------
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(detection_results_dir, "input_checksums.csv"))
write_csv(tibble(soil_core_min, n_permutations, base_seed = soil_core_seed, north_segment, lookout_segment,
  n_ranked_lineages, sediment_sites = nrow(sediment_sites), porewater_sites = nrow(porewater_sites),
  paired_sites = nrow(paired)), file.path(soil_core_results_dir, "settings.csv"))
session <- trimws(capture.output(sessionInfo()), which = "right")
writeLines(session, file.path(detection_results_dir, "session_info.txt"))
writeLines(session, file.path(soil_core_results_dir, "session_info.txt"))
print(full_summary %>% filter(category == "Soil"))
print(balanced_summary %>% filter(category == "Soil"))
message("Saved soil-sharing tables, Figure 3, and Figure 4.")
