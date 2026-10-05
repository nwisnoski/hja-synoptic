# Describe the earliest surveyed pool in which each ASV was detected.
# Adapt Ruiz-Gonzalez et al. (2015), doi:10.1111/ele.12499, to parallel
# aquatic compartments. Earliest detection is not a demonstrated source.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
target_depth <- 10000
detection_threshold <- 1
balanced_draws <- 1000
balanced_samples_per_cell <- 3
base_seed <- 20261005
habitats <- c("sediment", "hyporheic", "planktonic")
habitat_labels <- c(sediment = "Sediment", hyporheic = "Hyporheic", planktonic = "Planktonic")
source_pools <- c("Sediment", "Aquatic")
network_groups <- c("Headwater", "Intermediate", "Mainstem")
category_labels <- c("Soil", "Sediment", "Aquatic", "Sediment + Aquatic")
category_colors <- c("Soil" = "#D55E00", "Sediment" = "#009E73",
  "Aquatic" = "#0072B2", "Sediment + Aquatic" = "#CC79A7")
results_dir <- here("results", "diversity_2016", "first_detection_flowpaths")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    strip.background = element_blank(), legend.position = "bottom",
    legend.title = element_text(size = 10), legend.text = element_text(size = 10))
theme_set(plot_theme)

# 2. Align fixed-depth samples and retain every identity in an audit -------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
count_table <- read_csv(input_paths[1], show_col_types = FALSE)
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()), show_col_types = FALSE)
counts <- as.matrix(count_table[, -1])
rownames(counts) <- count_table$asv_id
metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
stopifnot(identical(metadata$sample_id, colnames(counts)), !anyDuplicated(metadata$sample_id),
  !anyDuplicated(rownames(counts)), all(metadata$included_10k),
  all(colSums(counts) == target_depth), all(rowSums(counts) > 0),
  all(counts >= 0 & counts == floor(counts)), detection_threshold == 1)
metadata <- metadata %>% mutate(network_group = case_when(
  habitat == "soil" ~ "Soil", site_stream_order %in% 1:2 ~ "Headwater",
  site_stream_order %in% 3:4 ~ "Intermediate", site_stream_order == 5 ~ "Mainstem"),
  source_pool = case_when(habitat == "soil" ~ "Soil", habitat == "sediment" ~ "Sediment",
    habitat %in% c("hyporheic", "planktonic") ~ "Aquatic"),
  cell = paste(habitat, network_group, sep = "_"))
stopifnot(nrow(metadata) == 100, sum(metadata$habitat == "soil") == 15,
  !anyNA(metadata$network_group), all(metadata$habitat %in% c("soil", habitats)))
write_csv(metadata %>% select(sample_id, site_code, habitat, mapping_status,
  site_stream_order, site_segment_number, network_group, source_pool, cell),
  file.path(results_dir, "sample_network_audit.csv"))
cell_design <- metadata %>% group_by(cell, habitat, source_pool, network_group) %>%
  summarize(samples = n(), segments = n_distinct(site_segment_number, na.rm = TRUE), .groups = "drop")
stopifnot(nrow(cell_design) == 10, min(cell_design$samples) == balanced_samples_per_cell)
write_csv(cell_design, file.path(results_dir, "sampling_cells.csv"))
presence <- counts >= detection_threshold
cell_indices <- split(seq_len(nrow(metadata)), metadata$cell)

# 3. Repeat the detection inventory with all samples and with balanced sampling --
# Soils have precedence. For remaining taxa, use the earliest network group
# containing a detection. Hyporheic and planktonic detections form one water
# pool. Sediment and water have no assumed flow order; ties remain shared.
# Equal-sample draws retain three profiles per original habitat-stage cell,
# including soils. The merged water pool therefore contains six per stage;
# this balances original habitat coverage, not the merged source pools.
# Original receiving habitats remain separate in the plotted panels. Each draw
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
    stage_profiles <- lapply(0:3, function(stage_id) {
      rows <- first_stage == stage_id
      metadata[aquatic_columns, ] %>% select(sample_id, habitat, network_group) %>%
        mutate(first_network_group = c("Soil", network_groups)[stage_id + 1L],
          richness_percent = 100 * colSums(presence[rows, aquatic_columns, drop = FALSE]) / observed_asvs)
    })
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

# 4. Show richness and read representation along the network ---------------------
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
    scale_fill_manual(values = category_colors, drop = FALSE, name = "First detection pool",
      guide = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = c("Headwater\n(order 1-2)", "Intermediate\n(order 3-4)", "Mainstem\n(order 5)")) +
    labs(x = NULL, y = y_label)
  ggsave(here("figures", paste0("2016_first_detection_", metric, ".pdf")),
    plots[[metric]], width = 10.5, height = 4.3, bg = "white")
}
combined <- (plots$richness_percent / plots$read_percent) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(here("figures", "2016_first_detection_flowpaths.pdf"), combined,
  width = 10.5, height = 7.4, bg = "white")
if (interactive()) {
  print(combined)
}

# 5. Retain settings and verify that the community inputs were not changed -------
write_csv(tibble(target_depth = target_depth, detection_threshold = detection_threshold,
  balanced_draws = balanced_draws, balanced_samples_per_cell = balanced_samples_per_cell,
  base_seed = base_seed, aquatic_pool = "Hyporheic + planktonic",
  assignment = "Soil first; earliest network group; retain sediment-water ties",
  summary = "Mean within-sample percentages"), file.path(results_dir, "analysis_settings.csv"))
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "input_checksums.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))
print(full_summary %>% filter(category == "Soil"))
print(balanced_summary %>% filter(category == "Soil"))
