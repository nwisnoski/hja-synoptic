# Extract exploratory stable-link groups from completed segment-removal refits.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)

# 1. Settings ------------------------------------------------------------------
test_run <- FALSE
minimum_support <- 0.8
minimum_nodes_per_type <- 3L
figures_dir <- here("figures")
figure_prefix <- "2016_sediment_bipartite"
display_group_rank <- 1L
node_colors <- c(Acidobacteriota = "#009E73", Pseudomonadota = "#0072B2",
  Gemmatimonadota = "#E69F00", "Other phyla" = "#999999", "Molecular feature" = "#333333")
plot_theme <- theme_void(base_size = 11) + theme(legend.position = "bottom")
results_dir <- here("results", "fticr_2016", "bipartite_networks")
if (test_run) {
  results_dir <- here("results", "fticr_2016", "bipartite_network_pilot")
  figure_prefix <- "pilot_2016_sediment_bipartite"
}
input_paths <- file.path(results_dir, c("candidate_edges.csv", "feature_pattern_map.csv",
  "asv_module_taxonomy.csv", "molecular_module_properties.csv",
  "leave_segment_out_memberships.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
stopifnot(!anyNA(checks$initial_md5))

# 2. Retain stable positive links with the same direction after adjustment --------
edges <- read_csv(input_paths[1], show_col_types = FALSE)
mapping <- read_csv(input_paths[2], show_col_types = FALSE)
asvs <- read_csv(input_paths[3], show_col_types = FALSE)
molecules <- read_csv(input_paths[4], show_col_types = FALSE)
refits <- read_csv(input_paths[5], show_col_types = FALSE)
stable_edges <- edges %>% filter(sign == "Positive",
  segment_threshold_support >= minimum_support,
  segment_coassignment >= minimum_support, partial_phi > 0)
stopifnot(nrow(stable_edges) > 0)
graph <- graph_from_data_frame(stable_edges %>% select(asv_pattern, molecular_pattern), directed = FALSE)
membership <- components(graph)$membership
nodes <- tibble(pattern_id = names(membership), core = unname(membership)) %>%
  left_join(mapping %>% distinct(node_type, pattern_id), by = "pattern_id")
stopifnot(!anyNA(nodes$node_type), !anyDuplicated(nodes$pattern_id))
stable_edges <- stable_edges %>% left_join(nodes %>% select(asv_pattern = pattern_id, core),
  by = "asv_pattern")
write_csv(stable_edges, file.path(results_dir, "stable_core_edges.csv"))
write_csv(nodes, file.path(results_dir, "stable_core_patterns.csv"))
feature_membership <- mapping %>% left_join(nodes %>% select(pattern_id, core), by = "pattern_id")
write_csv(feature_membership, file.path(results_dir, "stable_core_feature_membership.csv"))

# Connected stable links can form a chain. Check all cross-type pairs in each
# resulting component, rather than assuming every pair has 80% co-membership.
core_rows <- list()
for (core_id in sort(unique(nodes$core))) {
  core_asvs <- nodes$pattern_id[nodes$core == core_id & nodes$node_type == "ASV"]
  core_molecules <- nodes$pattern_id[nodes$core == core_id & nodes$node_type == "Molecular"]
  pair_support <- matrix(0, length(core_asvs), length(core_molecules))
  segments <- sort(unique(refits$omitted_segment))
  for (segment in segments) {
    current <- refits %>% filter(omitted_segment == segment)
    row_labels <- current$module[match(core_asvs, current$pattern_id)]
    column_labels <- current$module[match(core_molecules, current$pattern_id)]
    same <- outer(row_labels, column_labels, "==")
    same[is.na(same)] <- FALSE
    pair_support <- pair_support + same / length(segments)
  }
  selected <- stable_edges %>% filter(core == core_id)
  original_asvs <- asvs %>% filter(pattern_id %in% core_asvs)
  original_molecules <- molecules %>% filter(pattern_id %in% core_molecules)
  phyla <- original_asvs %>% count(Phylum, sort = TRUE)
  core_rows[[length(core_rows) + 1L]] <- tibble(core = core_id,
    asv_patterns = length(core_asvs), molecular_patterns = length(core_molecules),
    original_asvs = nrow(original_asvs), original_molecular_features = nrow(original_molecules),
    stable_edges = nrow(selected), median_phi = median(selected$phi),
    mean_partial_phi = mean(selected$partial_phi),
    mean_edge_threshold_support = mean(selected$segment_threshold_support),
    mean_all_cross_pair_coassignment = mean(pair_support),
    minimum_cross_pair_coassignment = min(pair_support),
    fraction_cross_pairs_support_80 = mean(pair_support >= minimum_support - 1e-10),
    dominant_phylum = phyla$Phylum[1], dominant_phylum_asvs = phyla$n[1],
    mean_h_c = mean(original_molecules$hydrogen_carbon_ratio),
    mean_o_c = mean(original_molecules$oxygen_carbon_ratio),
    fraction_n = mean(original_molecules$nitrogen_count > 0),
    fraction_s = mean(original_molecules$sulfur_count > 0),
    fraction_p = mean(original_molecules$phosphorus_count > 0),
    substantive_group = length(core_asvs) >= minimum_nodes_per_type &
      length(core_molecules) >= minimum_nodes_per_type)
}
summary <- bind_rows(core_rows) %>% arrange(desc(stable_edges))
write_csv(summary, file.path(results_dir, "stable_core_summary.csv"))
write_csv(asvs %>% left_join(nodes %>% select(pattern_id, core), by = "pattern_id") %>%
  filter(!is.na(core)), file.path(results_dir, "stable_core_asv_taxonomy.csv"))
write_csv(molecules %>% left_join(nodes %>% select(pattern_id, core), by = "pattern_id") %>%
  filter(!is.na(core)), file.path(results_dir, "stable_core_molecular_properties.csv"))

# 3. Record the exploratory rule and verify that inference inputs are unchanged --
write_csv(tibble(minimum_support, minimum_nodes_per_type,
  positive_partial_direction_required = TRUE, stable_links = nrow(stable_edges),
  connected_groups = nrow(summary), substantive_groups = sum(summary$substantive_group)),
  file.path(results_dir, "stable_core_settings.csv"))
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
write_csv(checks, file.path(results_dir, "stable_core_input_integrity.csv"))
stopifnot(all(checks$unchanged))
print(summary %>% filter(substantive_group), width = Inf)

# 4. Draw the largest substantive stable-link group with traceable node labels ---
display_groups <- summary %>% filter(substantive_group)
if (nrow(display_groups) >= display_group_rank) {
  display_core <- display_groups$core[display_group_rank]
  pattern_labels <- mapping %>% group_by(pattern_id) %>%
    summarize(representative_id = first(representative_id), original_features = n(), .groups = "drop")
  asv_nodes <- nodes %>% filter(core == display_core, node_type == "ASV") %>%
    left_join(pattern_labels, by = "pattern_id") %>%
    left_join(asvs %>% select(feature_id, Phylum), by = c("representative_id" = "feature_id")) %>%
    mutate(color_group = if_else(Phylum %in% names(node_colors), Phylum, "Other phyla")) %>%
    arrange(color_group, representative_id) %>%
    mutate(x = 0, y = row_number(), label = if_else(original_features > 1,
      paste0(representative_id, " (+", original_features - 1L, ")"), representative_id))
  selected_edges <- stable_edges %>% filter(core == display_core) %>%
    left_join(asv_nodes %>% select(asv_pattern = pattern_id, y), by = "asv_pattern")
  molecular_nodes <- nodes %>% filter(core == display_core, node_type == "Molecular") %>%
    left_join(pattern_labels, by = "pattern_id") %>%
    left_join(selected_edges %>% group_by(pattern_id = molecular_pattern) %>%
      summarize(mean_asv_position = mean(y), .groups = "drop"), by = "pattern_id") %>%
    arrange(mean_asv_position, representative_id) %>%
    mutate(x = 1, y = seq(1, nrow(asv_nodes), length.out = n()),
      color_group = "Molecular feature", label = representative_id)
  displayed_nodes <- bind_rows(asv_nodes, molecular_nodes)
  displayed_edges <- selected_edges %>% rename(y_start = y) %>%
    left_join(molecular_nodes %>% select(molecular_pattern = pattern_id, y_end = y),
      by = "molecular_pattern")
  core_plot <- ggplot() +
    geom_segment(data = displayed_edges, aes(x = 0, xend = 1, y = y_start, yend = y_end),
      color = "#0072B2", alpha = 0.25, linewidth = 0.35) +
    geom_point(data = displayed_nodes, aes(x, y, color = color_group, shape = node_type), size = 2.7) +
    geom_text(data = asv_nodes, aes(x = -0.025, y, label = label), hjust = 1, size = 2.6) +
    geom_text(data = molecular_nodes, aes(x = 1.025, y, label = label), hjust = 0, size = 2.9) +
    scale_color_manual(values = node_colors, name = NULL) +
    scale_shape_manual(values = c(ASV = 16, Molecular = 15), guide = "none") +
    coord_cartesian(xlim = c(-0.55, 1.5), ylim = c(0, nrow(asv_nodes) + 1)) + plot_theme
  ggsave(file.path(figures_dir, paste0(figure_prefix, "_stable_group.pdf")), core_plot,
    width = 9.2, height = 9.5, bg = "white")
  write_csv(displayed_nodes, file.path(results_dir, "stable_group_figure_nodes.csv"))
}
