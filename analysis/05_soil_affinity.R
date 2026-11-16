# Figure 4: where soil-shared taxa sit in the sediment network.
# Sediment ASVs are classed by how many soil samples they occur in, using the
# soil samples alone. The script asks where the soil-shared fraction is high,
# whether it is carried in porewater, and whether its composition is organized
# along the stream network like the fraction never found in soil.
# Run directly in Positron.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(vegan)
library(spdep)
library(ggplot2)
library(patchwork)
source(here("analysis", "figure_style.R"))

# 1. Settings and plotting theme -------------------------------------------------
# An ASV is soil-core if it occurs in at least soil_core_min of the 15 soils and
# soil-absent if it occurs in none. The classes use soil samples only, so the
# sediment data they are tested on play no part in defining them.
soil_core_min <- 5L
soil_absent_max <- 0L
target_depth <- 10000L
# Composition comparisons use ASVs present at 3 to 31 of the 34 sediment sites;
# an ASV at one or two sites, or at nearly all, carries no spatial pattern.
minimum_sites <- 3L
# Prevalence-matched draws: within each sediment prevalence, the same number of
# ASVs is drawn from each class, so the classes cannot differ through
# prevalence. The null reassigns class labels within each prevalence.
n_draws <- 199L
n_permutations <- 9999L
base_seed <- 20261116L
significance_threshold <- 0.05
# The major fork of the network: segment 134 drains the north branch and
# segment 39 the upper Lookout branch; they meet at 7.7 km from the outlet.
north_segment <- 134
lookout_segment <- 39
distance_tolerance_m <- 1
# Definitions for the class contrast. The first row is the primary analysis.
# Simpson dissimilarity removes richness differences and keeps turnover only.
contrast_design <- tribble(
  ~definition, ~core_min, ~absent_max, ~dissimilarity,
  "primary", 5L, 0L, "jaccard",
  "core_3_soils", 3L, 0L, "jaccard",
  "core_8_soils", 8L, 0L, "jaccard",
  "absent_1_soil", 5L, 1L, "jaccard",
  "turnover", 5L, 0L, "simpson")
contrast_design$seed <- base_seed + 100L * seq_len(nrow(contrast_design))
# Number of soil-core families shown in the ranked panel.
n_ranked_families <- 12L

results_dir <- here("results", "diversity_2016", "soil_affinity")
figures_dir <- here("figures", "supplement", "soil_affinity")
manuscript_dir <- here("figures", "manuscript")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Habitat colors and shapes match Figures 1 and 2.
habitat_labels <- c(soil = "Soil", sediment = "Sediment", hyporheic = "Porewater")
habitat_colors <- c(soil = "#D55E00", sediment = "#009E73", hyporheic = "#0072B2")
habitat_shapes <- c(soil = 17, sediment = 16, hyporheic = 15)
# The stable network group is called the sediment-specific group in figures:
# its ASVs are almost never found in soil or porewater.
group_label <- "Sediment-specific group\nASVs detected (% of 38)"
class_labels <- c(soil_core = "Soil-core ASVs", soil_absent = "Soil-absent ASVs")
class_colors <- c(soil_core = "#D55E00", soil_absent = "#009E73")
separation_labels <- c(watercourse_km = "Watercourse\ndistance", straight_line_km = "Straight-line\ndistance")
map_point_size <- 1.5 * figure_point_size
plot_theme <- figure_theme

# 2. Read the reviewed cohort and the stable group --------------------------------
input_paths <- c(
  here("data", "derived", "soil_porewater_sediment_2016", "asv_counts_10k.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sample_metadata.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sediment_porewater_pair_audit.csv"),
  here("results", "diversity_2016", "soil_porewater_sediment", "soil_inventory.csv"),
  here("results", "fticr_2016", "stable_group_distribution", "core_asv_identities.csv"),
  here("results", "fticr_2016", "stable_group_distribution", "core_site_profiles.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "matched_sediment_dispersion", "molecular_jaccard_distance_matrix.csv"),
  here("data", "StreamCenterline.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_metadata_screened.csv"))
input_md5 <- unname(tools::md5sum(input_paths))
stopifnot(!anyNA(input_md5))
count_table <- read_csv(input_paths[1], show_col_types = FALSE)
counts <- as.matrix(count_table[, -1])
rownames(counts) <- count_table$asv_id
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()), show_col_types = FALSE)
pairs <- read_csv(input_paths[3], col_types = cols(site_code = col_character()), show_col_types = FALSE) %>%
  filter(paired_manuscript)
stopifnot(identical(colnames(counts), metadata$sample_id), all(colSums(counts) == target_depth),
  sum(metadata$habitat == "soil") == 15, sum(metadata$habitat == "sediment") == 34,
  sum(metadata$habitat == "hyporheic") == 27, nrow(pairs) == 18)
soil_ids <- metadata$sample_id[metadata$habitat == "soil"]
sediment <- metadata %>% filter(habitat == "sediment")
receiving <- metadata %>% filter(habitat != "soil")
group_asvs <- read_csv(input_paths[5], show_col_types = FALSE)$feature_id
stopifnot(length(group_asvs) == 38, all(group_asvs %in% rownames(counts)))

# 3. Class every ASV by its occurrence in the soil samples -------------------------
soil_samples_detected <- rowSums(counts[, soil_ids] > 0)
# Script 02 records the same soil occurrence; the two must agree.
soil_inventory <- read_csv(input_paths[4], show_col_types = FALSE)
stopifnot(identical(unname(soil_samples_detected[soil_inventory$asv_id]), soil_inventory$soil_samples_detected))
soil_class <- case_when(
  soil_samples_detected >= soil_core_min ~ "soil_core",
  soil_samples_detected <= soil_absent_max ~ "soil_absent",
  TRUE ~ "intermediate")
names(soil_class) <- rownames(counts)
sediment_counts <- counts[, sediment$sample_id]
sediment_prevalence <- rowSums(sediment_counts > 0)
class_summary <- tibble(asv_id = rownames(counts), soil_class, sediment_prevalence,
  sediment_reads = rowSums(sediment_counts)) %>%
  filter(sediment_prevalence > 0) %>%
  group_by(soil_class) %>%
  summarize(sediment_asvs = n(),
    sediment_read_share = sum(sediment_reads) / (target_depth * nrow(sediment)),
    median_sediment_prevalence = median(sediment_prevalence),
    composition_asvs = sum(sediment_prevalence >= minimum_sites &
      sediment_prevalence <= nrow(sediment) - minimum_sites),
    stable_group_asvs = sum(asv_id %in% group_asvs), .groups = "drop")
write_csv(class_summary, file.path(results_dir, "class_summary.csv"))

# Rank the soil-core taxa by their abundance in sediment. ASVs are grouped by
# family; an ASV without a family assignment is grouped under its lowest
# assigned rank as "Unclassified <rank>". Abundances are mean percentages of
# each habitat's reads, summed over the family's soil-core ASVs.
taxonomy <- read_csv(input_paths[10], show_col_types = FALSE) %>% select(-sequence)
stopifnot(all(rownames(counts) %in% taxonomy$asv_id))
core_taxonomy <- taxonomy[match(rownames(counts)[soil_class == "soil_core"], taxonomy$asv_id), ] %>%
  mutate(family = case_when(
    !is.na(Family) ~ Family,
    !is.na(Order) ~ paste("Unclassified", Order),
    !is.na(Class) ~ paste("Unclassified", Class),
    TRUE ~ paste("Unclassified", Phylum)))
mean_percent <- function(habitat) {
  100 * rowMeans(counts[core_taxonomy$asv_id, metadata$habitat == habitat]) / target_depth
}
core_taxonomy$soil_percent <- mean_percent("soil")
core_taxonomy$sediment_percent <- mean_percent("sediment")
core_taxonomy$porewater_percent <- mean_percent("hyporheic")
core_families <- core_taxonomy %>% group_by(family, Phylum) %>%
  summarize(soil_core_asvs = n(), soil_percent = sum(soil_percent),
    sediment_percent = sum(sediment_percent), porewater_percent = sum(porewater_percent),
    .groups = "drop") %>%
  arrange(desc(sediment_percent))
stopifnot(!anyDuplicated(core_families$family))
write_csv(core_families, file.path(results_dir, "soil_core_families.csv"))

# 4. Place every sample on the stream network and assign its branch ---------------
# The graph follows script 04: consecutive centerline points within a segment
# are joined by their change in outlet distance, and segments join at shared
# confluence points.
geometry <- read_csv(input_paths[9], skip = 1,
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

# 5. Measure each sample's soil-core and soil-absent shares ------------------------
# Shares are fractions of each sample's 10,000 reads. The stable-group columns
# repeat the shares without the group's 38 ASVs, so the comparison with the
# group is not shared membership.
receiving_counts <- counts[, receiving$sample_id]
in_core <- soil_class == "soil_core"
in_absent <- soil_class == "soil_absent"
not_group <- !(rownames(counts) %in% group_asvs)
receiving <- receiving %>% mutate(
  soil_core_share = colSums(receiving_counts[in_core, ]) / target_depth,
  soil_absent_share = colSums(receiving_counts[in_absent, ]) / target_depth,
  soil_core_share_without_group = colSums(receiving_counts[in_core & not_group, ]) / target_depth,
  soil_absent_share_without_group = colSums(receiving_counts[in_absent & not_group, ]) / target_depth,
  stable_group_asv_fraction = colSums(receiving_counts[group_asvs, ] > 0) / length(group_asvs))
# The group fractions must reproduce the saved site profiles, which used the same draw.
profiles <- read_csv(input_paths[6], col_types = cols(site_code = col_character()), show_col_types = FALSE) %>%
  filter(microbial_profile_available)
sediment_rows <- receiving %>% filter(habitat == "sediment")
stopifnot(setequal(profiles$site_code, sediment_rows$site_code),
  max(abs(sediment_rows$stable_group_asv_fraction -
    profiles$core_asv_fraction[match(sediment_rows$site_code, profiles$site_code)])) < 1e-12)
sample_shares <- receiving %>% select(sample_id, site_code, habitat, branch, site_stream_order,
  site_drainage_area_ha, site_utm_x_m, site_utm_y_m, outlet_distance_m, soil_core_share,
  soil_absent_share, soil_core_share_without_group, soil_absent_share_without_group,
  stable_group_asv_fraction)
write_csv(sample_shares, file.path(results_dir, "sample_shares.csv"))

# 6. Test where the soil-core share is high ----------------------------------------
sediment_sites <- sample_shares %>% filter(habitat == "sediment")
porewater_sites <- sample_shares %>% filter(habitat == "hyporheic")
paired <- pairs %>% transmute(site_code,
  sediment = sample_shares$soil_core_share[match(sediment_id, sample_shares$sample_id)],
  porewater = sample_shares$soil_core_share[match(porewater_id, sample_shares$sample_id)])
stopifnot(!anyNA(paired))
# Paired difference: sediment minus porewater at the same site, tested by
# exchanging the two habitats within sites (random sign flips).
set.seed(base_seed + 1L)
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
set.seed(base_seed + 2L)
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
    n = nrow(sediment_sites), estimate = unname(branch_fit$statistic), p_value = branch_fit$p.value),
  spearman_row("Stable-group ASV fraction vs sediment soil-core share (group ASVs removed)",
    sediment_sites$stable_group_asv_fraction, sediment_sites$soil_core_share_without_group),
  spearman_row("Stable-group ASV fraction vs sediment soil-absent share (group ASVs removed)",
    sediment_sites$stable_group_asv_fraction, sediment_sites$soil_absent_share_without_group)) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH"))
write_csv(site_tests, file.path(results_dir, "site_tests.csv"))
stable_group_by_branch <- sediment_sites %>% group_by(branch) %>%
  summarize(sites = n(), sites_with_half_of_group = sum(stable_group_asv_fraction >= 0.5),
    median_group_fraction = median(stable_group_asv_fraction),
    median_soil_core_share = median(soil_core_share), .groups = "drop")
write_csv(stable_group_by_branch, file.path(results_dir, "branch_summary.csv"))

# 7. Compare how each class's composition is organized along the network -----------
# Separations among the 34 sediment sites. Molecular separations use the 33
# sites with FT-ICR profiles and are recorded, not plotted.
sediment_presence <- t(sediment_counts > 0) * 1L
rownames(sediment_presence) <- sediment$site_code
separations <- list(
  watercourse_km = watercourse_m / 1000,
  straight_line_km = as.matrix(dist(cbind(sediment_sites$site_utm_x_m, sediment_sites$site_utm_y_m))) / 1000)
molecular_jaccard <- read_csv(input_paths[8], show_col_types = FALSE)
molecular_sites <- molecular_jaccard[[1]]
molecular_jaccard <- as.matrix(molecular_jaccard[, -1])
dimnames(molecular_jaccard) <- list(molecular_sites, molecular_sites)
molecular_presence <- read_csv(input_paths[7], col_types = cols(site_code = col_character()), show_col_types = FALSE)
molecular_presence <- as.matrix(molecular_presence[match(molecular_sites, molecular_presence$site_code), -1])
rownames(molecular_presence) <- molecular_sites
molecular_simpson <- as.matrix(betadiver(molecular_presence, "sim"))
stopifnot(all(molecular_sites %in% sediment$site_code), length(molecular_sites) == 33)
sediment_lower <- lower.tri(separations$watercourse_km)
molecular_lower <- lower.tri(molecular_jaccard)
# Correlation of a class's dissimilarity with each separation, over site pairs.
composition_dissimilarity <- function(asvs, dissimilarity) {
  if (dissimilarity == "jaccard") {
    as.matrix(vegdist(sediment_presence[, asvs], method = "jaccard", binary = TRUE))
  } else {
    as.matrix(betadiver(sediment_presence[, asvs], "sim"))
  }
}
separation_correlations <- function(asvs, dissimilarity) {
  community <- composition_dissimilarity(asvs, dissimilarity)
  molecular <- if (dissimilarity == "jaccard") molecular_jaccard else molecular_simpson
  c(watercourse_km = cor(community[sediment_lower], separations$watercourse_km[sediment_lower]),
    straight_line_km = cor(community[sediment_lower], separations$straight_line_km[sediment_lower]),
    molecular = cor(community[molecular_sites, molecular_sites][molecular_lower], molecular[molecular_lower]))
}
eligible <- sediment_prevalence >= minimum_sites & sediment_prevalence <= nrow(sediment) - minimum_sites
draw_rows <- list()
for (design_row in seq_len(nrow(contrast_design))) {
  design <- contrast_design[design_row, ]
  core_pool <- rownames(counts)[eligible & soil_samples_detected >= design$core_min]
  absent_pool <- rownames(counts)[eligible & soil_samples_detected <= design$absent_max]
  set.seed(design$seed)
  for (draw in seq_len(n_draws)) {
    core_draw <- c()
    absent_draw <- c()
    for (prevalence in intersect(sediment_prevalence[core_pool], sediment_prevalence[absent_pool])) {
      core_here <- core_pool[sediment_prevalence[core_pool] == prevalence]
      absent_here <- absent_pool[sediment_prevalence[absent_pool] == prevalence]
      k <- min(length(core_here), length(absent_here))
      core_draw <- c(core_draw, core_here[sample.int(length(core_here), k)])
      absent_draw <- c(absent_draw, absent_here[sample.int(length(absent_here), k)])
    }
    # Null: pool the two matched sets and reassign labels within each prevalence.
    pooled <- c(core_draw, absent_draw)
    labels <- rep(c("a", "b"), each = length(core_draw))
    for (prevalence in unique(sediment_prevalence[pooled])) {
      positions <- which(sediment_prevalence[pooled] == prevalence)
      labels[positions] <- labels[positions][sample.int(length(positions))]
    }
    draw_rows[[length(draw_rows) + 1L]] <- tibble(definition = design$definition,
      dissimilarity = design$dissimilarity, draw, asvs_per_class = length(core_draw),
      separation = c("watercourse_km", "straight_line_km", "molecular"),
      soil_core = separation_correlations(core_draw, design$dissimilarity),
      soil_absent = separation_correlations(absent_draw, design$dissimilarity),
      null_a = separation_correlations(pooled[labels == "a"], design$dissimilarity),
      null_b = separation_correlations(pooled[labels == "b"], design$dissimilarity))
  }
}
contrast_draws <- bind_rows(draw_rows)
write_csv(contrast_draws, file.path(results_dir, "class_contrast_draws.csv"))
# The difference test compares the median observed difference with the null
# differences from every draw.
contrast_summary <- contrast_draws %>%
  group_by(definition, dissimilarity, separation) %>%
  summarize(draws = n(), asvs_per_class = median(asvs_per_class),
    soil_core_median = median(soil_core), soil_core_q025 = quantile(soil_core, 0.025),
    soil_core_q975 = quantile(soil_core, 0.975),
    soil_absent_median = median(soil_absent), soil_absent_q025 = quantile(soil_absent, 0.025),
    soil_absent_q975 = quantile(soil_absent, 0.975),
    random_set_q025 = quantile(c(null_a, null_b), 0.025), random_set_q975 = quantile(c(null_a, null_b), 0.975),
    difference_median = median(soil_absent - soil_core),
    draws_absent_greater = mean(soil_absent > soil_core),
    p_difference = (1 + sum(abs(null_b - null_a) >= abs(difference_median))) / (n() + 1),
    .groups = "drop") %>%
  group_by(definition) %>% mutate(p_bh = p.adjust(p_difference, method = "BH")) %>% ungroup()
write_csv(contrast_summary, file.path(results_dir, "class_contrast_summary.csv"))
# Within-class Mantel tests on every eligible ASV of each primary class. Two
# separations by two dissimilarities by two classes form one BH family; the
# molecular tests form a second.
mantel_rows <- list()
for (class_name in names(class_labels)) {
  asvs <- rownames(counts)[eligible & soil_class == class_name]
  for (dissimilarity in c("jaccard", "simpson")) {
    community <- composition_dissimilarity(asvs, dissimilarity)
    molecular <- if (dissimilarity == "jaccard") molecular_jaccard else molecular_simpson
    targets <- list(
      watercourse_km = list(community, separations$watercourse_km),
      straight_line_km = list(community, separations$straight_line_km),
      molecular = list(community[molecular_sites, molecular_sites], molecular))
    for (separation in names(targets)) {
      set.seed(base_seed + 3L)
      fit <- mantel(as.dist(targets[[separation]][[1]]), as.dist(targets[[separation]][[2]]),
        permutations = n_permutations)
      mantel_rows[[length(mantel_rows) + 1L]] <- tibble(soil_class = class_name, asvs = length(asvs),
        dissimilarity, separation, sites = nrow(targets[[separation]][[1]]),
        mantel_r = fit$statistic, p_value = fit$signif)
    }
  }
}
mantel_tests <- bind_rows(mantel_rows) %>%
  mutate(family = if_else(separation == "molecular", "molecular", "spatial")) %>%
  group_by(family) %>% mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
write_csv(mantel_tests, file.path(results_dir, "class_mantel_tests.csv"))

# 8. Draw the figure ---------------------------------------------------------------
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% 10 == 1 | row_number() == n()) %>% ungroup()
map_extent_x <- c(map_lines$x, sediment_sites$site_utm_x_m)
map_extent_y <- c(map_lines$y, sediment_sites$site_utm_y_m)
panels <- list()
# A: soil-core share of sediment reads.
panels$soil_map <- ggplot(map_lines, aes(x, y, group = segment)) +
  geom_path(color = "#C0C0C0", linewidth = 0.35) +
  geom_point(data = sediment_sites %>% arrange(soil_core_share),
    aes(site_utm_x_m, site_utm_y_m, fill = 100 * soil_core_share),
    inherit.aes = FALSE, shape = 21, size = map_point_size, color = "black", stroke = 0.4) +
  scale_fill_viridis_c(option = "D", name = "Soil-core reads (%)") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = grid::unit(3, "cm"),
    barheight = grid::unit(0.25, "cm"))) +
  coord_equal() + utm_km_scales + plot_theme +
  fixed_equal_panels(map_extent_x, map_extent_y)
# B: the sediment-specific group (the stable network group). As in panel A,
# yellow marks the high end.
panels$group_map <- ggplot(map_lines, aes(x, y, group = segment)) +
  geom_path(color = "#C0C0C0", linewidth = 0.35) +
  geom_point(data = sediment_sites %>% arrange(stable_group_asv_fraction),
    aes(site_utm_x_m, site_utm_y_m, fill = 100 * stable_group_asv_fraction),
    inherit.aes = FALSE, shape = 21, size = map_point_size, color = "black", stroke = 0.4) +
  scale_fill_viridis_c(option = "C", limits = c(0, 100), name = group_label) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = grid::unit(3, "cm"),
    barheight = grid::unit(0.25, "cm"))) +
  coord_equal() + utm_km_scales + plot_theme +
  fixed_equal_panels(map_extent_x, map_extent_y)
# C: soil-core share against drainage area in both habitats, paired sites joined.
# No line is drawn: neither drainage-area slope passes the threshold.
pair_segments <- pairs %>% transmute(
  site_drainage_area_ha = sample_shares$site_drainage_area_ha[match(sediment_id, sample_shares$sample_id)],
  sediment = 100 * paired$sediment, porewater = 100 * paired$porewater)
drainage_tests <- site_tests %>% filter(grepl("vs drainage area", test))
stopifnot(all(drainage_tests$p_bh > significance_threshold))
panels$drainage <- ggplot(sample_shares %>% mutate(habitat = factor(habitat, levels = names(habitat_labels))),
    aes(site_drainage_area_ha, 100 * soil_core_share,
    color = habitat, shape = habitat)) +
  geom_segment(data = pair_segments, aes(x = site_drainage_area_ha, xend = site_drainage_area_ha,
    y = porewater, yend = sediment), inherit.aes = FALSE, color = "#BBBBBB", linewidth = 0.3) +
  geom_point(size = figure_point_size, alpha = 0.85) +
  scale_x_log10(limits = c(5, 10^4), breaks = 10^(1:4), labels = scales::label_log()) +
  scale_color_manual(values = habitat_colors, labels = habitat_labels, name = NULL) +
  scale_shape_manual(values = habitat_shapes, labels = habitat_labels, name = NULL) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Drainage area (ha)", y = "Soil-core reads (%)") +
  plot_theme + fixed_panels(figure_large_panel_width, figure_large_panel_height)
# D: the soil-core families most abundant in sediment, with their abundance in
# soil and porewater. Abundances span three orders of magnitude, so the axis is
# logarithmic; a family absent from a habitat is not drawn there. Family names
# sit on the right so they do not push the map above away from panel A.
ranked <- core_families %>% slice_head(n = n_ranked_families) %>%
  mutate(family = factor(family, levels = rev(family))) %>%
  pivot_longer(c(soil_percent, sediment_percent, porewater_percent),
    names_to = "habitat", values_to = "percent") %>%
  mutate(habitat = factor(recode(habitat, soil_percent = "soil", sediment_percent = "sediment",
    porewater_percent = "hyporheic"), levels = names(habitat_labels))) %>%
  filter(percent > 0)
panels$families <- ggplot(ranked, aes(percent, family)) +
  geom_line(aes(group = family), color = "#BBBBBB", linewidth = 0.3) +
  geom_point(aes(color = habitat, shape = habitat), size = figure_point_size) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_discrete(position = "right") +
  scale_color_manual(values = habitat_colors, labels = habitat_labels, name = NULL) +
  scale_shape_manual(values = habitat_shapes, labels = habitat_labels, name = NULL) +
  labs(x = "Mean reads from soil-core ASVs (%)", y = NULL) +
  plot_theme + fixed_panels(figure_large_panel_width, figure_large_panel_height)
# Supplement: how strongly each class's composition tracks distance, primary
# definition. Points are draw medians with 95% draw ranges; the gray band is the
# 95% range for random prevalence-matched sets drawn from both classes.
primary <- contrast_summary %>%
  filter(definition == "primary", separation %in% names(separation_labels)) %>%
  mutate(separation = factor(separation, levels = names(separation_labels)))
primary_long <- bind_rows(
  primary %>% transmute(separation, soil_class = "soil_core", median = soil_core_median,
    low = soil_core_q025, high = soil_core_q975),
  primary %>% transmute(separation, soil_class = "soil_absent", median = soil_absent_median,
    low = soil_absent_q025, high = soil_absent_q975)) %>%
  mutate(soil_class = factor(soil_class, levels = names(class_labels)))
panels$contrast <- ggplot(primary_long, aes(as.numeric(separation), median, color = soil_class)) +
  geom_hline(yintercept = 0, color = "#999999", linewidth = 0.3) +
  geom_rect(data = primary, aes(xmin = as.numeric(separation) - 0.3, xmax = as.numeric(separation) + 0.3,
    ymin = random_set_q025, ymax = random_set_q975), inherit.aes = FALSE, fill = "#E5E5E5") +
  geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.4),
    size = 0.25, linewidth = figure_line_width) +
  scale_x_continuous(breaks = seq_along(separation_labels), labels = separation_labels,
    limits = c(0.5, length(separation_labels) + 0.5)) +
  scale_color_manual(values = class_colors, labels = class_labels, name = NULL) +
  labs(x = NULL, y = "Correlation of Jaccard\ndissimilarity with distance") +
  plot_theme + fixed_panels(figure_large_panel_width, figure_large_panel_height)
panel_files <- c(soil_map = "soil_core_map", group_map = "specific_group_map",
  drainage = "soil_core_drainage", families = "soil_core_families",
  contrast = "class_distance_contrast")
for (panel in names(panel_files)) {
  save_figure(panels[[panel]], file.path(figures_dir, paste0("2016_sediment_", panel_files[[panel]], ".pdf")))
}
figure <- wrap_plots(panels[c("soil_map", "group_map", "drainage", "families")], design = "AB\nCD") +
  plot_annotation(tag_levels = "A")
save_figure(figure, file.path(manuscript_dir, "2016_sediment_soil_affinity.pdf"))

# 9. Verify inputs were not changed and record the settings ------------------------
stopifnot(identical(input_md5, unname(tools::md5sum(input_paths))))
write_csv(tibble(soil_core_min, soil_absent_max, minimum_sites, n_draws, n_permutations, base_seed,
  north_segment, lookout_segment, sediment_sites = nrow(sediment_sites),
  porewater_sites = nrow(porewater_sites), paired_sites = nrow(paired)),
  file.path(results_dir, "settings.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))
message("Saved soil-affinity tables and Figure 4.")
