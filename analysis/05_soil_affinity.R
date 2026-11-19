# Figure 4: diffuse soil-shared taxa in sediment across the stream network.
# Sediment ASVs found in at least five of the 15 soils are called soil-core,
# using the soil samples alone. The script asks where their share of sediment
# reads is high, whether porewater carries the same share, and which taxa they
# are.
# Run directly in Positron.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(spdep)
library(ggplot2)
library(patchwork)
source(here("analysis", "figure_style.R"))

# 1. Settings and plotting theme -------------------------------------------------
# An ASV is soil-core if it occurs in at least soil_core_min of the 15 soils. The
# class uses soil samples only, so the sediment data it is tested on play no
# part in defining it.
soil_core_min <- 5L
target_depth <- 10000L
n_permutations <- 9999L
base_seed <- 20261116L
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

results_dir <- here("results", "diversity_2016", "soil_affinity")
figures_dir <- here("figures", "supplement", "soil_affinity")
manuscript_dir <- here("figures", "manuscript")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Habitat colors and shapes match Figures 1 and 3.
habitat_labels <- c(soil = "Soil", sediment = "Sediment", hyporheic = "Porewater")
habitat_colors <- c(soil = "#D55E00", sediment = "#009E73", hyporheic = "#0072B2")
habitat_shapes <- c(soil = 17, sediment = 16, hyporheic = 15)
map_point_size <- 1.5 * figure_point_size
plot_theme <- figure_theme

# 2. Read the reviewed cohort -------------------------------------------------------
input_paths <- c(
  here("data", "derived", "soil_porewater_sediment_2016", "asv_counts_10k.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sample_metadata.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sediment_porewater_pair_audit.csv"),
  here("results", "diversity_2016", "soil_porewater_sediment", "soil_inventory.csv"),
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

# 3. Identify the soil-core ASVs ----------------------------------------------------
soil_samples_detected <- rowSums(counts[, soil_ids] > 0)
# Script 02 records the same soil occurrence; the two must agree.
soil_inventory <- read_csv(input_paths[4], show_col_types = FALSE)
stopifnot(identical(unname(soil_samples_detected[soil_inventory$asv_id]), soil_inventory$soil_samples_detected))
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
write_csv(class_summary, file.path(results_dir, "class_summary.csv"))

# 4. Place every sample on the stream network and assign its branch ---------------
# The graph follows script 03: consecutive centerline points within a segment
# are joined by their change in outlet distance, and segments join at shared
# confluence points.
geometry <- read_csv(input_paths[5], skip = 1,
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

# 5. Measure each sample's soil-core share ----------------------------------------
# Shares are fractions of each sample's 10,000 reads.
receiving$soil_core_share <- colSums(counts[in_core, receiving$sample_id]) / target_depth
sample_shares <- receiving %>% select(sample_id, site_code, habitat, branch, site_stream_order,
  site_drainage_area_ha, site_utm_x_m, site_utm_y_m, outlet_distance_m, soil_core_share)
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
    n = nrow(sediment_sites), estimate = unname(branch_fit$statistic), p_value = branch_fit$p.value)) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH"))
write_csv(site_tests, file.path(results_dir, "site_tests.csv"))

# 7. Rank the soil-core taxa by their abundance in soil -------------------------------
# ASVs are grouped by family; an ASV without a family assignment is grouped
# under its lowest assigned rank as "Unclassified <rank>". SILVA's placeholder
# families are named by their order. Abundances are mean percentages of
# each habitat's reads, summed over the lineage's soil-core ASVs.
placeholder_families <- c("Unknown Family", "Incertae Sedis")
taxonomy <- read_csv(input_paths[6], show_col_types = FALSE) %>% select(-sequence)
stopifnot(all(rownames(counts) %in% taxonomy$asv_id))
core_taxonomy <- taxonomy[match(rownames(counts)[in_core], taxonomy$asv_id), ] %>%
  mutate(lineage = case_when(
    Family %in% placeholder_families ~ Order,
    !is.na(Family) ~ Family,
    !is.na(Order) ~ paste("Unclassified", Order),
    !is.na(Class) ~ paste("Unclassified", Class),
    TRUE ~ paste("Unclassified", Phylum)))
for (habitat in names(habitat_labels)) {
  core_taxonomy[[paste0(habitat, "_percent")]] <-
    100 * rowMeans(counts[core_taxonomy$asv_id, metadata$habitat == habitat]) / target_depth
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
write_csv(core_lineages, file.path(results_dir, "soil_core_lineages.csv"))

# 8. Draw the figure ---------------------------------------------------------------
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
  save_figure(panels[[panel]], file.path(figures_dir, paste0("2016_sediment_", panel_files[[panel]], ".pdf")))
}
figure <- wrap_plots(panels[c("map", "drainage", "lineages")], design = "AC\nBC") +
  plot_annotation(tag_levels = "A")
save_figure(figure, file.path(manuscript_dir, "fig4_soil_affinity.pdf"))

# 9. Verify inputs were not changed and record the settings ------------------------
stopifnot(identical(input_md5, unname(tools::md5sum(input_paths))))
write_csv(tibble(soil_core_min, n_permutations, base_seed, north_segment, lookout_segment,
  n_ranked_lineages, sediment_sites = nrow(sediment_sites), porewater_sites = nrow(porewater_sites),
  paired_sites = nrow(paired)), file.path(results_dir, "settings.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))
message("Saved soil-affinity tables and Figure 4.")
