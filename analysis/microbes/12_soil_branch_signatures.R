# Compare soil signatures from two disjoint tributary drainage areas.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(ggrepel)

# 1. Settings and plotting theme ------------------------------------------------
branch_outlet_sites <- c("47", "66")
minimum_source_soils <- 2
minimum_source_reads <- 5
branch_abundance_ratio <- 3
target_depth <- 10000
make_branch_figures <- FALSE
habitats <- c("planktonic", "hyporheic", "sediment")
source_colors <- c(`47` = "#0072B2", `66` = "#D55E00")
group_levels <- c("Branch 47", "Branch 66", "Downstream of both", "Other")
label_seed <- 20261005
results_dir <- here("results", "diversity_2016", "soil_stream_localization")
routing_dir <- here("data", "derived", "soil_drainage_2016")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
theme_set(plot_theme)

# 2. Audit the two drainage polygons and their confidently mapped soil sources ---
input_paths <- c(
  file.path(routing_dir, "aquatic_catchments_30m.geojson"),
  file.path(routing_dir, "soil_site_membership.csv"),
  file.path(routing_dir, "soil_coordinate_audit.csv"),
  file.path(routing_dir, "aquatic_catchment_audit.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  file.path(results_dir, "soil_asv_identities.csv"),
  here("data", "StreamCenterline.csv"))
input_checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
catchments <- st_read(input_paths[1], quiet = TRUE) %>% group_by(site_code) %>% summarize()
branches <- catchments %>% filter(site_code %in% branch_outlet_sites)
branches <- branches[match(branch_outlet_sites, branches$site_code), ]
intersection_area <- sum(as.numeric(st_area(suppressWarnings(st_intersection(branches[1, ], branches[2, ])))))
stopifnot(nrow(branches) == 2, intersection_area == 0)
write_csv(tibble(branch_outlet_site = branches$site_code,
  drainage_area_ha = as.numeric(st_area(branches)) / 10000,
  overlap_area_ha = intersection_area / 10000), file.path(results_dir, "branch_polygon_audit.csv"))
membership <- read_csv(input_paths[2], col_types = cols(aquatic_site_code = col_character()), show_col_types = FALSE)
soils <- read_csv(input_paths[3], show_col_types = FALSE)
source_soils <- membership %>% filter(aquatic_site_code %in% branch_outlet_sites,
  contributing, robust_spatial_included) %>% select(branch = aquatic_site_code,
    sample_id = soil_sample_id, soil_site_code) %>% left_join(soils, by = "sample_id")
stopifnot(!anyDuplicated(source_soils$sample_id), all(table(source_soils$branch) >= minimum_source_soils))
write_csv(source_soils, file.path(results_dir, "branch_source_soils.csv"))

# 3. Classify aquatic samples with the full drainage footprints ------------------
audit <- read_csv(input_paths[4], col_types = cols(site_code = col_character()), show_col_types = FALSE)
site_qc <- audit %>% group_by(site_code) %>% summarize(all_radii_qc_pass = all(routing_qc_pass), .groups = "drop")
outlets <- audit %>% filter(snap_radius_m == 30) %>%
  st_as_sf(coords = c("outlet_x_m", "outlet_y_m"), crs = 26910, remove = FALSE)
source_outlets <- outlets[outlets$site_code %in% branch_outlet_sites, ]
site_rows <- list()
for (i in seq_len(nrow(outlets))) {
  inside <- st_intersects(outlets[i, ], branches)[[1]]
  current_catchment <- catchments[catchments$site_code == outlets$site_code[i], ]
  drains_both <- length(st_intersects(current_catchment, source_outlets)[[1]]) == 2
  group <- "Other"
  if (length(inside) == 1) {
    group <- paste("Branch", branches$site_code[inside])
  } else if (drains_both) {
    group <- "Downstream of both"
  }
  # Keep an explicit outlet-sensitivity flag for branch territory assignment.
  alternatives <- audit %>% filter(site_code == outlets$site_code[i]) %>%
    st_as_sf(coords = c("outlet_x_m", "outlet_y_m"), crs = 26910)
  alternative_inside <- st_intersects(alternatives, branches)
  alternative_groups <- rep("Outside both branches", length(alternative_inside))
  for (j in seq_along(alternative_inside)) {
    if (length(alternative_inside[[j]]) == 1) {
      alternative_groups[j] <- paste("Branch", branches$site_code[alternative_inside[[j]]])
    }
  }
  site_rows[[i]] <- tibble(site_code = outlets$site_code[i], branch_context = group,
    branch_territory_stable = n_distinct(alternative_groups) == 1)
}
site_context <- bind_rows(site_rows) %>% left_join(site_qc, by = "site_code") %>%
  mutate(included = branch_territory_stable & all_radii_qc_pass)
write_csv(site_context, file.path(results_dir, "branch_site_context_audit.csv"))
metadata <- read_csv(input_paths[6], col_types = cols(site_code = col_character()), show_col_types = FALSE) %>%
  filter(included_10k, habitat %in% habitats) %>% left_join(site_context, by = "site_code")
write_csv(metadata %>% count(habitat, branch_context, included, name = "samples"),
  file.path(results_dir, "branch_sample_counts.csv"))
count_table <- read_csv(input_paths[5], show_col_types = FALSE)
counts <- as.matrix(count_table[, -1])
rownames(counts) <- count_table$asv_id
stopifnot(all(colSums(counts) == target_depth))
taxonomy <- read_csv(input_paths[7], show_col_types = FALSE)

# 4. Define ASV signatures using soils alone, then describe their aquatic reads ---
# "Concentrated" compares the two sampled source pools; it does not establish
# exclusivity to either watershed or use aquatic responses to select taxa.
source_counts <- list()
for (branch_name in branch_outlet_sites) {
  source_ids <- source_soils$sample_id[source_soils$branch == branch_name]
  source_counts[[branch_name]] <- counts[, source_ids, drop = FALSE]
}
feature_rows <- list()
profile_rows <- list()
for (branch_name in branch_outlet_sites) {
  other_branch <- setdiff(branch_outlet_sites, branch_name)
  current_counts <- source_counts[[branch_name]]
  other_counts <- source_counts[[other_branch]]
  current <- tibble(asv_id = rownames(counts), branch = branch_name,
    source_soils_detected = rowSums(current_counts > 0), source_reads = rowSums(current_counts),
    source_mean_read_percent = 100 * rowMeans(current_counts) / target_depth,
    other_branch_soils_detected = rowSums(other_counts > 0),
    other_branch_mean_read_percent = 100 * rowMeans(other_counts) / target_depth) %>%
    mutate(branch_mean_ratio = source_mean_read_percent / other_branch_mean_read_percent,
      branch_detected = source_soils_detected > 0,
      branch_concentrated = source_soils_detected >= minimum_source_soils &
        source_reads >= minimum_source_reads & branch_mean_ratio >= branch_abundance_ratio) %>%
    left_join(taxonomy %>% select(asv_id, Phylum, Genus, soil_enriched, soil_mean_read_percent), by = "asv_id")
  feature_rows[[length(feature_rows) + 1]] <- current
  for (pool_name in c("branch_detected", "branch_concentrated")) {
    selected <- current$asv_id[current[[pool_name]]]
    profile_rows[[length(profile_rows) + 1]] <- metadata %>% select(sample_id, site_code, habitat,
      branch_context, included, site_stream_order, site_segment_number, site_drainage_area_ha,
      site_utm_x_m, site_utm_y_m) %>%
      mutate(source_branch = branch_name, pool = pool_name, source_pool_asvs = length(selected),
        source_asvs_detected = colSums(counts[selected, sample_id, drop = FALSE] > 0),
        source_detection_fraction = source_asvs_detected / source_pool_asvs,
        source_read_percent = 100 * colSums(counts[selected, sample_id, drop = FALSE]) / target_depth)
  }
}
feature_table <- bind_rows(feature_rows)
profiles <- bind_rows(profile_rows)
write_csv(feature_table, file.path(results_dir, "branch_soil_asv_signatures.csv"))
write_csv(profiles, file.path(results_dir, "branch_signature_site_profiles.csv"))
summary <- profiles %>% filter(included) %>% group_by(habitat, source_branch, pool, branch_context) %>%
  summarize(samples = n(), segments = n_distinct(site_segment_number),
    source_pool_asvs = first(source_pool_asvs), median_read_percent = median(source_read_percent),
    median_detected_asvs = median(source_asvs_detected), .groups = "drop")
write_csv(summary, file.path(results_dir, "branch_signature_summary.csv"))

# Individual ASVs can depart from the aggregate response. Keep nondetections,
# and do not treat the small branch samples as independent calibrated tests.
asv_rows <- list()
for (branch_name in branch_outlet_sites) {
  selected <- feature_table$asv_id[feature_table$branch == branch_name & feature_table$branch_concentrated]
  for (habitat_name in habitats) {
    for (context in group_levels) {
      sample_ids <- metadata$sample_id[metadata$included & metadata$habitat == habitat_name &
        metadata$branch_context == context]
      if (length(sample_ids) > 0) {
        current_counts <- counts[selected, sample_ids, drop = FALSE]
        asv_rows[[length(asv_rows) + 1]] <- tibble(asv_id = selected, source_branch = branch_name,
          habitat = habitat_name, branch_context = context, samples = length(sample_ids),
          samples_detected = rowSums(current_counts > 0), detection_fraction = rowMeans(current_counts > 0),
          mean_read_percent = 100 * rowMeans(current_counts) / target_depth)
      }
    }
  }
}
branch_asvs <- bind_rows(asv_rows) %>% left_join(taxonomy %>% select(asv_id, Phylum, Genus), by = "asv_id")
write_csv(branch_asvs, file.path(results_dir, "branch_signature_asv_contexts.csv"))
local_contrasts <- branch_asvs %>% filter(branch_context %in% c("Branch 47", "Branch 66")) %>%
  mutate(local = branch_context == paste("Branch", source_branch)) %>%
  select(-branch_context) %>% pivot_wider(names_from = local,
    values_from = c(samples, samples_detected, detection_fraction, mean_read_percent), names_prefix = "local_") %>%
  mutate(detection_difference = detection_fraction_local_TRUE - detection_fraction_local_FALSE,
    read_percent_difference = mean_read_percent_local_TRUE - mean_read_percent_local_FALSE)
write_csv(local_contrasts, file.path(results_dir, "branch_signature_asv_contrasts.csv"))
write_csv(local_contrasts %>% group_by(habitat, source_branch) %>% summarize(asvs = n(),
  detected_either_branch = sum(samples_detected_local_TRUE > 0 | samples_detected_local_FALSE > 0),
  more_detection_own_branch = sum(detection_difference > 0),
  more_detection_other_branch = sum(detection_difference < 0),
  equal_detection = sum(detection_difference == 0),
  more_abundance_own_branch = sum(read_percent_difference > 0), .groups = "drop"),
  file.path(results_dir, "branch_signature_localization_summary.csv"))

# 5. Plot drainage footprints and the independently soil-defined signatures ------
if (make_branch_figures) {
  geometry <- read_csv(input_paths[8], skip = 1,
    col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
    col_types = cols(.default = col_double()), show_col_types = FALSE)
  lines <- geometry %>% group_by(segment) %>% filter(row_number() %% 12 == 1 | row_number() == n()) %>% ungroup()
  mapped_soils <- soils %>% filter(spatial_included) %>%
    left_join(source_soils %>% select(sample_id, branch), by = "sample_id") %>%
    mutate(branch = ifelse(is.na(branch), "Other / uncertain", branch))
  plot <- ggplot() + geom_sf(data = st_simplify(branches, dTolerance = 25),
    aes(fill = site_code), color = NA, alpha = 0.2) +
    geom_path(data = lines, aes(x, y, group = segment), color = "#AAAAAA", linewidth = 0.25) +
    geom_point(data = metadata %>% distinct(site_code, .keep_all = TRUE),
      aes(site_utm_x_m, site_utm_y_m), shape = 21, fill = "white", size = 1.7) +
    geom_point(data = mapped_soils, aes(soil_utm_x_m, soil_utm_y_m, color = branch), shape = 17, size = 3) +
    geom_text_repel(data = mapped_soils, aes(soil_utm_x_m, soil_utm_y_m,
      label = sub("HJA-2016_", "", site_code)), seed = label_seed, size = 3, max.overlaps = Inf) +
    geom_text_repel(data = metadata %>% filter(site_code %in% branch_outlet_sites) %>% distinct(site_code, .keep_all = TRUE),
      aes(site_utm_x_m, site_utm_y_m, label = site_code), seed = label_seed, size = 3.5) +
    scale_color_manual(values = c(source_colors, `Other / uncertain` = "#555555"), name = "Soil source branch") +
    scale_fill_manual(values = source_colors, name = "Drainage area to site") +
    coord_sf(crs = st_crs(26910), datum = st_crs(26910)) +
    labs(x = "UTM Zone 10N easting (m)", y = "Northing (m)")
  ggsave(here("figures", "2016_soil_branch_drainage_map.pdf"), plot, width = 9, height = 7, bg = "white")
  for (habitat_name in habitats) {
    plot <- profiles %>% filter(included, habitat == habitat_name, pool == "branch_concentrated") %>%
      mutate(branch_context = factor(branch_context, levels = group_levels)) %>%
      ggplot(aes(branch_context, source_read_percent, color = source_branch)) +
      geom_point(position = position_dodge(width = 0.35), size = 2.3, alpha = 0.8) +
      scale_color_manual(values = source_colors, name = "Soil source branch") +
      labs(x = "Aquatic drainage context", y = "Branch-concentrated soil ASV reads (%)") +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave(here("figures", paste0("2016_soil_branch_signatures_", habitat_name, ".pdf")),
      plot, width = 8, height = 5, bg = "white")
  }
}

input_checks$final_md5 <- unname(tools::md5sum(input_paths))
input_checks$unchanged <- input_checks$initial_md5 == input_checks$final_md5
stopifnot(all(input_checks$unchanged))
write_csv(input_checks, file.path(results_dir, "branch_input_checksums.csv"))
print(feature_table %>% group_by(branch) %>% summarize(
  branch_detected_asvs = sum(branch_detected), branch_concentrated_asvs = sum(branch_concentrated)))
print(summary %>% filter(pool == "branch_concentrated"), n = Inf)
