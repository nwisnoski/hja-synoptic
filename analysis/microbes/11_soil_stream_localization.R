# Describe soil-detected ASVs and their localization in contributing drainage areas.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(sf)

# 1. Settings and plotting theme ------------------------------------------------
target_depth <- 10000
make_asv_example_figures <- FALSE
habitats <- c("planktonic", "hyporheic", "sediment")
detection_thresholds <- c(1, 3)
minimum_contrast_samples <- 3
soil_enrichment_ratio <- 5
minimum_soil_reads_for_enrichment <- 5
label_seed <- 20261005
results_dir <- here("results", "diversity_2016", "soil_stream_localization")
routing_dir <- here("data", "derived", "soil_drainage_2016")
figures_dir <- here("figures")
habitat_colors <- c(planktonic = "#0072B2", hyporheic = "#E69F00", sediment = "#009E73")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
theme_set(plot_theme)

# 2. Validate fixed-10K counts and explicit soil identity exclusions -------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_metadata_screened.csv"),
  file.path(routing_dir, "soil_site_membership.csv"),
  file.path(routing_dir, "soil_coordinate_audit.csv"),
  file.path(routing_dir, "aquatic_catchment_audit.csv"),
  here("data", "StreamCenterline.csv"),
  here("data", "spatial", "soil_routing", "hja_gaged_watersheds.geojson"))
input_checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
count_table <- read_csv(input_paths[1], show_col_types = FALSE)
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()), show_col_types = FALSE)
taxonomy <- read_csv(input_paths[3], show_col_types = FALSE)
counts <- as.matrix(count_table[, -1])
rownames(counts) <- count_table$asv_id
metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
stopifnot(!anyNA(metadata$sample_id), identical(metadata$sample_id, colnames(counts)),
  all(metadata$included_10k), all(colSums(counts) == target_depth),
  !anyDuplicated(rownames(counts)), all(counts >= 0), all(counts == floor(counts)))
soil_metadata <- metadata %>% filter(habitat == "soil")
aquatic_metadata <- metadata %>% filter(habitat %in% habitats)
soils <- read_csv(input_paths[5], show_col_types = FALSE)
mapped_soils <- soils %>% filter(spatial_included)
membership <- read_csv(input_paths[4], col_types = cols(aquatic_site_code = col_character()), show_col_types = FALSE)
catchment_audit <- read_csv(input_paths[6], col_types = cols(site_code = col_character()), show_col_types = FALSE)
stopifnot(nrow(soils) == nrow(soil_metadata), !anyDuplicated(mapped_soils$sample_id),
  all(mapped_soils$sample_id %in% soil_metadata$sample_id),
  !anyDuplicated(paste(membership$soil_sample_id, membership$aquatic_site_code)))
soil_counts <- counts[, soil_metadata$sample_id, drop = FALSE]
mapped_counts <- counts[, mapped_soils$sample_id, drop = FALSE]
aquatic_counts <- counts[, aquatic_metadata$sample_id, drop = FALSE]
write_csv(metadata %>% transmute(sample_id, site_code, habitat,
  general_sharing_included = TRUE,
  spatial_coordinate_available = ifelse(habitat == "soil", sample_id %in% mapped_soils$sample_id, TRUE)) %>%
  left_join(soils %>% select(sample_id, exclusion_reason), by = "sample_id"),
  file.path(results_dir, "sample_inclusion_audit.csv"))

# Independent point-in-polygon checks establish the known WS1 soil subset.
gaged <- st_read(input_paths[8], quiet = TRUE) %>% st_transform(26910)
soil_points <- st_as_sf(mapped_soils, coords = c("soil_utm_x_m", "soil_utm_y_m"), crs = 26910, remove = FALSE)
polygon_audit <- st_join(soil_points, gaged %>% select(WS_)) %>% st_drop_geometry()
write_csv(polygon_audit, file.path(results_dir, "gaged_watershed_soil_audit.csv"))
ws1_soils <- polygon_audit$sample_id[!is.na(polygon_audit$WS_) & polygon_audit$WS_ == "1"]
ws1_sites <- aquatic_metadata %>% filter(grepl("^WS1-", site_code)) %>% distinct(site_code)
ws1_checks <- membership %>% filter(soil_sample_id %in% ws1_soils,
  aquatic_site_code %in% ws1_sites$site_code)
write_csv(ws1_checks, file.path(results_dir, "ws1_membership_crosscheck.csv"))

# 3. General sharing includes all 15 soils, even those without mapped locations ---
asv_summary <- tibble(asv_id = rownames(counts), soil_samples_detected = rowSums(soil_counts > 0),
  mapped_soil_samples_detected = rowSums(mapped_counts > 0),
  soil_reads = rowSums(soil_counts), soil_mean_read_percent = 100 * rowMeans(soil_counts) / target_depth,
  aquatic_samples_detected = rowSums(aquatic_counts > 0),
  aquatic_mean_read_percent = 100 * rowMeans(aquatic_counts) / target_depth,
  mainstem_samples_detected = rowSums(aquatic_counts[, aquatic_metadata$site_stream_order == 5, drop = FALSE] > 0)) %>%
  mutate(soil_detected = soil_samples_detected > 0,
    soil_to_aquatic_mean_ratio = soil_mean_read_percent / aquatic_mean_read_percent,
    soil_enriched = soil_reads >= minimum_soil_reads_for_enrichment &
      soil_to_aquatic_mean_ratio >= soil_enrichment_ratio) %>%
  left_join(taxonomy %>% select(asv_id, Kingdom, Phylum, Class, Order, Family, Genus), by = "asv_id")
write_csv(asv_summary, file.path(results_dir, "soil_asv_identities.csv"))
overview <- tibble(total_asvs = nrow(counts), soil_samples = ncol(soil_counts),
  mapped_soil_samples = ncol(mapped_counts), aquatic_samples = ncol(aquatic_counts),
  soil_detected_asvs = sum(asv_summary$soil_detected),
  soil_detected_asvs_in_aquatic = sum(asv_summary$soil_detected & asv_summary$aquatic_samples_detected > 0),
  soil_detected_asvs_in_mainstem = sum(asv_summary$soil_detected & asv_summary$mainstem_samples_detected > 0),
  soil_enriched_asvs = sum(asv_summary$soil_enriched),
  soil_enriched_asvs_in_aquatic = sum(asv_summary$soil_enriched & asv_summary$aquatic_samples_detected > 0))
write_csv(overview, file.path(results_dir, "sharing_overview.csv"))
site_rows <- list()
habitat_rows <- list()
for (habitat_name in habitats) {
  current_metadata <- aquatic_metadata %>% filter(habitat == habitat_name)
  current_counts <- counts[, current_metadata$sample_id, drop = FALSE]
  for (pool_name in c("soil_detected", "soil_enriched")) {
    keep <- asv_summary[[pool_name]]
    habitat_rows[[length(habitat_rows) + 1]] <- tibble(habitat = habitat_name, pool = pool_name,
      samples = nrow(current_metadata), soil_pool_asvs = sum(keep),
      soil_pool_asvs_detected = sum(rowSums(current_counts[keep, , drop = FALSE]) > 0),
      median_read_percent = median(100 * colSums(current_counts[keep, , drop = FALSE]) / target_depth))
    site_rows[[length(site_rows) + 1]] <- current_metadata %>% select(sample_id, site_code, habitat,
      site_stream_order, site_segment_number, site_drainage_area_ha, site_utm_x_m, site_utm_y_m) %>%
      mutate(pool = pool_name, soil_pool_asvs_detected = colSums(current_counts[keep, , drop = FALSE] > 0),
        soil_pool_read_percent = 100 * colSums(current_counts[keep, , drop = FALSE]) / target_depth)
  }
}
site_sharing <- bind_rows(site_rows)
write_csv(bind_rows(habitat_rows), file.path(results_dir, "habitat_sharing_summary.csv"))
write_csv(site_sharing, file.path(results_dir, "site_soil_asv_sharing.csv"))
write_csv(site_sharing %>% filter(site_stream_order == 5), file.path(results_dir, "mainstem_soil_asv_sharing.csv"))

# 4. Compare each individual soil profile with all aquatic samples ---------------
# Threshold 1 is primary. Threshold 3 checks sensitivity to low-count detections.
pair_rows <- list()
asv_contrast_rows <- list()
for (threshold in detection_thresholds) {
  soil_presence <- mapped_counts >= threshold
  water_presence <- aquatic_counts >= threshold
  for (j in seq_len(nrow(mapped_soils))) {
    soil_id <- mapped_soils$sample_id[j]
    for (pool_name in c("soil_detected", "soil_enriched")) {
      keep <- soil_presence[, j] & asv_summary[[pool_name]]
      current <- aquatic_metadata %>% select(sample_id, site_code, habitat, site_stream_order,
        site_segment_number, site_drainage_area_ha, site_utm_x_m, site_utm_y_m) %>%
        mutate(soil_sample_id = soil_id, soil_site_code = mapped_soils$site_code[j],
          detection_threshold = threshold, pool = pool_name, soil_asvs = sum(keep),
          shared_asvs = colSums(water_presence[keep, , drop = FALSE]),
          soil_asv_detection_fraction = shared_asvs / soil_asvs,
          aquatic_read_percent = 100 * colSums(aquatic_counts[keep, , drop = FALSE]) / target_depth,
          euclidean_distance_m = sqrt((site_utm_x_m - mapped_soils$soil_utm_x_m[j])^2 +
            (site_utm_y_m - mapped_soils$soil_utm_y_m[j])^2)) %>%
        left_join(membership %>% select(soil_sample_id, aquatic_site_code, contributing,
          contributing_fraction, membership_stable, routing_qc_pass, robust_spatial_included),
          by = c("soil_sample_id", "site_code" = "aquatic_site_code"))
      pair_rows[[length(pair_rows) + 1]] <- current
    }
  }
  # For each ASV, a contributing source means detection in at least one mapped
  # soil within that site's catchment. Multiple possible soils remain possible.
  connection <- matrix(FALSE, nrow(mapped_soils), nrow(aquatic_metadata))
  uncertain <- connection
  for (j in seq_len(nrow(mapped_soils))) {
    current <- membership %>% filter(soil_sample_id == mapped_soils$sample_id[j])
    current <- current[match(aquatic_metadata$site_code, current$aquatic_site_code), ]
    stopifnot(!anyNA(current$contributing))
    connection[j, ] <- current$contributing
    uncertain[j, ] <- !current$robust_spatial_included
  }
  local_source <- (soil_presence %*% connection) > 0
  ambiguous_source <- (soil_presence %*% uncertain) > 0
  # Only contrasts with stable membership for every detected candidate source
  # are used. This conservatively preserves spatial uncertainty.
  for (habitat_name in habitats) {
    water_index <- which(aquatic_metadata$habitat == habitat_name)
    for (pool_name in c("soil_detected", "soil_enriched")) {
      keep <- asv_summary[[pool_name]] & rowSums(soil_presence) > 0
      selected <- which(keep)
      for (asv_index in selected) {
        usable <- !ambiguous_source[asv_index, water_index]
        local <- local_source[asv_index, water_index] & usable
        other <- !local_source[asv_index, water_index] & usable
        n_local <- sum(local)
        n_other <- sum(other)
        if (n_local >= minimum_contrast_samples & n_other >= minimum_contrast_samples) {
          reads <- aquatic_counts[asv_index, water_index]
          detected <- water_presence[asv_index, water_index]
          asv_contrast_rows[[length(asv_contrast_rows) + 1]] <- tibble(
            asv_id = rownames(counts)[asv_index], habitat = habitat_name, pool = pool_name,
            detection_threshold = threshold, mapped_soil_sources = sum(soil_presence[asv_index, ]),
            contributing_samples = n_local, noncontributing_samples = n_other,
            contributing_detection_fraction = mean(detected[local]),
            noncontributing_detection_fraction = mean(detected[other]),
            detection_difference = mean(detected[local]) - mean(detected[other]),
            contributing_mean_read_percent = 100 * mean(reads[local]) / target_depth,
            noncontributing_mean_read_percent = 100 * mean(reads[other]) / target_depth,
            read_percent_difference = 100 * (mean(reads[local]) - mean(reads[other])) / target_depth)
        }
      }
    }
  }
}
pairs <- bind_rows(pair_rows)
write_csv(pairs, file.path(results_dir, "soil_aquatic_pair_sharing.csv"))
source_contrasts <- pairs %>% filter(robust_spatial_included) %>%
  group_by(detection_threshold, pool, habitat, soil_sample_id, soil_site_code, contributing) %>%
  summarize(samples = n(), segments = n_distinct(site_segment_number),
    detection_fraction = mean(soil_asv_detection_fraction), read_percent = mean(aquatic_read_percent),
    median_distance_m = median(euclidean_distance_m), .groups = "drop") %>%
  pivot_wider(names_from = contributing, values_from = c(samples, segments, detection_fraction,
    read_percent, median_distance_m), names_prefix = "contributing_") %>%
  mutate(contrast_available = !is.na(samples_contributing_TRUE) & !is.na(samples_contributing_FALSE),
    detection_difference = detection_fraction_contributing_TRUE - detection_fraction_contributing_FALSE,
    read_percent_difference = read_percent_contributing_TRUE - read_percent_contributing_FALSE)
write_csv(source_contrasts, file.path(results_dir, "individual_soil_contrasts.csv"))

# Compare within stream order to avoid treating mainstem accumulation as a
# headwater locality signal. Strata without both categories are explicitly absent.
order_strata <- pairs %>% filter(robust_spatial_included) %>%
  group_by(detection_threshold, pool, habitat, soil_sample_id, site_stream_order, contributing) %>%
  summarize(samples = n(), detection_fraction = mean(soil_asv_detection_fraction),
    read_percent = mean(aquatic_read_percent), .groups = "drop") %>%
  pivot_wider(names_from = contributing, values_from = c(samples, detection_fraction, read_percent),
    names_prefix = "contributing_") %>%
  mutate(contrast_available = !is.na(samples_contributing_TRUE) & !is.na(samples_contributing_FALSE),
    detection_difference = detection_fraction_contributing_TRUE - detection_fraction_contributing_FALSE,
    read_percent_difference = read_percent_contributing_TRUE - read_percent_contributing_FALSE)
write_csv(order_strata, file.path(results_dir, "within_order_soil_contrasts.csv"))
write_csv(order_strata %>% filter(contrast_available) %>%
  group_by(detection_threshold, pool, habitat) %>% summarize(soil_order_strata = n(),
    distinct_soils = n_distinct(soil_sample_id),
    mean_detection_difference = mean(detection_difference),
    mean_read_percent_difference = mean(read_percent_difference), .groups = "drop"),
  file.path(results_dir, "within_order_contrast_summary.csv"))
asv_contrasts <- bind_rows(asv_contrast_rows) %>%
  left_join(asv_summary %>% select(asv_id, Phylum, Genus, soil_samples_detected, soil_mean_read_percent), by = "asv_id")
write_csv(asv_contrasts, file.path(results_dir, "asv_localization_contrasts.csv"))
write_csv(asv_contrasts %>% group_by(detection_threshold, habitat, pool) %>%
  summarize(eligible_asvs = n(), asvs_detected_anywhere = sum(contributing_detection_fraction > 0 | noncontributing_detection_fraction > 0),
    more_detection_contributing = sum(detection_difference > 0),
    equal_detection = sum(detection_difference == 0),
    more_detection_noncontributing = sum(detection_difference < 0),
    more_abundant_contributing = sum(read_percent_difference > 0),
    median_detection_difference = median(detection_difference), .groups = "drop"),
  file.path(results_dir, "asv_localization_summary.csv"))

# 5. Descriptive adjustment for sample identity, source identity, and proximity ---
# These coefficients are exploratory effect estimates. No pair-independent P
# values are reported: pairs share sources, aquatic samples, ASVs, and segments.
model_rows <- list()
for (habitat_name in habitats) {
  for (pool_name in c("soil_detected", "soil_enriched")) {
    current <- pairs %>% filter(habitat == habitat_name, pool == pool_name,
      detection_threshold == 1, robust_spatial_included)
    for (response in c("soil_asv_detection_fraction", "aquatic_read_percent")) {
      model <- lm(reformulate(c("factor(sample_id)", "factor(soil_sample_id)",
        "log1p(euclidean_distance_m)", "contributing"), response = response), data = current)
      estimate <- unname(coef(model)["contributingTRUE"])
      model_rows[[length(model_rows) + 1]] <- tibble(habitat = habitat_name, pool = pool_name, response = response,
        pairs = nrow(current), soils = n_distinct(current$soil_sample_id),
        samples = n_distinct(current$sample_id), contributing_effect = estimate,
        full_rank = model$rank == ncol(model.matrix(model)),
        diagnostic = ifelse(is.na(estimate), "Nonestimable contributing effect", "Descriptive; no calibrated P value"))
    }
  }
}
write_csv(bind_rows(model_rows), file.path(results_dir, "adjusted_descriptive_effects.csv"))

# 6. Save maps, source contrasts, and individual-ASV examples ---------------------
geometry <- read_csv(input_paths[7], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()), show_col_types = FALSE)
map_lines <- geometry %>% group_by(segment) %>% filter(row_number() %% 12 == 1 | row_number() == n()) %>% ungroup()
for (habitat_name in habitats) {
  current <- site_sharing %>% filter(habitat == habitat_name, pool == "soil_detected")
  plot <- ggplot(map_lines, aes(x / 1000, y / 1000, group = segment)) +
    geom_path(color = "#B8B8B8", linewidth = 0.3) +
    geom_point(data = current, aes(site_utm_x_m / 1000, site_utm_y_m / 1000,
      fill = soil_pool_read_percent), inherit.aes = FALSE, shape = 21, size = 3) +
    geom_point(data = mapped_soils, aes(soil_utm_x_m / 1000, soil_utm_y_m / 1000),
      inherit.aes = FALSE, shape = 17, size = 2.5, color = "#D55E00") +
    geom_text_repel(data = mapped_soils, aes(soil_utm_x_m / 1000, soil_utm_y_m / 1000,
      label = sub("HJA-2016_", "", site_code)), inherit.aes = FALSE,
      size = 3, seed = label_seed, max.overlaps = Inf) +
    scale_fill_viridis_c(name = "Soil-detected\nASV reads (%)") +
    coord_equal() + labs(x = "UTM Zone 10N easting (km)", y = "Northing (km)")
  ggsave(file.path(figures_dir, paste0("2016_soil_asv_sharing_map_", habitat_name, ".pdf")),
    plot, width = 8, height = 7, bg = "white")
}
plot <- source_contrasts %>% filter(detection_threshold == 1, pool == "soil_detected", contrast_available) %>%
  mutate(soil_site_code = sub("HJA-2016_", "", soil_site_code)) %>%
  ggplot(aes(detection_difference * 100, soil_site_code, color = habitat)) +
  geom_vline(xintercept = 0, color = "#999999", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.55), size = 2.5) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "Contributing minus other sites: shared soil ASVs (percentage points)", y = "Soil profile", color = "Habitat")
ggsave(file.path(figures_dir, "2016_individual_soil_localization.pdf"), plot, width = 9, height = 6, bg = "white")

# Select examples by the contributing-minus-other abundance contrast, retaining
# both directions. These examples are selected descriptions, not tested taxa.
examples <- asv_contrasts %>% filter(detection_threshold == 1, pool == "soil_enriched") %>%
  group_by(habitat) %>% arrange(desc(abs(read_percent_difference)), .by_group = TRUE) %>%
  slice_head(n = 4) %>% ungroup()
write_csv(examples, file.path(results_dir, "selected_asv_examples.csv"))
example_rows <- list()
for (i in seq_len(nrow(examples))) {
  asv_id <- examples$asv_id[i]
  source_ids <- mapped_soils$sample_id[mapped_counts[asv_id, ] > 0]
  current <- pairs %>% filter(detection_threshold == 1, pool == "soil_detected", habitat == examples$habitat[i],
    soil_sample_id %in% source_ids) %>% group_by(sample_id, site_code, habitat,
      site_stream_order, site_drainage_area_ha, site_utm_x_m, site_utm_y_m) %>%
    summarize(contributing_source = any(contributing), robust = all(robust_spatial_included), .groups = "drop") %>%
    mutate(asv_id = .env$asv_id, reads = counts[.env$asv_id, sample_id], read_percent = 100 * reads / target_depth)
  example_rows[[length(example_rows) + 1]] <- current
}
example_sites <- bind_rows(example_rows)
write_csv(example_sites, file.path(results_dir, "selected_asv_site_profiles.csv"))
if (make_asv_example_figures) {
  for (habitat_name in habitats) {
    plot <- example_sites %>% filter(habitat == habitat_name, robust) %>%
      mutate(source = ifelse(contributing_source, "Contributing soil detected", "Other sampled soils")) %>%
      ggplot(aes(log10(site_drainage_area_ha), read_percent, color = source, shape = factor(site_stream_order))) +
      geom_point(size = 2.3) + facet_wrap(~asv_id, scales = "free_y", ncol = 2) +
      scale_color_manual(values = c("Contributing soil detected" = "#0072B2", "Other sampled soils" = "#D55E00")) +
      labs(x = "Log10 drainage area (ha)", y = "ASV reads (%)", color = "Source context", shape = "Stream order")
    ggsave(file.path(figures_dir, paste0("2016_soil_asv_examples_", habitat_name, ".pdf")),
      plot, width = 9, height = 7, bg = "white")
  }
}

# 7. Retain input-integrity and execution records --------------------------------
input_checks$final_md5 <- unname(tools::md5sum(input_paths))
input_checks$unchanged <- input_checks$initial_md5 == input_checks$final_md5
stopifnot(all(input_checks$unchanged))
write_csv(input_checks, file.path(results_dir, "analysis_input_checksums.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))
print(overview)
print(bind_rows(habitat_rows))
print(bind_rows(model_rows))
