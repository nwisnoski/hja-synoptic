# Describe soil–aquatic ASV sharing along geographic distance, without routing.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

# 1. Settings and plotting theme ------------------------------------------------
target_depth <- 10000
make_distance_correlation_figure <- FALSE
habitats <- c("planktonic", "hyporheic", "sediment")
pools <- c("soil_detected", "soil_enriched")
detection_thresholds <- c(1, 3)
minimum_aquatic_detections <- 3
soil_enrichment_ratio <- 5
minimum_soil_reads_for_enrichment <- 5
projected_crs <- 26910
results_dir <- here("results", "diversity_2016", "soil_stream_distance")
derived_dir <- here("data", "derived", "soil_stream_distance_2016")
habitat_colors <- c(planktonic = "#0072B2", hyporheic = "#E69F00", sediment = "#009E73")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
theme_set(plot_theme)

# 2. Read the fixed table and GPS metadata; no drainage input is required --------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_metadata_screened.csv"),
  here("data", "StreamCenterline.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
count_table <- read_csv(input_paths[1], show_col_types = FALSE)
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()), show_col_types = FALSE)
taxonomy <- read_csv(input_paths[3], show_col_types = FALSE)
counts <- as.matrix(count_table[, -1])
rownames(counts) <- count_table$asv_id
metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
stopifnot(identical(metadata$sample_id, colnames(counts)), all(metadata$included_10k),
  !anyDuplicated(metadata$sample_id), !anyDuplicated(rownames(counts)),
  all(colSums(counts) == target_depth), all(counts >= 0), all(counts == floor(counts)))
soil_metadata <- metadata %>% filter(habitat == "soil") %>%
  select(sample_id, site_code, mapping_status, soil_longitude, soil_latitude) %>%
  mutate(included = !is.na(soil_longitude) & !is.na(soil_latitude) & mapping_status == "MATCHED",
    exclusion_reason = case_when(mapping_status != "MATCHED" ~ "Unresolved soil identity",
      is.na(soil_longitude) | is.na(soil_latitude) ~ "Missing GPS coordinates", TRUE ~ ""))
write_csv(soil_metadata, file.path(results_dir, "soil_inclusion_audit.csv"))
mapped_soils <- soil_metadata %>% filter(included)
soil_points <- st_as_sf(mapped_soils, coords = c("soil_longitude", "soil_latitude"),
  crs = 4326, remove = FALSE) %>% st_transform(projected_crs)
mapped_soils$soil_utm_x_m <- st_coordinates(soil_points)[, 1]
mapped_soils$soil_utm_y_m <- st_coordinates(soil_points)[, 2]
aquatic_metadata <- metadata %>% filter(habitat %in% habitats) %>%
  select(sample_id, site_code, habitat, site_stream_order, site_segment_number,
    site_drainage_area_ha, site_utm_x_m, site_utm_y_m)
stopifnot(!anyNA(aquatic_metadata), nrow(mapped_soils) == 13, nrow(aquatic_metadata) == 85)

# Distinguish distance to any mapped channel from distance to a sampled site.
# The nearest densely sampled centerline point approximates nearest-channel
# distance; it supplies an audit, not the predictor for sample-level sharing.
geometry <- read_csv(input_paths[4], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()), show_col_types = FALSE)
mapped_soils$nearest_network_point_distance_m <- NA_real_
mapped_soils$nearest_network_segment <- NA_real_
distance <- matrix(NA_real_, nrow(aquatic_metadata), nrow(mapped_soils),
  dimnames = list(aquatic_metadata$sample_id, mapped_soils$sample_id))
for (j in seq_len(nrow(mapped_soils))) {
  network_distance <- sqrt((geometry$x - mapped_soils$soil_utm_x_m[j])^2 +
    (geometry$y - mapped_soils$soil_utm_y_m[j])^2)
  closest <- which.min(network_distance)
  mapped_soils$nearest_network_point_distance_m[j] <- network_distance[closest]
  mapped_soils$nearest_network_segment[j] <- geometry$segment[closest]
  distance[, j] <- sqrt((aquatic_metadata$site_utm_x_m - mapped_soils$soil_utm_x_m[j])^2 +
    (aquatic_metadata$site_utm_y_m - mapped_soils$soil_utm_y_m[j])^2)
}
stopifnot(all(is.finite(distance)), all(distance > 0))
write_csv(mapped_soils, file.path(derived_dir, "soil_coordinate_nearness_audit.csv"))
write_csv(as.data.frame(distance) %>% mutate(sample_id = rownames(distance), .before = 1),
  file.path(derived_dir, "soil_aquatic_distances_m.csv"))
write_csv(aquatic_metadata, file.path(results_dir, "aquatic_inclusion_audit.csv"))

# 3. Calculate per-soil sharing and nearest candidate-soil distance for each ASV --
soil_counts <- counts[, metadata$sample_id[metadata$habitat == "soil"], drop = FALSE]
mapped_counts <- counts[, mapped_soils$sample_id, drop = FALSE]
aquatic_counts <- counts[, aquatic_metadata$sample_id, drop = FALSE]
soil_mean <- rowMeans(soil_counts)
aquatic_mean <- rowMeans(aquatic_counts)
pool_flags <- tibble(asv_id = rownames(counts), soil_detected = rowSums(soil_counts) > 0,
  soil_enriched = rowSums(soil_counts) >= minimum_soil_reads_for_enrichment &
    soil_mean / aquatic_mean >= soil_enrichment_ratio)
write_csv(pool_flags, file.path(results_dir, "asv_pool_definitions.csv"))
pair_rows <- list()
asv_rows <- list()
for (threshold in detection_thresholds) {
  soil_presence <- mapped_counts >= threshold
  water_presence <- aquatic_counts >= threshold
  nearest <- matrix(Inf, nrow(counts), ncol(aquatic_counts))
  for (j in seq_len(nrow(mapped_soils))) {
    source_asvs <- which(soil_presence[, j])
    candidate_distance <- matrix(rep(distance[, j], each = length(source_asvs)), nrow = length(source_asvs))
    nearest[source_asvs, ] <- pmin(nearest[source_asvs, , drop = FALSE], candidate_distance)
    for (pool_name in pools) {
      keep <- soil_presence[, j] & pool_flags[[pool_name]]
      stopifnot(sum(keep) > 0)
      pair_rows[[length(pair_rows) + 1]] <- aquatic_metadata %>%
        mutate(soil_sample_id = mapped_soils$sample_id[j], soil_site_code = mapped_soils$site_code[j],
          pool = pool_name, detection_threshold = threshold,
          distance_m = distance[, j], log2_distance_km = log2(distance_m / 1000),
          soil_asvs = sum(keep), shared_asvs = colSums(water_presence[keep, , drop = FALSE]),
          soil_asv_detection_fraction = shared_asvs / soil_asvs,
          aquatic_read_percent = 100 * colSums(aquatic_counts[keep, , drop = FALSE]) / target_depth)
    }
  }
  # A nearest candidate soil is any mapped soil in which this ASV is detected.
  # Require three aquatic detections, retain zeros elsewhere, and export flags
  # for taxa whose abundance or nearest-distance values cannot support a slope.
  for (habitat_name in habitats) {
    water_index <- which(aquatic_metadata$habitat == habitat_name)
    eligible <- which(rowSums(soil_presence) > 0 &
      rowSums(water_presence[, water_index, drop = FALSE]) >= minimum_aquatic_detections)
    for (asv_index in eligible) {
      current_distance <- nearest[asv_index, water_index]
      current_reads <- aquatic_counts[asv_index, water_index]
      diagnostic <- ""
      rho <- NA_real_
      if (length(unique(current_distance)) < 2 | length(unique(current_reads)) < 2) {
        diagnostic <- "Constant distance or abundance"
      } else {
        rho <- cor(current_distance, current_reads, method = "spearman")
      }
      asv_rows[[length(asv_rows) + 1]] <- tibble(asv_id = rownames(counts)[asv_index],
        habitat = habitat_name, detection_threshold = threshold,
        soil_enriched = pool_flags$soil_enriched[asv_index],
        mapped_soil_sources = sum(soil_presence[asv_index, ]), aquatic_samples = length(water_index),
        aquatic_samples_detected = sum(water_presence[asv_index, water_index]),
        spearman_nearest_distance_abundance = rho, diagnostic = diagnostic)
    }
  }
}
pairs <- bind_rows(pair_rows)
write_csv(pairs, file.path(results_dir, "soil_aquatic_distance_sharing.csv"))
asv_results <- bind_rows(asv_rows) %>% left_join(taxonomy %>% select(asv_id, Phylum, Genus), by = "asv_id")
write_csv(asv_results, file.path(results_dir, "asv_nearest_soil_distance_correlations.csv"))
write_csv(asv_results %>% group_by(habitat, detection_threshold, soil_enriched) %>%
  summarize(asvs = n(), estimable = sum(!is.na(spearman_nearest_distance_abundance)),
    negative_distance_correlations = sum(spearman_nearest_distance_abundance < 0, na.rm = TRUE),
    positive_distance_correlations = sum(spearman_nearest_distance_abundance > 0, na.rm = TRUE),
    median_rho = median(spearman_nearest_distance_abundance, na.rm = TRUE), .groups = "drop"),
  file.path(results_dir, "asv_distance_summary.csv"))

# 4. Describe source-specific correlations and distance slopes ------------------
# No independent-pair P values: sources, aquatic samples, and ASVs are reused.
responses <- c("soil_asv_detection_fraction", "aquatic_read_percent")
correlation_rows <- list()
model_rows <- list()
design <- expand_grid(habitat = habitats, pool = pools, detection_threshold = detection_thresholds,
  response = responses) %>% mutate(model_id = row_number())
write_csv(design, file.path(results_dir, "model_design.csv"))
for (i in seq_len(nrow(design))) {
  current <- pairs %>% filter(habitat == design$habitat[i], pool == design$pool[i],
    detection_threshold == design$detection_threshold[i])
  response <- design$response[i]
  # Two fixed intercept sets absorb source-pool and aquatic-site differences.
  # Report only the descriptive change per doubling of geographic distance.
  model <- lm(reformulate(c("factor(sample_id)", "factor(soil_sample_id)", "log2_distance_km"),
    response = response), data = current)
  effect <- unname(coef(model)["log2_distance_km"])
  model_rows[[i]] <- design[i, ] %>% mutate(pairs = nrow(current), soils = n_distinct(current$soil_sample_id),
    aquatic_samples = n_distinct(current$sample_id), effect_per_distance_doubling = effect,
    full_rank = model$rank == ncol(model.matrix(model)),
    diagnostic = ifelse(is.na(effect), "Nonestimable distance effect", "Descriptive; no calibrated P value"))
  for (soil_id in mapped_soils$sample_id) {
    source <- current %>% filter(soil_sample_id == soil_id)
    variable <- length(unique(source[[response]])) > 1
    rho <- NA_real_
    if (variable) {
      rho <- cor(source$distance_m, source[[response]], method = "spearman")
    }
    correlation_rows[[length(correlation_rows) + 1]] <- design[i, ] %>%
      mutate(soil_sample_id = soil_id, soil_site_code = source$soil_site_code[1],
        aquatic_samples = nrow(source), spearman_rho = rho,
        diagnostic = ifelse(variable, "", "Constant response"))
  }
}
source_correlations <- bind_rows(correlation_rows)
effects <- bind_rows(model_rows)
write_csv(source_correlations, file.path(results_dir, "individual_soil_distance_correlations.csv"))
write_csv(effects, file.path(results_dir, "adjusted_distance_effects.csv"))

# 5. Show variation among soils and preserve provenance --------------------------
if (make_distance_correlation_figure) {
  plot <- source_correlations %>% filter(pool == "soil_enriched", detection_threshold == 1,
    response == "soil_asv_detection_fraction") %>%
    mutate(soil_site_code = sub("HJA-2016_", "", soil_site_code)) %>%
    ggplot(aes(spearman_rho, soil_site_code, color = habitat)) +
    geom_vline(xintercept = 0, color = "#999999", linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.55), size = 2.5) +
    scale_color_manual(values = habitat_colors) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(x = "Spearman correlation: distance vs. soil ASVs detected", y = "Soil profile", color = "Habitat")
  ggsave(here("figures", "2016_soil_stream_distance_correlations.pdf"), plot,
    width = 9, height = 6, bg = "white")
}

checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
stopifnot(all(checks$unchanged))
write_csv(checks, file.path(results_dir, "input_checksums.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))
print(effects %>% filter(pool == "soil_enriched", detection_threshold == 1), width = Inf)
