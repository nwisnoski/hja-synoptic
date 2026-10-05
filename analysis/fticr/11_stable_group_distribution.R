# Describe the identities and catchment distribution of a completed stable group.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
core_id <- 2L
target_depth <- 10000L
label_seed <- 20261044L
expected_asvs <- 38L
expected_asv_patterns <- 37L
expected_molecules <- 8L
results_dir <- here("results", "fticr_2016", "stable_group_distribution")
figures_dir <- here("figures")
network_dir <- here("results", "fticr_2016", "bipartite_networks")
figure_prefix <- "2016_sediment_stable_group"
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
    legend.position = "bottom")
input_paths <- c(
  file.path(network_dir, "stable_core_asv_taxonomy.csv"),
  file.path(network_dir, "stable_core_molecular_properties.csv"),
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "network_dispersion", "site_network_position_audit.csv"),
  here("data", "StreamCenterline.csv"))
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
stopifnot(!anyNA(checks$initial_md5))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Retain the existing identities; do not infer a new group --------------------
asvs <- read_csv(input_paths[1], show_col_types = FALSE) %>% filter(core == core_id)
molecules <- read_csv(input_paths[2], show_col_types = FALSE) %>% filter(core == core_id)
metadata <- read_csv(input_paths[3], col_types = cols(site_code = col_character()),
  show_col_types = FALSE) %>% filter(included_10k, habitat == "sediment")
counts <- read_csv(input_paths[4], show_col_types = FALSE)
molecular <- read_csv(input_paths[5], col_types = cols(site_code = col_character()),
  show_col_types = FALSE)
stopifnot(nrow(asvs) == expected_asvs,
  n_distinct(asvs$pattern_id) == expected_asv_patterns,
  nrow(molecules) == expected_molecules,
  !anyDuplicated(asvs$feature_id), !anyDuplicated(molecules$feature_id),
  !anyDuplicated(metadata$site_code), !anyDuplicated(molecular$site_code),
  all(asvs$feature_id %in% counts$asv_id),
  all(molecules$feature_id %in% names(molecular)),
  all(metadata$sample_id %in% names(counts)))
stopifnot(all(colSums(as.matrix(counts[metadata$sample_id])) == target_depth))
core_counts <- as.matrix(counts[match(asvs$feature_id, counts$asv_id), metadata$sample_id])
rownames(core_counts) <- asvs$feature_id
representatives <- asvs %>% distinct(pattern_id, .keep_all = TRUE)
core_pattern_counts <- core_counts[representatives$feature_id, , drop = FALSE]
for (pattern in unique(asvs$pattern_id)) {
  members <- asvs$feature_id[asvs$pattern_id == pattern]
  # Pattern equivalence was established on 33 paired sites. Check that subset;
  # the extra microbial-only site need not share the same collapsed pattern.
  paired_samples <- metadata$sample_id[metadata$site_code %in% molecular$site_code]
  occurrence <- core_counts[members, paired_samples, drop = FALSE] > 0
  stopifnot(all(occurrence == occurrence[rep(1, nrow(occurrence)), , drop = FALSE]))
}
asvs <- asvs %>% mutate(detected_sediment_sites_34 = rowSums(core_counts > 0),
  reads_sediment_sites_34 = rowSums(core_counts),
  fraction_core_reads_34 = reads_sediment_sites_34 / sum(core_counts))
write_csv(asvs, file.path(results_dir, "core_asv_identities.csv"))
write_csv(asvs %>% count(Phylum, name = "asvs", sort = TRUE),
  file.path(results_dir, "core_phylum_summary.csv"))
write_csv(asvs %>% count(Genus, name = "asvs", sort = TRUE),
  file.path(results_dir, "core_genus_summary.csv"))

# Formulas are the workbook assignments, not identified structures. Preserve
# its Candidates field rather than treating every assignment as unique.
molecules$formula <- ""
element_columns <- c(C = "carbon_count", H = "hydrogen_count", N = "nitrogen_count",
  O = "oxygen_count", P = "phosphorus_count", S = "sulfur_count", Na = "sodium_count")
for (element in names(element_columns)) {
  number <- molecules[[element_columns[[element]]]]
  stopifnot(!anyNA(number), all(number >= 0), all(number == round(number)))
  molecules$formula <- paste0(molecules$formula,
    ifelse(number == 0, "", paste0(element, ifelse(number == 1, "", number))))
}
molecules <- molecules %>% mutate(
  broad_signature = recode(molecular_class_source, Protein = "Protein-like",
    Lipid = "Lipid-like", ConHC = "Condensed-hydrocarbon-like"),
  multiple_formula_candidates = candidate_count > 1)
write_csv(molecules, file.path(results_dir, "core_molecular_identities.csv"))

# 3. Describe all available sites, preserving unmeasured profiles as missing ----
asv_site <- metadata %>% transmute(site_code, sample_id,
  core_asvs_detected = colSums(core_counts > 0),
  core_asv_reads = colSums(core_counts),
  core_asv_read_percent = 100 * core_asv_reads / target_depth)
pattern_site <- tibble(site_code = metadata$site_code,
  core_asv_patterns_detected = colSums(core_pattern_counts > 0))
# Retain patterns only on the inference sites. For the microbial-only site,
# report original ASV detections and reads without extrapolating equivalence.
pattern_site$core_asv_patterns_detected[!pattern_site$site_code %in% molecular$site_code] <- NA
asv_site <- left_join(asv_site, pattern_site, by = "site_code")
mol_incidence <- as.matrix(molecular[molecules$feature_id])
stopifnot(all(mol_incidence %in% c(0, 1)))
molecular_site <- molecular %>% transmute(site_code,
  core_molecules_detected = rowSums(mol_incidence))
network <- read_csv(input_paths[6], col_types = cols(site_code = col_character()),
  show_col_types = FALSE) %>% select(site_code, network_group, site_stream_order,
    site_segment_number, site_utm_x_m, site_utm_y_m, site_drainage_area_ha)
extra <- metadata %>% anti_join(network, by = "site_code") %>%
  select(site_code, site_stream_order, site_segment_number, site_utm_x_m,
    site_utm_y_m, site_drainage_area_ha) %>%
  mutate(network_group = case_when(site_stream_order <= 2 ~ "Headwater",
    site_stream_order <= 4 ~ "Intermediate", site_stream_order == 5 ~ "Mainstem"))
network <- bind_rows(network, extra)
geometry <- read_csv(input_paths[7], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
stopifnot(nrow(problems(geometry)) == 0, !anyDuplicated(network$site_code),
  !anyNA(network), all(is.finite(geometry$x)), all(is.finite(geometry$y)))
network$centerline_match_error_m <- NA_real_
network$basin_outlet_stream_distance_m <- NA_real_
network$centerline_stream_order <- NA_real_
for (i in seq_len(nrow(network))) {
  current <- geometry %>% filter(segment == network$site_segment_number[i])
  stopifnot(nrow(current) > 0)
  distance <- sqrt((current$x - network$site_utm_x_m[i])^2 +
    (current$y - network$site_utm_y_m[i])^2)
  closest <- which.min(distance)
  network$centerline_match_error_m[i] <- distance[closest]
  network$basin_outlet_stream_distance_m[i] <- current$outlet_distance[closest]
  network$centerline_stream_order[i] <- current$order[closest]
}
stopifnot(all(network$centerline_match_error_m < 1),
  all(network$centerline_stream_order == network$site_stream_order),
  all(is.finite(network$basin_outlet_stream_distance_m)))
write_csv(network, file.path(results_dir, "site_network_audit.csv"))
site <- full_join(asv_site, molecular_site, by = "site_code") %>%
  left_join(network, by = "site_code") %>%
  mutate(microbial_profile_available = !is.na(core_asvs_detected),
    molecular_profile_available = !is.na(core_molecules_detected),
    inference_site = microbial_profile_available & molecular_profile_available,
    core_asv_fraction = core_asvs_detected / expected_asvs,
    core_asv_pattern_fraction = core_asv_patterns_detected / expected_asv_patterns,
    core_molecular_fraction = core_molecules_detected / expected_molecules,
    network_group = factor(network_group, levels = c("Headwater", "Intermediate", "Mainstem")))
stopifnot(nrow(site) == 45, sum(site$inference_site) == 33,
  sum(site$microbial_profile_available) == 34,
  sum(site$molecular_profile_available) == 44, !anyNA(site$site_utm_x_m))
write_csv(site, file.path(results_dir, "core_site_profiles.csv"))
summary <- site %>% group_by(network_group) %>% summarize(
  microbial_sites = sum(microbial_profile_available),
  molecular_sites = sum(molecular_profile_available), paired_sites = sum(inference_site),
  median_asvs_detected = median(core_asvs_detected, na.rm = TRUE),
  median_asv_read_percent = median(core_asv_read_percent, na.rm = TRUE),
  median_molecules_detected = median(core_molecules_detected, na.rm = TRUE),
  paired_median_asvs_detected = median(core_asvs_detected[inference_site]),
  paired_median_asv_read_percent = median(core_asv_read_percent[inference_site]),
  paired_median_molecules_detected = median(core_molecules_detected[inference_site]),
  .groups = "drop")
write_csv(summary, file.path(results_dir, "core_network_group_summary.csv"))
asv_long <- as.data.frame(t(core_counts)) %>% mutate(sample_id = rownames(.)) %>%
  pivot_longer(-sample_id, names_to = "feature_id", values_to = "reads") %>%
  left_join(metadata %>% select(sample_id, site_code), by = "sample_id") %>%
  mutate(detected = reads > 0)
write_csv(asv_long, file.path(results_dir, "core_asv_site_detections.csv"))
write_csv(molecular %>% select(site_code, all_of(molecules$feature_id)) %>%
  pivot_longer(-site_code, names_to = "feature_id", values_to = "detected"),
  file.path(results_dir, "core_molecular_site_detections.csv"))

# 4. Map the microbial and molecular components separately ----------------------
# Use original ASV fractions on all 34 microbial sites. This avoids imputing a
# paired-site occurrence-pattern equivalence at the unmatched site 43.
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% 10 == 1 | row_number() == n()) %>% ungroup()
map_rows <- bind_rows(
  site %>% mutate(component = "ASVs", detected_fraction = core_asv_fraction),
  site %>% mutate(component = "Molecular features", detected_fraction = core_molecular_fraction))
maps <- list()
for (component_name in c("ASVs", "Molecular features")) {
  current <- map_rows %>% filter(component == component_name)
  map <- ggplot(map_lines, aes(x, y, group = segment)) +
    geom_path(color = "#C0C0C0", linewidth = 0.35) +
    geom_point(data = current %>% filter(!is.na(detected_fraction)),
      aes(site_utm_x_m, site_utm_y_m, fill = detected_fraction),
      inherit.aes = FALSE, shape = 21, size = 3.3, color = "black", stroke = 0.4) +
    geom_point(data = current %>% filter(is.na(detected_fraction)),
      aes(site_utm_x_m, site_utm_y_m, shape = "Profile unavailable"),
      inherit.aes = FALSE, color = "#666666", size = 2.5) +
    geom_text_repel(data = current, aes(site_utm_x_m, site_utm_y_m, label = site_code),
      inherit.aes = FALSE, seed = label_seed, size = 2.8, max.overlaps = Inf,
      box.padding = 0.3, point.padding = 0.25, min.segment.length = 0,
      segment.color = "#999999", segment.size = 0.2) +
    scale_fill_viridis_c(option = "C", direction = -1, limits = c(0, 1),
      breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"),
      name = ifelse(component_name == "ASVs", "Fraction of 38 ASVs detected",
        "Fraction of 8 molecular features detected")) +
    scale_shape_manual(values = c("Profile unavailable" = 4), name = NULL) +
    coord_equal() + labs(x = "UTM easting (m)", y = "UTM northing (m)") +
    plot_theme + theme(panel.grid = element_blank())
  maps[[component_name]] <- map
  suffix <- ifelse(component_name == "ASVs", "asv_map", "molecular_map")
  ggsave(file.path(figures_dir, paste0(figure_prefix, "_", suffix, ".pdf")), map,
    width = 8, height = 7.5, bg = "white")
}
combined <- maps[["ASVs"]] + maps[["Molecular features"]] +
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")
ggsave(file.path(figures_dir, paste0(figure_prefix, "_catchment_maps.pdf")), combined,
  width = 15, height = 7.5, bg = "white")

# 5. Verify input integrity and record the scope --------------------------------
write_csv(tibble(core = core_id, original_asvs = expected_asvs,
  paired_asv_patterns = expected_asv_patterns, molecular_features = expected_molecules,
  microbial_sites = sum(site$microbial_profile_available),
  molecular_sites = sum(site$molecular_profile_available), paired_sites = sum(site$inference_site),
  inferential_tests_run = FALSE, label_seed = label_seed),
  file.path(results_dir, "settings.csv"))
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
write_csv(checks, file.path(results_dir, "input_integrity.csv"))
stopifnot(all(checks$unchanged))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(results_dir, "session_info.txt"))
print(summary, width = Inf)
print(site %>% filter(inference_site) %>% arrange(desc(core_asv_read_percent)) %>%
  select(site_code, network_group, core_asvs_detected, core_asv_read_percent,
    core_molecules_detected), n = 33)
