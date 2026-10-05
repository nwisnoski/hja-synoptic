# Audit the data available for sediment ASV-molecular bipartite networks.
# This script does not infer edges, fit modules, or run network null models.
library(here)
library(readr)
library(dplyr)
library(tidyr)

# 1. Settings ------------------------------------------------------------------
results_dir <- here("results", "fticr_2016", "bipartite_feasibility")
target_depth <- 10000
prevalence_design <- tibble(minimum_present_and_absent = c(5L, 7L, 10L, 11L),
  asvs = NA_integer_, molecular_features = NA_integer_, candidate_pairs = NA_integer_,
  distinct_asv_patterns = NA_integer_, distinct_molecular_patterns = NA_integer_)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"))

# 2. Record the intersection and preserve all input files ------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
stopifnot(!anyNA(checks$initial_md5))
metadata <- read_csv(input_paths[1],
  col_types = cols(site_code = col_character()), show_col_types = FALSE) %>%
  filter(included_10k, habitat == "sediment") %>%
  select(sample_id, site_code, site_stream_order, site_segment_number,
    site_drainage_area_ha, site_utm_x_m, site_utm_y_m)
molecular <- read_csv(input_paths[3],
  col_types = cols(site_code = col_character()), show_col_types = FALSE)
stopifnot(!anyDuplicated(metadata$site_code), !anyDuplicated(metadata$sample_id),
  !anyDuplicated(molecular$site_code))
inclusion <- full_join(
  metadata %>% transmute(site_code, sample_id, sediment_10k = TRUE),
  molecular %>% transmute(site_code, molecular_profile = TRUE), by = "site_code") %>%
  mutate(sediment_10k = replace_na(sediment_10k, FALSE),
    molecular_profile = replace_na(molecular_profile, FALSE),
    included = sediment_10k & molecular_profile,
    reason = case_when(included ~ "paired sediment observations",
      !molecular_profile ~ "no molecular profile",
      !sediment_10k ~ "no sediment sample retained at 10K"))
write_csv(inclusion, file.path(results_dir, "site_inclusion_audit.csv"))
site <- metadata %>% filter(site_code %in% molecular$site_code) %>% arrange(site_code)
stopifnot(!anyNA(site))

# 3. Audit incidence rather than interpret millions of candidate edges -----------
count_table <- read_csv(input_paths[2], show_col_types = FALSE)
stopifnot(!anyDuplicated(count_table$asv_id), all(site$sample_id %in% names(count_table)))
counts <- t(as.matrix(count_table[, site$sample_id]))
colnames(counts) <- count_table$asv_id
rownames(counts) <- site$site_code
pa <- as.matrix(molecular[match(site$site_code, molecular$site_code), -1])
rownames(pa) <- site$site_code
stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)),
  all(rowSums(counts) == target_depth), all(pa %in% c(0, 1)),
  identical(rownames(counts), rownames(pa)))
asv_prevalence <- colSums(counts > 0)
molecular_prevalence <- colSums(pa)
prevalence <- bind_rows(
  tibble(node_type = "ASV", feature_id = colnames(counts), present = asv_prevalence),
  tibble(node_type = "Molecular feature", feature_id = colnames(pa),
    present = molecular_prevalence)) %>%
  mutate(absent = nrow(site) - present, prevalence_fraction = present / nrow(site))
write_csv(prevalence, file.path(results_dir, "feature_prevalence.csv"))
site <- site %>% mutate(asv_richness = rowSums(counts > 0),
  detected_molecular_features = rowSums(pa))
write_csv(site, file.path(results_dir, "matched_site_manifest.csv"))
write_csv(site %>% count(site_stream_order, name = "sites"),
  file.path(results_dir, "stream_order_counts.csv"))
write_csv(site %>% count(site_segment_number, name = "sites"),
  file.path(results_dir, "segment_counts.csv"))
summary <- tibble(paired_sites = nrow(site),
  distinct_segments = n_distinct(site$site_segment_number),
  detected_asvs = sum(asv_prevalence > 0),
  detected_molecular_features = sum(molecular_prevalence > 0),
  all_detected_candidate_pairs = sum(asv_prevalence > 0) * sum(molecular_prevalence > 0),
  minimum_molecular_detections = min(rowSums(pa)),
  maximum_molecular_detections = max(rowSums(pa)),
  richness_spearman_rho = cor(site$asv_richness, site$detected_molecular_features,
    method = "spearman"))
write_csv(summary, file.path(results_dir, "data_summary.csv"))

# These thresholds describe feasibility for a binary-binary network.
# They are not a fitted-model grid or a selected final analysis threshold.
for (run in seq_len(nrow(prevalence_design))) {
  minimum <- prevalence_design$minimum_present_and_absent[run]
  keep_asv <- asv_prevalence >= minimum & asv_prevalence <= nrow(site) - minimum
  keep_molecular <- molecular_prevalence >= minimum &
    molecular_prevalence <= nrow(site) - minimum
  prevalence_design$asvs[run] <- sum(keep_asv)
  prevalence_design$molecular_features[run] <- sum(keep_molecular)
  prevalence_design$candidate_pairs[run] <- sum(keep_asv) * sum(keep_molecular)
  # Identical detection patterns cannot support distinct incidence-based links.
  prevalence_design$distinct_asv_patterns[run] <- nrow(unique(t(counts[, keep_asv] > 0)))
  prevalence_design$distinct_molecular_patterns[run] <- nrow(unique(t(pa[, keep_molecular])))
}
write_csv(prevalence_design, file.path(results_dir, "prevalence_feasibility.csv"))

# 4. Verify unchanged inputs and retain the software record -----------------------
checks$final_md5 <- unname(tools::md5sum(input_paths))
checks$unchanged <- checks$initial_md5 == checks$final_md5
write_csv(checks, file.path(results_dir, "input_integrity.csv"))
stopifnot(all(checks$unchanged))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(results_dir, "session_info.txt"))
print(summary)
print(prevalence_design)
