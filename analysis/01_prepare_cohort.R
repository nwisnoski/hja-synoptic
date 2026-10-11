# Construct the reviewed soil-porewater-sediment cohort without new rarefaction.
library(here)
library(readr)
library(dplyr)

# 1. Settings ------------------------------------------------------------------
target_depth <- 10000L
output_dir <- here("data", "derived", "soil_porewater_sediment_2016")
habitats <- c("soil", "hyporheic", "sediment")
# Exact libraries identified in the October 10 review, not an ID-range filter.
flagged <- tibble::tribble(
  ~sample_id, ~expected_site, ~expected_habitat, ~expected_lane_sample,
  "hja2016_099", "162", "planktonic", "S1",
  "hja2016_100", "64", "planktonic", "S2",
  "hja2016_101", "58", "planktonic", "S3",
  "hja2016_102", "73", "hyporheic", "S4",
  "hja2016_103", "260", "planktonic", "S5",
  "hja2016_105", "178", "planktonic", "S6",
  "hja2016_106", "170", "planktonic", "S7",
  "hja2016_107", "57", "hyporheic", "S8",
  "hja2016_109", "169", "planktonic", "S9",
  "hja2016_110", "76", "hyporheic", "S10",
  "hja2016_111", "160", "planktonic", "S11",
  "hja2016_114", "WS1-5", "planktonic", "S60",
  "hja2016_116", "53", "planktonic", "S13",
  "hja2016_119", "56", "planktonic", "S14",
  "hja2016_121", "84", "hyporheic", "S15"
)

# 2. Validate identities and distinguish QC exclusion from manuscript scope -----
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv"),
  here("data", "derived", "microbial_diversity_2016", "asv_counts_10k.csv")
)
checks <- tibble(path = input_paths, initial_md5 = unname(tools::md5sum(input_paths)))
metadata <- read_csv(input_paths[1], col_types = cols(site_code = col_character()))
original <- read_csv(input_paths[2], show_col_types = FALSE)
stopifnot(!anyDuplicated(metadata$sample_id), !anyDuplicated(original$asv_id))
flag_metadata <- metadata[match(flagged$sample_id, metadata$sample_id), ]
stopifnot(identical(flag_metadata$sample_id, flagged$sample_id),
  identical(flag_metadata$site_code, flagged$expected_site),
  identical(flag_metadata$habitat, flagged$expected_habitat),
  identical(flag_metadata$lane_sample, flagged$expected_lane_sample),
  all(flag_metadata$plate == "Plate2"), all(flag_metadata$included_10k))
audit <- metadata %>% mutate(
  in_fixed_10k_table = sample_id %in% names(original)[-1],
  excluded_suspect_library = sample_id %in% flagged$sample_id,
  manuscript_habitat = habitat %in% habitats,
  included_manuscript = in_fixed_10k_table & !excluded_suspect_library & manuscript_habitat,
  exclusion_reason = case_when(
    excluded_suspect_library ~ "Suspected contamination or sample-identity problem; cause unresolved",
    !in_fixed_10k_table ~ "Below established 10K target-read threshold",
    !manuscript_habitat ~ "Planktonic community outside manuscript scope",
    TRUE ~ "Included"))
stopifnot(identical(audit$in_fixed_10k_table, audit$included_10k))
write_csv(audit, file.path(output_dir, "sample_inclusion_audit.csv"))
write_csv(audit %>% filter(excluded_suspect_library), file.path(output_dir, "excluded_libraries.csv"))

# 3. Subset the existing draw and remove only ASVs absent from this cohort -------
sample_ids <- names(original)[-1]
sample_ids <- sample_ids[sample_ids %in% audit$sample_id[audit$included_manuscript]]
counts <- as.matrix(original[, sample_ids])
stopifnot(all(colSums(counts) == target_depth), all(counts >= 0),
  all(counts == floor(counts)), !anyNA(counts))
present <- rowSums(counts) > 0
cohort <- original[present, c("asv_id", sample_ids)]
metadata <- audit[match(sample_ids, audit$sample_id), ] %>% mutate(
  habitat_label = recode(habitat, soil = "Soil", hyporheic = "Porewater", sediment = "Sediment"),
  network_group = case_when(habitat == "soil" ~ NA_character_,
    site_stream_order %in% 1:2 ~ "Headwater",
    site_stream_order %in% 3:4 ~ "Intermediate",
    site_stream_order == 5 ~ "Mainstem"))
stopifnot(nrow(metadata) == 76, sum(metadata$habitat == "soil") == 15,
  sum(metadata$habitat == "hyporheic") == 27, sum(metadata$habitat == "sediment") == 34,
  identical(as.matrix(cohort[, -1]), counts[present, , drop = FALSE]))
write_csv(cohort, file.path(output_dir, "asv_counts_10k.csv"))
write_csv(metadata, file.path(output_dir, "sample_metadata.csv"))
write_csv(metadata %>% count(habitat, name = "samples"), file.path(output_dir, "cohort_summary.csv"))

# 4. Audit all candidate pairs, including those lost to depth and QC filtering ---
sediment <- audit %>% filter(habitat == "sediment") %>%
  select(site_code, sediment_id = sample_id, sediment_10k = in_fixed_10k_table,
    sediment_included = included_manuscript)
porewater <- audit %>% filter(habitat == "hyporheic") %>%
  select(site_code, porewater_id = sample_id, porewater_10k = in_fixed_10k_table,
    porewater_flagged = excluded_suspect_library, porewater_included = included_manuscript)
stopifnot(!anyDuplicated(sediment$site_code), !anyDuplicated(porewater$site_code))
pairs <- full_join(sediment, porewater, by = "site_code") %>% mutate(
  paired_10k_before_qc = coalesce(sediment_10k & porewater_10k, FALSE),
  paired_manuscript = coalesce(sediment_included & porewater_included, FALSE),
  pairing_status = case_when(
    is.na(sediment_id) ~ "No sediment library",
    is.na(porewater_id) ~ "No porewater library",
    !sediment_10k | !porewater_10k ~ "One or both below 10K",
    porewater_flagged ~ "Porewater library excluded by review",
    TRUE ~ "Paired"))
stopifnot(sum(pairs$paired_10k_before_qc) == 20, sum(pairs$paired_manuscript) == 18)
write_csv(pairs, file.path(output_dir, "sediment_porewater_pair_audit.csv"))
checks$final_md5 <- unname(tools::md5sum(checks$path))
stopifnot(identical(checks$initial_md5, checks$final_md5))
write_csv(checks, file.path(output_dir, "input_checksums.csv"))
message("Prepared 76 samples and ", nrow(cohort), " ASVs; 18 sediment-porewater pairs.")
