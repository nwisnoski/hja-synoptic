# Within-bin phylogenetic signal in observed habitat affiliation.
# This is a habitat-based diagnostic, not a test of all environmental niches.

library(here)
library(ape)
library(vegan)
library(dplyr)
library(readr)

# 1. Settings ---------------------------------------------------------------
input_dir <- here("results", "icamp_2016", "catchment_castor")
output_dir <- file.path(input_dir, "phylogenetic_signal")
tree_file <- here("results", "phylogeny_2016", "asv_tree_screened_fasttree.nwk")
metadata_file <- here("data", "derived", "microbial_diversity_2016", "sample_metadata.csv")
habitats <- c("planktonic", "hyporheic", "sediment", "soil")
min_occurrences <- 3L
min_asvs <- 6L
permutations <- 999L
random_seed <- 20161009L
correlation_method <- "pearson"
r_cutoff <- 0.1
p_cutoff <- 0.05
profile_methods <- c("habitat_balanced", "observed_sampling")

# 2. Reuse the completed community and bins, without changing the result ------
result_file <- file.path(input_dir, "icamp_result.rds")
input_files <- c(result_file, tree_file, metadata_file)
input_hashes <- tools::md5sum(input_files)
result <- readRDS(result_file)
community <- result$detail$comm
bins <- result$detail$taxabin$sp.bin
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
metadata <- metadata[match(rownames(community), metadata$sample_id), ]
stopifnot(identical(rownames(community), metadata$sample_id))
stopifnot(all(metadata$habitat %in% habitats), all(rowSums(community) == 10000))
stopifnot(setequal(rownames(bins), colnames(community)))
bins <- bins[colnames(community), , drop = FALSE]
tree <- keep.tip(read.tree(tree_file), colnames(community))

relative_abundance <- community / rowSums(community)
asv_abundance <- colMeans(relative_abundance)
occurrences <- colSums(community > 0)
eligible <- occurrences >= min_occurrences
habitat <- factor(metadata$habitat, levels = habitats)
habitat_matrix <- model.matrix(~ habitat - 1)
colnames(habitat_matrix) <- habitats
habitat_n <- colSums(habitat_matrix)
bin_ids <- sort(unique(bins$bin.id.new))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Estimate ASV affiliation as abundance-weighted means of habitat indicators
# This is the niche.value calculation in iCAMP::dniche. Dividing each sample
# by its habitat's sample count gives each habitat equal sampling weight.
# Profiles describe realized associations, not measured physiological traits.
profiles <- list()
for (profile_method in profile_methods) {
  profile_community <- relative_abundance
  if (profile_method == "habitat_balanced") {
    profile_community <- relative_abundance / habitat_n[as.character(habitat)]
  }
  profile <- crossprod(profile_community, habitat_matrix)
  profile <- profile / colSums(profile_community)
  stopifnot(max(abs(rowSums(profile) - 1)) < 1e-10)
  profiles[[profile_method]] <- profile
}

# 4. Test each existing bin with the Mantel test used by iCAMP::ps.bin ---------
# Use only bin-sized tree distances: no full ASV-by-ASV matrix is needed.
# Permutations shuffle ASV identities, not sample or habitat labels.
test_rows <- list()
row_index <- 1L
started <- Sys.time()
for (bin_id in bin_ids) {
  bin_asvs <- rownames(bins)[bins$bin.id.new == bin_id]
  test_asvs <- bin_asvs[eligible[bin_asvs]]
  bin_abundance <- sum(asv_abundance[bin_asvs])
  retained_abundance <- sum(asv_abundance[test_asvs])
  if (length(test_asvs) >= min_asvs) {
    phylo <- cophenetic.phylo(keep.tip(tree, test_asvs))
    phylo <- as.dist(phylo[test_asvs, test_asvs])
  }
  for (profile_method in profile_methods) {
    seed <- random_seed + bin_id
    r <- p <- NA_real_
    status <- "too_few_asvs"
    if (length(test_asvs) >= min_asvs) {
      niche <- dist(profiles[[profile_method]][test_asvs, , drop = FALSE])
      status <- "constant_distance"
      if (sd(as.vector(phylo)) > 0 && sd(as.vector(niche)) > 0) {
        set.seed(seed)
        fit <- mantel(phylo, niche, method = correlation_method,
          permutations = permutations, parallel = 1)
        r <- unname(fit$statistic)
        p <- fit$signif
        status <- "tested"
      }
    }
    test_rows[[row_index]] <- tibble(profile_method, bin_id,
      n_asvs = length(bin_asvs), n_tested_asvs = length(test_asvs),
      bin_abundance, retained_abundance, seed, permutations, status, r, p)
    row_index <- row_index + 1L
  }
  if (bin_id %% 100L == 0L) message("Checked ", bin_id, "/", length(bin_ids), " bins.")
}
tests <- bind_rows(test_rows) %>%
  group_by(profile_method) %>%
  mutate(p_bh = p.adjust(p, method = "BH"),
    positive_nominal = coalesce(r >= r_cutoff & p <= p_cutoff, FALSE),
    positive_bh = coalesce(r >= r_cutoff & p_bh <= p_cutoff, FALSE)) %>%
  ungroup()

# 5. Report coverage and signal without renormalizing away excluded reads ----
signal_summary <- tests %>%
  group_by(profile_method) %>%
  summarize(bins_total = n(), bins_tested = sum(status == "tested"),
    asvs_tested = sum(n_tested_asvs[status == "tested"]),
    tested_abundance = sum(retained_abundance[status == "tested"]),
    positive_bins_nominal = sum(positive_nominal),
    positive_bins_bh = sum(positive_bh),
    positive_abundance_nominal = sum(retained_abundance[positive_nominal]),
    positive_abundance_bh = sum(retained_abundance[positive_bh]),
    weighted_mean_r = weighted.mean(r[status == "tested"],
      retained_abundance[status == "tested"]), .groups = "drop")
profile_rows <- list()
for (profile_method in profile_methods) {
  profile_rows[[profile_method]] <- data.frame(
    asv_id = colnames(community), bin_id = bins$bin.id.new,
    occurrences, mean_relative_abundance = asv_abundance, eligible,
    profile_method, profiles[[profile_method]], row.names = NULL)
}
write_csv(tests, file.path(output_dir, "bin_habitat_signal.csv"))
write_csv(signal_summary, file.path(output_dir, "habitat_signal_summary.csv"))
write_csv(bind_rows(profile_rows), file.path(output_dir, "asv_habitat_affiliation.csv"))
write_csv(tibble(sample_id = metadata$sample_id, habitat = metadata$habitat,
  habitat_sample_count = unname(habitat_n[as.character(habitat)])),
  file.path(output_dir, "sample_habitat_audit.csv"))
write_csv(tibble(input = c("icamp_result.rds", "asv_tree_screened_fasttree.nwk",
  "sample_metadata.csv"), md5 = unname(input_hashes)),
  file.path(output_dir, "input_checksums.csv"))
write_csv(tibble(min_occurrences, min_asvs, permutations, random_seed,
  correlation_method, r_cutoff, p_cutoff,
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))),
  file.path(output_dir, "signal_settings.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
stopifnot(identical(input_hashes, tools::md5sum(input_files)))
print(signal_summary, width = Inf)
