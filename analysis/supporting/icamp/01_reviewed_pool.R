# Regional-pool sensitivity: reviewed soil, porewater, and sediment communities.
# A fresh pooled null model; sediment and porewater figures use only pairs within
# their own habitat. Reopening a completed run regenerates exports and figures.

library(ape)
library(phangorn)
library(iCAMP)
library(here)
library(castor)
library(bigmemory)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1. Settings ------------------------------------------------------------------
data_dir <- here("data", "derived", "soil_porewater_sediment_2016")
taxonomy_file <- here("data", "derived", "microbial_diversity_2016", "asv_metadata_screened.csv")
tree_file <- here("results", "phylogeny_2016", "asv_tree_screened_fasttree.nwk")
output_dir <- here("results", "icamp_2016", "reviewed_pool_castor")
pd_dir <- file.path(output_dir, "phylogenetic_distance")
figures_dir <- here("figures", "supplement", "icamp", "reviewed_pool")
panels_dir <- file.path(figures_dir, "panels")
diagnostics_dir <- file.path(figures_dir, "diagnostics")
rand_time <- 1000L
random_seed <- 2016L
nworker <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
memory_gb <- 150
distance_block_size <- 256L
distance_tolerance <- 1e-8
required_icamp_version <- "1.9.1"
# This workflow always starts fresh randomizations. Do not load historical,
# sediment-only, or pilot checkpoints. Same-run recovery requires a separate audit.
temp_saved_file <- NULL
group_levels <- c("Headwater", "Intermediate", "Mainstem")
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
n_balanced_draws <- 1000L
segments_per_group <- 5L
base_seed <- 20261012L
centerline_stride <- 10L
map_width <- 6.8
map_height <- 5.6
comparison_width <- 5.4
comparison_height <- 4.5
process_design <- tribble(
  ~process, ~code, ~label,
  "Heterogeneous.Selection", "HeS", "Heterogeneous selection",
  "Homogeneous.Selection", "HoS", "Homogeneous selection",
  "Dispersal.Limitation", "DL", "Dispersal limitation",
  "Homogenizing.Dispersal", "HD", "Homogenizing dispersal",
  "Drift.and.Others", "DR", "Drift and other processes"
)
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
process_colors <- c(HeS = "#D55E00", HoS = "#009E73", DL = "#0072B2",
  HD = "#CC79A7", DR = "#999999")
plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 10),
    plot.margin = margin(8, 20, 8, 8)
  )

# 2. Validate the reviewed fixed-10K cohort and its master-tree coverage ---------
preparation_start <- proc.time()[["elapsed"]]
input_paths <- c(
  file.path(data_dir, "asv_counts_10k.csv"),
  file.path(data_dir, "sample_metadata.csv"), taxonomy_file, tree_file
)
stopifnot(all(file.exists(input_paths)))
input_checksums <- data.frame(
  input = sub(paste0(here(), "/"), "", input_paths, fixed = TRUE),
  md5 = unname(tools::md5sum(input_paths))
)
metadata <- read_csv(input_paths[2], col_types = cols(site_code = col_character()))
counts_table <- read_csv(input_paths[1], show_col_types = FALSE)
taxonomy <- read_csv(taxonomy_file, show_col_types = FALSE)
stopifnot(!anyDuplicated(metadata$sample_id), !anyDuplicated(counts_table$asv_id),
  !anyDuplicated(taxonomy$asv_id), identical(names(counts_table)[-1], metadata$sample_id))
community <- t(as.matrix(counts_table[, -1]))
colnames(community) <- counts_table$asv_id
stopifnot(all(is.finite(community)), all(community >= 0),
  all(community == round(community)), all(rowSums(community) == 10000),
  all(colSums(community) > 0), ncol(community) == 33974L,
  nrow(community) == 76L, sum(metadata$habitat == "soil") == 15L,
  sum(metadata$habitat == "sediment") == 34L, sum(metadata$habitat == "hyporheic") == 27L,
  all(metadata$habitat %in% c("soil", "sediment", "hyporheic")),
  !any(metadata$excluded_suspect_library), all(metadata$included_manuscript))
storage.mode(community) <- "integer"
rm(counts_table)
stream_metadata <- metadata %>% filter(habitat != "soil")
stopifnot(!anyNA(stream_metadata[, c("site_code", "site_segment_number",
  "site_stream_order", "site_utm_x_m", "site_utm_y_m")]),
  all(stream_metadata$site_stream_order %in% c(headwater_orders, intermediate_orders, mainstem_order)),
  !anyDuplicated(paste(stream_metadata$habitat, stream_metadata$site_code)))
# Missing soil GPS/identity flags stay in the manifest; soils are a regional pool.
tree <- read.tree(tree_file)
stopifnot(all(colnames(community) %in% tree$tip.label),
  all(colnames(community) %in% taxonomy$asv_id))
tree <- midpoint(keep.tip(tree, colnames(community)))
tree$node.label <- NULL
stopifnot(all(is.finite(tree$edge.length)), all(tree$edge.length >= 0), is.rooted(tree))
message("Validated ", nrow(community), " samples and ", ncol(community), " ASVs in ",
  round(proc.time()[["elapsed"]] - preparation_start, 2), " seconds. Full distance matrix: ",
  round(8 * Ntip(tree)^2 / 1024^3, 2), " GiB before temporary copies.")

# Record the analytical identity before any expensive work. Completed runs may
# be plotted with newer local packages; unfinished runs require matching versions.
analysis_identity <- list(inputs = input_checksums, samples = rownames(community),
  asvs = colnames(community), tree_tips = tree$tip.label, randomizations = rand_time,
  seed = random_seed, phylo_metric = "bMPD", sig_index = "Confidence",
  phylo_rand_scale = "within.bin", taxa_rand_scale = "across.all", taxo_metric = "bray",
  bin_size = 24L, ds = 0.2, rooting = "midpoint after pruning",
  pool_weighting = "mean relative abundance across all 76 samples",
  icamp_version = required_icamp_version)
inference_packages <- c("ape", "phangorn", "iCAMP", "castor", "bigmemory", "vegan", "parallel")
inference_versions <- c(R = as.character(getRversion()),
  setNames(vapply(inference_packages, function(x) as.character(packageVersion(x)), ""),
    inference_packages))
result_file <- file.path(output_dir, "icamp_result.rds")
manifest_file <- file.path(output_dir, "run_manifest.rds")
if (file.exists(manifest_file)) {
  saved_manifest <- readRDS(manifest_file)
  stopifnot(identical(saved_manifest$analysis, analysis_identity))
  if (!file.exists(result_file)) {
    stopifnot(identical(saved_manifest$versions, inference_versions))
  }
} else if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE))) {
  stop("Unidentified outputs exist. Preserve them and choose a fresh output directory.")
}

# 3. Run a fresh pooled null model only when no completed result exists -----------
if (!file.exists(result_file)) {
  stopifnot(as.character(packageVersion("iCAMP")) == required_icamp_version,
    is.finite(nworker), nworker >= 1L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pd_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  pd_dir <- normalizePath(pd_dir)
  if (!file.exists(manifest_file)) {
    saveRDS(list(analysis = analysis_identity, versions = inference_versions), manifest_file)
  }
  write_csv(input_checksums, file.path(output_dir, "input_checksums.csv"))
  write_csv(metadata, file.path(output_dir, "pooled_samples.csv"), na = "")
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "inference_session_info.txt"))
  distance_inputs <- data.frame(
    input = c("asv_counts_10k.csv", basename(tree_file)),
    md5 = unname(tools::md5sum(c(file.path(data_dir, "asv_counts_10k.csv"), tree_file))),
    method = "castor branch-length tip distances; midpoint after pruning v1",
    versions = paste(inference_versions[c("ape", "phangorn", "castor", "bigmemory")], collapse = ";")
  )
  completion_file <- file.path(pd_dir, "distance_complete.csv")
  descriptor_file <- file.path(pd_dir, "pd.desc")
  taxon_file <- file.path(pd_dir, "pd.taxon.name.csv")
  distance_start <- Sys.time()

  if (file.exists(completion_file)) {
    saved_inputs <- read.csv(completion_file, stringsAsFactors = FALSE)
    saved_tips <- read.csv(taxon_file, stringsAsFactors = FALSE)$taxon
    if (!identical(saved_inputs, distance_inputs) ||
        !identical(saved_tips, tree$tip.label)) {
      stop("Distance inputs or tip order changed. Use a new distance directory.")
    }
    pd <- attach.big.matrix(descriptor_file)
    message("Reusing completed tip distances; checking the saved matrix.")
  } else {
    if (file.exists(descriptor_file) || file.exists(file.path(pd_dir, "pd.bin"))) {
      stop("An incomplete distance matrix exists in: ", pd_dir,
           ". Preserve it and choose a fresh directory; do not reuse it.")
    }
    message("Calculating branch-length distances for ", Ntip(tree), " ASVs.")
    distances <- get_all_pairwise_distances(
      tree, only_clades = seq_len(Ntip(tree)), as_edge_counts = FALSE
    )
    stopifnot(identical(dim(distances), c(Ntip(tree), Ntip(tree))))
    pd <- filebacked.big.matrix(
      nrow = Ntip(tree), ncol = Ntip(tree), type = "double",
      backingpath = pd_dir, backingfile = "pd.bin", descriptorfile = "pd.desc"
    )
    for (first in seq(1, Ntip(tree), by = distance_block_size)) {
      columns <- first:min(first + distance_block_size - 1, Ntip(tree))
      pd[, columns] <- distances[, columns, drop = FALSE]
    }
    flush(pd)
    write.csv(data.frame(taxon = tree$tip.label), taxon_file, row.names = TRUE)
    rm(distances)
    gc()
  }

  # Check the entire saved matrix in small blocks, not another full-size copy.
  stopifnot(all(dim(pd) == c(Ntip(tree), Ntip(tree))))
  for (first in seq(1, Ntip(tree), by = distance_block_size)) {
    columns <- first:min(first + distance_block_size - 1, Ntip(tree))
    block <- pd[, columns, drop = FALSE]
    stopifnot(all(is.finite(block)), all(block >= 0))
    stopifnot(max(abs(block - t(pd[columns, , drop = FALSE]))) < distance_tolerance)
    stopifnot(all(block[cbind(columns, seq_along(columns))] == 0))
  }
  # Independently compare 64 spread-out tips with ape after pruning.
  check_indices <- unique(round(seq(1, Ntip(tree), length.out = min(64, Ntip(tree)))))
  check_tips <- tree$tip.label[check_indices]
  reference <- cophenetic.phylo(keep.tip(tree, check_tips))[check_tips, check_tips]
  max_error <- max(abs(pd[check_indices, check_indices] - reference))
  stopifnot(max_error < distance_tolerance)
  write.csv(
    data.frame(asvs = Ntip(tree), reused = file.exists(completion_file),
               elapsed_seconds = as.numeric(difftime(Sys.time(), distance_start, units = "secs")),
               max_ape_error = max_error, castor_version = as.character(packageVersion("castor"))),
    file.path(output_dir, "distance_checks.csv"), row.names = FALSE
  )
  # This record is written LAST: the descriptor alone never means completion.
  write.csv(distance_inputs, completion_file, row.names = FALSE)
  rm(pd, block, reference)
  gc()

  write.csv(
    data.frame(
      setting = c(
        "samples", "asvs", "reads_per_sample", "randomizations", "random_seed",
        "workers", "phylogenetic_metric", "significance_index",
        "phylogenetic_randomization", "taxonomic_randomization",
        "minimum_bin_size", "phylogenetic_signal_distance", "rooting", "pool_weighting", "iCAMP_version"
      ),
      value = c(
        nrow(community), ncol(community), 10000, rand_time, random_seed,
        nworker, "bMPD", "Confidence", "within.bin", "across.all",
        24, 0.2, "midpoint root after pruning to reviewed-cohort ASVs",
        "sample-weighted mean relative abundance across 76 samples", required_icamp_version
      )
    ),
    file.path(output_dir, "icamp_settings.csv"),
    row.names = FALSE
  )

  set.seed(random_seed)
  icamp_result <- icamp.big(
    comm = community,
    tree = tree,
    pd.desc = "pd.desc",
    pd.spname = tree$tip.label,
    pd.wd = pd_dir,
    rand = rand_time,
    prefix = "hja2016_reviewed_pool",
    ds = 0.2,
    phylo.rand.scale = "within.bin",
    taxa.rand.scale = "across.all",
    taxo.metric = "bray",
    phylo.metric = "bMPD",
    sig.index = "Confidence",
    bin.size.limit = 24,
    nworker = nworker,
    memory.G = memory_gb,
    rtree.save = FALSE,
    detail.save = TRUE,
    qp.save = FALSE,
    detail.null = FALSE,
    ignore.zero = TRUE,
    output.wd = output_dir,
    unit.sum = rowSums(community),
    # Make the default regional weighting explicit; habitats are not given equal mass.
    meta.ab = colMeans(community / rowSums(community)),
    temp.save = TRUE,
    temp.saved.file = temp_saved_file
  )
  saveRDS(icamp_result, file.path(output_dir, "icamp_result.rds"))

}

# 4. Export and validate the completed pooled assignments -----------------------
icamp_result <- readRDS(result_file)
pairwise <- as_tibble(icamp_result$CbMPDiCBraya)
stopifnot(nrow(pairwise) == choose(nrow(community), 2),
  all(pairwise$sample1 %in% metadata$sample_id), all(pairwise$sample2 %in% metadata$sample_id),
  all(pairwise$sample1 != pairwise$sample2),
  !anyDuplicated(paste(pmin(pairwise$sample1, pairwise$sample2),
    pmax(pairwise$sample1, pairwise$sample2), sep = "|")))
fractions <- as.matrix(pairwise[, process_design$process])
stopifnot(all(is.finite(fractions)), all(fractions >= 0 & fractions <= 1),
  max(abs(rowSums(fractions) - 1)) < 1e-10)
icamp_summary <- icamp.bins(icamp_result$detail, silent = TRUE)
write.csv(
  icamp_result$CbMPDiCBraya,
  file.path(output_dir, "pairwise_process_importance.csv"),
  row.names = FALSE
)
write.csv(
  icamp_summary$Pt,
  file.path(output_dir, "overall_process_importance.csv"),
  row.names = FALSE
)
write.csv(
  icamp_summary$Ptk,
  file.path(output_dir, "process_importance_by_bin.csv"),
  row.names = FALSE
)
write.csv(
  icamp_summary$Binwt,
  file.path(output_dir, "bin_relative_abundance.csv"),
  row.names = FALSE
)

taxa_bins <- data.frame(
  asv_id = rownames(icamp_result$detail$taxabin$sp.bin),
  icamp_result$detail$taxabin$sp.bin,
  row.names = NULL,
  check.names = FALSE
)
stopifnot(setequal(taxa_bins$asv_id, taxonomy$asv_id[taxonomy$asv_id %in% colnames(community)]))
taxonomy <- taxonomy[match(taxa_bins$asv_id, taxonomy$asv_id), ]
taxa_bins <- cbind(taxa_bins, taxonomy[setdiff(names(taxonomy), "asv_id")])
write.csv(
  taxa_bins,
  file.path(output_dir, "asv_bin_taxonomy.csv"),
  row.names = FALSE,
  na = ""
)


# 5. Identical habitat-specific maps and network summaries -----------------------
plot_habitat <- function(pairwise, site, prefix, output_dir, figures_dir, map_lines) {
  panels_dir <- file.path(figures_dir, "panels")
  diagnostics_dir <- file.path(figures_dir, "diagnostics")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(nrow(pairwise) == choose(nrow(site), 2), !anyDuplicated(site$sample_id),
    !anyDuplicated(site$site_code))
  site <- site %>% select(sample_id, site_code, site_stream_order,
    site_drainage_area_ha, site_distance_to_outlet_m, site_utm_x_m, site_utm_y_m,
    site_segment_number)
  site <- site %>% mutate(network_group = case_when(
    site_stream_order %in% headwater_orders ~ "Headwater",
    site_stream_order %in% intermediate_orders ~ "Intermediate",
    site_stream_order == mainstem_order ~ "Mainstem"),
    group_index = match(network_group, group_levels))
  stopifnot(!anyNA(site$group_index))
  sizes <- site %>% group_by(network_group) %>%
    summarize(sites = n(), segments = n_distinct(site_segment_number), .groups = "drop")
  stopifnot(nrow(sizes) == length(group_levels), all(sizes$segments >= segments_per_group))
  segment_groups <- site %>% distinct(site_segment_number, network_group)
  stopifnot(!anyDuplicated(segment_groups$site_segment_number))
  # 5a. Describe sites relative to other sampled stream segments.
  # Average over sites within each partner segment, then weight partner segments
  # equally. Excluding the focal segment prevents nearby replicates from defining
  # a site's mapped value. These remain pairwise summaries, not local process rates.
  first <- match(pairwise$sample1, site$sample_id)
  second <- match(pairwise$sample2, site$sample_id)
  pairwise <- pairwise %>%
    mutate(
      segment1 = site$site_segment_number[first],
      segment2 = site$site_segment_number[second],
      group1 = site$group_index[first],
      group2 = site$group_index[second],
      same_segment = segment1 == segment2,
      segment_low = pmin(segment1, segment2),
      segment_high = pmax(segment1, segment2),
      group_low = pmin(group1, group2),
      group_high = pmax(group1, group2)
    )
  long_pairs <- pairwise %>% pivot_longer(all_of(process_design$process),
    names_to = "process", values_to = "importance")
  incident <- bind_rows(
    long_pairs %>% transmute(sample_id = sample1, partner_id = sample2,
      partner_segment = segment2, same_segment, process, importance),
    long_pairs %>% transmute(sample_id = sample2, partner_id = sample1,
      partner_segment = segment1, same_segment, process, importance)
  )
  site_profiles <- incident %>% filter(!same_segment) %>%
    group_by(sample_id, partner_segment, process) %>%
    summarize(segment_mean = mean(importance), partner_sites = n(), .groups = "drop") %>%
    group_by(sample_id, process) %>%
    summarize(segment_balanced_mean = mean(segment_mean),
      partner_segments = n(), partner_sites = sum(partner_sites), .groups = "drop") %>%
    left_join(incident %>% group_by(sample_id, process) %>%
      summarize(all_partner_mean = mean(importance), .groups = "drop"),
      by = c("sample_id", "process")) %>%
    left_join(site, by = "sample_id")
  stopifnot(nrow(site_profiles) == nrow(site) * nrow(process_design),
    all(site_profiles$partner_segments == nrow(segment_groups) - 1L))
  site_sums <- site_profiles %>% group_by(sample_id) %>%
    summarize(total = sum(segment_balanced_mean), .groups = "drop")
  stopifnot(max(abs(site_sums$total - 1)) < 1e-10)
  write_csv(site_profiles, file.path(output_dir, "site_process_importance.csv"))

  # 5b. Give distinct segment pairs equal weight in each network comparison.
  # Within-segment pairs are retained in the audit but excluded from comparisons.
  # Six cells distinguish within-group turnover from turnover between positions.
  segment_pairs <- long_pairs %>% filter(!same_segment) %>%
    group_by(group_low, group_high, segment_low, segment_high, process) %>%
    summarize(importance = mean(importance), sample_pairs = n(), .groups = "drop")
  group_profiles <- segment_pairs %>% group_by(group_low, group_high, process) %>%
    summarize(segment_balanced_mean = mean(importance),
      segment_pairs = n(), sample_pairs = sum(sample_pairs), .groups = "drop") %>%
    left_join(long_pairs %>% filter(!same_segment) %>%
      group_by(group_low, group_high, process) %>%
      summarize(sample_pair_mean = mean(importance), .groups = "drop"),
      by = c("group_low", "group_high", "process")) %>%
    mutate(group1 = group_levels[group_low], group2 = group_levels[group_high])
  write_csv(group_profiles, file.path(output_dir, "network_pair_process_importance.csv"))
  pair_audit <- pairwise %>%
    transmute(sample1, sample2, site1 = site$site_code[first], site2 = site$site_code[second],
      segment1, segment2, group1 = group_levels[group1], group2 = group_levels[group2],
      same_segment, included_network_comparison = !same_segment)
  write_csv(pair_audit, file.path(output_dir, "pair_structure_audit.csv"))

  # 5c. Select one site in each of five segments per group.
  # The same 15 sites supply every process and every group comparison in a draw.
  # This samples the completed assignments; it does not rerun the iCAMP null model.
  draw_design <- tibble(run = seq_len(n_balanced_draws), seed = base_seed + run)
  draw_rows <- vector("list", n_balanced_draws)
  selection_rows <- vector("list", n_balanced_draws)
  start_time <- proc.time()[["elapsed"]]
  for (run in draw_design$run) {
    set.seed(draw_design$seed[run])
    selected <- integer(0)
    for (group in seq_along(group_levels)) {
      candidates <- which(site$group_index == group)
      segments <- unique(site$site_segment_number[candidates])
      selected_segments <- segments[sample.int(length(segments), segments_per_group)]
      for (segment in selected_segments) {
        indices <- candidates[site$site_segment_number[candidates] == segment]
        selected <- c(selected, indices[sample.int(length(indices), 1)])
      }
    }
    selected_ids <- site$sample_id[selected]
    selection_rows[[run]] <- site[selected, ] %>%
      transmute(run, seed = draw_design$seed[run], sample_id, site_code,
        network_group, segment = site_segment_number)
    draw_rows[[run]] <- long_pairs %>%
      filter(sample1 %in% selected_ids, sample2 %in% selected_ids) %>%
      group_by(group_low, group_high, process) %>%
      summarize(importance = mean(importance), sample_pairs = n(), .groups = "drop") %>%
      mutate(run, seed = draw_design$seed[run])
    if (run == 1L) {
      message("Representative balanced draw: ",
        round(proc.time()[["elapsed"]] - start_time, 3), " s")
    }
  }
  draws <- bind_rows(draw_rows)
  selections <- bind_rows(selection_rows)
  selection_audit <- selections %>% group_by(run, network_group) %>%
    summarize(sites = n(), segments = n_distinct(segment), .groups = "drop")
  stopifnot(nrow(draws) == n_balanced_draws * 6L * nrow(process_design),
    all(selection_audit$sites == segments_per_group),
    all(selection_audit$segments == segments_per_group),
    all(draws$sample_pairs[draws$group_low == draws$group_high] == choose(segments_per_group, 2)),
    all(draws$sample_pairs[draws$group_low != draws$group_high] == segments_per_group^2))
  balanced_summary <- draws %>% group_by(group_low, group_high, process) %>%
    summarize(draws = n(), median = median(importance),
      lower = quantile(importance, 0.025), upper = quantile(importance, 0.975),
      .groups = "drop") %>%
    mutate(group1 = group_levels[group_low], group2 = group_levels[group_high])
  write_csv(draws, file.path(output_dir, "balanced_segment_draws.csv"))
  write_csv(selections, file.path(output_dir, "balanced_site_selections.csv"))
  # Keep the draw-level contrasts so the sensitivity of differences is inspectable.
  contrasts <- tibble(group_a = c(1L, 1L, 2L), group_b = c(2L, 3L, 3L))
  contrast_rows <- vector("list", nrow(contrasts))
  within_draws <- draws %>% filter(group_low == group_high)
  for (i in seq_len(nrow(contrasts))) {
    a <- within_draws %>% filter(group_low == contrasts$group_a[i])
    b <- within_draws %>% filter(group_low == contrasts$group_b[i])
    stopifnot(identical(a$run, b$run), identical(a$process, b$process))
    contrast_rows[[i]] <- tibble(run = a$run, seed = a$seed, process = a$process,
      contrast = paste(group_levels[contrasts$group_a[i]], group_levels[contrasts$group_b[i]], sep = " - "),
      difference = a$importance - b$importance)
  }
  write_csv(bind_rows(contrast_rows), file.path(output_dir, "balanced_within_group_contrasts.csv"))

  # 5d. Map segment-balanced pairwise importance at sampled locations.
  maps <- vector("list", nrow(process_design))
  comparisons <- vector("list", nrow(process_design))
  size_reference <- sizes[match(group_levels, sizes$network_group), ]
  group_labels <- paste0(group_levels, "\n", size_reference$sites, " sites\n",
    size_reference$segments, " segments")
  for (i in seq_len(nrow(process_design))) {
    current_process <- process_design$process[i]
    mapped <- site_profiles %>% filter(process == current_process)
    map <- ggplot(map_lines, aes(x, y, group = segment)) +
      geom_path(color = "#BBBBBB", linewidth = 0.3) +
      geom_point(data = mapped,
        aes(site_utm_x_m, site_utm_y_m, fill = 100 * segment_balanced_mean),
        inherit.aes = FALSE, shape = 21, size = 3.1, color = "black", stroke = 0.35) +
      # Each process has its own color range to reveal spatial differences;
      # its numerical legend remains in percentages for direct interpretation.
      scale_fill_viridis_c(option = "C",
      breaks = scales::breaks_pretty(n = 3),
      labels = function(x) sprintf("%.1f", x), name = "Mean pairwise importance (%)",
      guide = guide_colorbar(barwidth = grid::unit(3.2, "cm"))) +
      coord_equal() + labs(x = "UTM easting (m)", y = "UTM northing (m)") +
      facet_wrap(~process, labeller = as_labeller(setNames(process_design$label, process_design$process))) +
      plot_theme
    maps[[i]] <- map
    ggsave(file.path(panels_dir, paste0(prefix, "_map_", process_design$code[i], ".pdf")),
      map, width = map_width, height = map_height, bg = "white")

    # Each range is the middle 95% of balanced draws, with their median marked.
    # Hollow diamonds retain the mean using every observed segment pair.
    panel <- balanced_summary %>%
      filter(process == current_process, group_low == group_high) %>%
      mutate(network_group = factor(group1, levels = group_levels))
    full <- group_profiles %>%
      filter(process == current_process, group_low == group_high) %>%
      mutate(network_group = factor(group1, levels = group_levels))
    comparison <- ggplot(panel,
      aes(group_low - 0.06, 100 * median, color = network_group)) +
      geom_linerange(aes(ymin = 100 * lower, ymax = 100 * upper), linewidth = 1) +
      geom_point(aes(shape = "Balanced draws: median"), size = 2.8) +
      geom_point(data = full,
        aes(x = group_low + 0.06, y = 100 * segment_balanced_mean,
          shape = "All segment pairs: mean"),
        fill = "white", size = 2.8, stroke = 0.8) +
      scale_color_manual(values = group_colors, guide = "none") +
      scale_shape_manual(values = c("Balanced draws: median" = 16,
        "All segment pairs: mean" = 23), name = NULL) +
      scale_x_continuous(breaks = seq_along(group_levels), labels = group_labels,
        limits = c(0.5, 3.5)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.04, 0.08))) +
      labs(x = NULL, y = "Within-group pairwise importance (%)") +
      facet_wrap(~process, labeller = as_labeller(setNames(process_design$label, process_design$process))) +
      plot_theme + theme(aspect.ratio = 0.8, axis.text.x = element_text(size = 9))
    comparisons[[i]] <- comparison
    ggsave(file.path(panels_dir, paste0(prefix, "_network_", process_design$code[i], ".pdf")),
      comparison, width = comparison_width, height = comparison_height, bg = "white")
  }
  # The sixth map identifies the network positions used in the comparisons.
  position_map <- ggplot(map_lines, aes(x, y, group = segment)) +
    geom_path(color = "#BBBBBB", linewidth = 0.3) +
    geom_point(data = site %>% mutate(process = "Network position"),
      aes(site_utm_x_m, site_utm_y_m, fill = factor(network_group, levels = group_levels)),
      inherit.aes = FALSE, shape = 21, size = 3.1, color = "black", stroke = 0.35) +
    scale_fill_manual(values = group_colors, name = NULL) +
    coord_equal() + labs(x = "UTM easting (m)", y = "UTM northing (m)") +
    facet_wrap(~process) + plot_theme
  ggsave(file.path(figures_dir, paste0(prefix, "_catchment_maps.pdf")),
    wrap_plots(c(maps, list(position_map)), ncol = 3),
    width = 15.6, height = 10.4, bg = "white")
  ggsave(file.path(diagnostics_dir, paste0(prefix, "_network_comparison.pdf")),
    wrap_plots(comparisons, ncol = 3, guides = "collect") & theme(legend.position = "bottom"),
    width = 15.6, height = 9, bg = "white")

  # A compact composition view keeps all processes on the common 0-100% scale.
  composition <- group_profiles %>% filter(group_low == group_high) %>%
    mutate(network_group = factor(group1, levels = group_levels),
      code = factor(process_design$code[match(process, process_design$process)],
        levels = process_design$code))
  composition_plot <- ggplot(composition,
    aes(network_group, 100 * segment_balanced_mean, fill = code)) +
    geom_col(width = 0.65, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = process_colors,
      labels = setNames(process_design$label, process_design$code), name = NULL) +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Within-group pairwise importance (%)") + plot_theme +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    theme(axis.text.x = element_text(size = 9))
  ggsave(file.path(figures_dir, paste0(prefix, "_network_composition.pdf")),
    composition_plot, width = 7.2, height = 5.5, bg = "white")

  # 5e. Compare within-position and between-position assignments.
  heatmap <- group_profiles %>%
    mutate(process = factor(process, levels = process_design$process),
      group1 = factor(group1, levels = group_levels),
      group2 = factor(group2, levels = rev(group_levels)),
      label = sprintf("%.1f", 100 * segment_balanced_mean))
  pair_plot <- ggplot(heatmap, aes(group1, group2, fill = 100 * segment_balanced_mean)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label, color = segment_balanced_mean < 0.5), size = 3.6) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_viridis_c(option = "C", limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100), name = "Mean pairwise importance (%)") +
    facet_wrap(~process, ncol = 3,
      labeller = as_labeller(setNames(process_design$label, process_design$process))) +
    coord_equal() + labs(x = NULL, y = NULL) + plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(diagnostics_dir, paste0(prefix, "_network_pair_comparison.pdf")),
    pair_plot, width = 13, height = 8.5, bg = "white")

  # 5f. Show every within-position sample pair with equal weight.
  # Include same-segment pairs here, marking them separately. The comparison is
  # all within-position pairs, each using all eligible habitat sites.
  # Individual pair rows share sites and are not independent replicates.
  within_pairs <- long_pairs %>% filter(group1 == group2) %>%
    mutate(network_group = factor(group_levels[group1], levels = group_levels),
      pair_type = if_else(same_segment, "Same segment", "Different segments"))
  within_summary <- within_pairs %>% group_by(network_group, process) %>%
    summarize(sites = n_distinct(c(sample1, sample2)), sample_pairs = n(),
      same_segment_pairs = sum(same_segment),
      mean = mean(importance), median = median(importance),
      lower_quartile = quantile(importance, 0.25),
      upper_quartile = quantile(importance, 0.75),
      minimum = min(importance), maximum = max(importance), .groups = "drop")
  within_sizes <- within_summary %>% distinct(network_group, sites, sample_pairs)
  stopifnot(all(within_sizes$sample_pairs == choose(within_sizes$sites, 2)))
  within_labels <- paste0(group_levels, "\n", within_sizes$sites, " sites\n",
    within_sizes$sample_pairs, " pairs")
  write_csv(within_summary, file.path(output_dir, "all_within_group_pair_summary.csv"))
  all_pair_panels <- vector("list", nrow(process_design))
  for (i in seq_len(nrow(process_design))) {
    current_process <- process_design$process[i]
    panel <- within_pairs %>% filter(process == current_process)
    summary <- within_summary %>% filter(process == current_process)
    pair_distribution <- ggplot(panel,
      aes(network_group, 100 * importance, color = network_group)) +
      geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.6) +
      geom_point(aes(shape = pair_type),
        position = position_jitter(width = 0.12, height = 0, seed = base_seed),
        alpha = 0.5, size = 1.4) +
      geom_point(data = summary, aes(y = 100 * mean, shape = "Group mean"),
        fill = "white", color = "black", size = 2.8, stroke = 0.8) +
      scale_color_manual(values = group_colors, guide = "none") +
      scale_shape_manual(values = c("Different segments" = 16,
        "Same segment" = 4, "Group mean" = 23), name = NULL) +
      guides(shape = guide_legend(override.aes = list(color = "black", alpha = 1))) +
      scale_x_discrete(labels = within_labels) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.06))) +
      labs(x = NULL, y = "Pairwise importance (%)") +
      facet_wrap(~process,
        labeller = as_labeller(setNames(process_design$label, process_design$process))) +
      plot_theme + theme(aspect.ratio = 0.8, axis.text.x = element_text(size = 9))
    all_pair_panels[[i]] <- pair_distribution
    ggsave(file.path(panels_dir, paste0(prefix, "_all_within_pairs_",
      process_design$code[i], ".pdf")), pair_distribution,
      width = comparison_width, height = comparison_height, bg = "white")
  }
  ggsave(file.path(diagnostics_dir, paste0(prefix, "_all_within_pairs.pdf")),
    wrap_plots(all_pair_panels, ncol = 3, guides = "collect") & theme(legend.position = "bottom"),
    width = 15.6, height = 9, bg = "white")
  all_pair_composition <- within_summary %>%
    mutate(code = factor(process_design$code[match(process, process_design$process)],
      levels = process_design$code))
  all_pair_composition_plot <- ggplot(all_pair_composition,
    aes(network_group, 100 * mean, fill = code)) +
    geom_col(width = 0.65, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = process_colors,
      labels = setNames(process_design$label, process_design$code), name = NULL) +
    scale_x_discrete(labels = within_labels) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Mean within-group pairwise importance (%)") + plot_theme +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    theme(axis.text.x = element_text(size = 9))
  ggsave(file.path(diagnostics_dir, paste0(prefix, "_all_within_pair_composition.pdf")),
    all_pair_composition_plot, width = 7.2, height = 5.5, bg = "white")

  invisible(list(within_summary = within_summary, site_profiles = site_profiles))
}

# 6. Draw each habitat using only its own completed pairwise assignments ----------
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panels_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
geometry <- read_csv(here("data", "StreamCenterline.csv"), skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
stopifnot(nrow(problems(geometry)) == 0L,
  all(is.finite(geometry$x)), all(is.finite(geometry$y)))
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% centerline_stride == 1L | row_number() == n()) %>% ungroup()
habitat_design <- tribble(~habitat, ~label, ~filename,
  "sediment", "Sediment", "sediment", "hyporheic", "Porewater", "porewater")
for (h in seq_len(nrow(habitat_design))) {
  site <- metadata %>% filter(habitat == habitat_design$habitat[h])
  habitat_pairs <- pairwise %>% filter(sample1 %in% site$sample_id, sample2 %in% site$sample_id)
  plot_habitat(habitat_pairs, site,
    paste0("2016_reviewed_pool_icamp_", habitat_design$filename[h]),
    file.path(output_dir, "network_figures", habitat_design$filename[h]), figures_dir, map_lines)
}

# 7. Match habitat coverage by using the same 18 sites and unordered site pairs ---
# This is a coverage sensitivity, using the pooled assignments without refitting.
# Three intermediate sites remain, so no five-segment resampling is attempted here.
shared_sites <- intersect(metadata$site_code[metadata$habitat == "sediment"],
  metadata$site_code[metadata$habitat == "hyporheic"])
stopifnot(length(shared_sites) == 18L)
matched_samples <- metadata %>% filter(habitat %in% habitat_design$habitat, site_code %in% shared_sites)
matched_geometry <- matched_samples %>% group_by(site_code) %>%
  summarize(orders = n_distinct(site_stream_order), segments = n_distinct(site_segment_number),
    eastings = n_distinct(site_utm_x_m), northings = n_distinct(site_utm_y_m), .groups = "drop")
stopifnot(all(matched_geometry$orders == 1L), all(matched_geometry$segments == 1L),
  all(matched_geometry$eastings == 1L), all(matched_geometry$northings == 1L))
write_csv(matched_samples %>% select(sample_id, site_code, habitat, network_group,
  site_segment_number, site_stream_order), file.path(output_dir, "matched_stream_samples.csv"))
matched_rows <- vector("list", nrow(habitat_design))
for (h in seq_len(nrow(habitat_design))) {
  site <- matched_samples %>% filter(habitat == habitat_design$habitat[h])
  pairs <- pairwise %>% filter(sample1 %in% site$sample_id, sample2 %in% site$sample_id)
  first <- match(pairs$sample1, site$sample_id)
  second <- match(pairs$sample2, site$sample_id)
  stopifnot(nrow(pairs) == choose(length(shared_sites), 2))
  matched_rows[[h]] <- pairs %>% mutate(
    site_low = pmin(site$site_code[first], site$site_code[second]),
    site_high = pmax(site$site_code[first], site$site_code[second]),
    group1 = site$network_group[first], group2 = site$network_group[second],
    habitat = habitat_design$label[h]) %>%
    filter(group1 == group2) %>%
    pivot_longer(all_of(process_design$process), names_to = "process", values_to = "importance") %>%
    select(site_low, site_high, network_group = group1, habitat, process, importance)
}
matched_pairs <- bind_rows(matched_rows) %>%
  pivot_wider(names_from = habitat, values_from = importance) %>%
  mutate(porewater_minus_sediment = Porewater - Sediment)
stopifnot(!anyNA(matched_pairs), nrow(matched_pairs) == (choose(11, 2) + choose(3, 2) + choose(4, 2)) * 5)
write_csv(matched_pairs, file.path(output_dir, "matched_within_group_pair_contrasts.csv"))
matched_summary <- matched_pairs %>% group_by(network_group, process) %>%
  summarize(pairs = n(), sediment_mean = mean(Sediment), porewater_mean = mean(Porewater),
    mean_difference = mean(porewater_minus_sediment), .groups = "drop")
write_csv(matched_summary, file.path(output_dir, "matched_within_group_summary.csv"))
matched_plot_data <- matched_summary %>%
  pivot_longer(c(sediment_mean, porewater_mean), names_to = "habitat", values_to = "importance") %>%
  mutate(network_group = factor(network_group, levels = group_levels),
    habitat = if_else(habitat == "sediment_mean", "Sediment", "Porewater"))
matched_panels <- vector("list", nrow(process_design))
for (i in seq_len(nrow(process_design))) {
  panel <- matched_plot_data %>% filter(process == process_design$process[i])
  matched_plot <- ggplot(panel, aes(network_group, 100 * importance, color = habitat, group = habitat)) +
    geom_line(linewidth = 0.6) + geom_point(aes(shape = habitat), size = 2.8) +
    scale_color_manual(values = c(Sediment = "#0072B2", Porewater = "#D55E00"), name = NULL) +
    scale_shape_manual(values = c(Sediment = 16, Porewater = 17), name = NULL) +
    scale_x_discrete(labels = c("Headwater\n11 sites\n55 pairs", "Intermediate\n3 sites\n3 pairs", "Mainstem\n4 sites\n6 pairs")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.06))) +
    facet_wrap(~process, labeller = as_labeller(setNames(process_design$label, process_design$process))) +
    labs(x = NULL, y = "Mean within-group pairwise importance (%)") + plot_theme +
    theme(aspect.ratio = 0.8, axis.text.x = element_text(size = 9))
  matched_panels[[i]] <- matched_plot
  ggsave(file.path(panels_dir, paste0("2016_reviewed_pool_icamp_matched_habitats_", process_design$code[i], ".pdf")),
    matched_plot, width = comparison_width, height = comparison_height, bg = "white")
}
ggsave(file.path(figures_dir, "2016_reviewed_pool_icamp_matched_habitats.pdf"),
  wrap_plots(matched_panels, ncol = 3, guides = "collect") & theme(legend.position = "bottom"),
  width = 15.6, height = 9, bg = "white")

# 8. Confirm input preservation and record the plotting environment ---------------
stopifnot(identical(input_checksums$md5, unname(tools::md5sum(input_paths))))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "figure_session_info.txt"))
message("Reviewed-pool iCAMP and within-habitat figures complete: ", output_dir)
