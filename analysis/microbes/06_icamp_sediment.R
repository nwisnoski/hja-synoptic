# iCAMP analysis of the 2016 stream-sediment communities.
# Run from the repository root. This is intended for the computing cluster.

library(ape)
library(phangorn)
library(iCAMP)
library(here)
library(castor)
library(bigmemory)

# 1. Settings. A pilot keeps the full community but uses only 100 null draws.
test_run <- FALSE
data_dir <- here("data", "derived", "microbial_diversity_2016")
tree_file <- here("results", "phylogeny_2016", "asv_tree_screened_fasttree.nwk")
output_dir <- here("results", "icamp_2016", "sediment_castor")
pd_dir <- file.path(output_dir, "phylogenetic_distance")
if (test_run) output_dir <- file.path(output_dir, "pilot")

rand_time <- 1000
if (test_run) rand_time <- 100
random_seed <- 2016
nworker <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
memory_gb <- 150
distance_block_size <- 256
distance_tolerance <- 1e-8

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pd_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
pd_dir <- normalizePath(pd_dir)

if (file.exists(file.path(output_dir, "icamp_result.rds"))) {
  stop("The completed iCAMP result already exists: ", output_dir)
}

# 2. Load the fixed 10K community table and prune/root the sediment tree.
metadata <- read.csv(
  file.path(data_dir, "sample_metadata.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
counts_table <- read.csv(
  file.path(data_dir, "asv_counts_10k.csv"),
  check.names = FALSE
)
taxonomy <- read.csv(
  file.path(data_dir, "asv_metadata_screened.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(counts_table) <- counts_table$asv_id
counts_table$asv_id <- NULL
community <- t(as.matrix(counts_table))
storage.mode(community) <- "integer"

metadata <- metadata[match(rownames(community), metadata$sample_id), ]
if (anyNA(metadata$sample_id) || anyNA(metadata$habitat)) {
  stop("Some 10K samples lack sample or habitat metadata.")
}
sediment_rows <- metadata$habitat == "sediment"
community <- community[sediment_rows, , drop = FALSE]
metadata <- metadata[sediment_rows, , drop = FALSE]
community <- community[, colSums(community) > 0, drop = FALSE]

if (anyNA(metadata$sample_id) ||
    !identical(rownames(community), metadata$sample_id) ||
    !all(rowSums(community) == 10000)) {
  stop("The sediment 10K community table and metadata are not aligned.")
}

tree <- read.tree(tree_file)
if (!all(colnames(community) %in% tree$tip.label)) {
  stop("Some sediment ASVs are absent from the phylogeny.")
}
tree <- keep.tip(tree, colnames(community))
tree <- midpoint(tree)
tree$node.label <- NULL

# 3. Calculate tip distances once and save a validated disk-backed matrix.
distance_inputs <- data.frame(
  input = c("asv_counts_10k.csv", basename(tree_file)),
  md5 = unname(tools::md5sum(c(file.path(data_dir, "asv_counts_10k.csv"), tree_file))),
  method = "castor branch-length tip distances v1"
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

# 4. Record samples/settings and run the unchanged iCAMP null model.
write.csv(
  metadata[, c(
    "sample_id", "site_code", "site_stream_order", "site_drainage_area_ha",
    "site_distance_to_outlet_m", "site_utm_x_m", "site_utm_y_m"
  )],
  file.path(output_dir, "sediment_samples.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    setting = c(
      "samples", "asvs", "reads_per_sample", "randomizations", "random_seed",
      "workers", "phylogenetic_metric", "significance_index",
      "phylogenetic_randomization", "taxonomic_randomization",
      "minimum_bin_size", "phylogenetic_signal_distance", "rooting"
    ),
    value = c(
      nrow(community), ncol(community), 10000, rand_time, random_seed,
      nworker, "bMPD", "Confidence", "within.bin", "across.all",
      24, 0.2, "midpoint root after pruning to sediment ASVs"
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
  prefix = "hja2016_sediment",
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
  qp.save = TRUE,
  detail.null = FALSE,
  ignore.zero = TRUE,
  output.wd = output_dir,
  unit.sum = rowSums(community),
  temp.save = TRUE
)
saveRDS(icamp_result, file.path(output_dir, "icamp_result.rds"))

# 5. Export process fractions and bin taxonomy.
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
taxonomy <- taxonomy[match(taxa_bins$asv_id, taxonomy$asv_id), ]
taxa_bins <- cbind(taxa_bins, taxonomy[setdiff(names(taxonomy), "asv_id")])
write.csv(
  taxa_bins,
  file.path(output_dir, "asv_bin_taxonomy.csv"),
  row.names = FALSE,
  na = ""
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)
message("Sediment iCAMP analysis complete: ", output_dir)
