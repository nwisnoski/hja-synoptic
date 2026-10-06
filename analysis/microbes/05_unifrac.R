# UniFrac analysis of the 2016 microbial communities.
# Use the completed master tree; this script does not rebuild the phylogeny.

library(ape)
library(phangorn)
library(phyloseq)
library(here)
library(vegan)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme.
data_dir <- here("data", "derived", "microbial_diversity_2016")
tree_file <- here("results", "phylogeny_2016", "asv_tree_screened_fasttree.nwk")
results_dir <- here("results", "diversity_2016", "tables")
figures_dir <- here("figures")

habitat_order <- c("planktonic", "hyporheic", "sediment", "soil")
habitat_labels <- c(
  planktonic = "Planktonic streamwater",
  hyporheic = "Hyporheic porewater",
  sediment = "Stream sediment",
  soil = "Terrestrial soil"
)
habitat_colors <- c(
  planktonic = "#0072B2",
  hyporheic = "#E69F00",
  sediment = "#009E73",
  soil = "#D55E00"
)
rarefaction_depth <- 10000L

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# 2. Read and align the fixed-depth counts and sample metadata.
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  file.path(data_dir, "sample_metadata.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
counts_table <- read.csv(
  file.path(data_dir, "asv_counts_10k.csv"),
  check.names = FALSE
)
rownames(counts_table) <- counts_table$asv_id
counts_table$asv_id <- NULL
counts <- t(as.matrix(counts_table))
storage.mode(counts) <- "integer"

metadata <- metadata[match(rownames(counts), metadata$sample_id), ]
if (anyNA(metadata$sample_id) || !all(rowSums(counts) == rarefaction_depth)) {
  stop("The 10K count table and sample metadata are not aligned.")
}

# 3. Prune and midpoint-root the completed tree in memory.
# FastTree estimates an unrooted tree, whereas UniFrac requires a root. Prune
# to the ASVs in the 10K table and apply a midpoint root in memory. The full
# tree on disk remains unchanged and is retained for later iCAMP analyses.
tree <- read.tree(tree_file)
if (!all(colnames(counts) %in% tree$tip.label)) {
  stop("Some ASVs in the 10K table are absent from the phylogeny.")
}
tree <- keep.tip(tree, colnames(counts))
tree <- midpoint(tree)

# 4. Calculate weighted and unweighted UniFrac.
physeq <- phyloseq(
  otu_table(counts, taxa_are_rows = FALSE),
  phy_tree(tree)
)
unweighted_unifrac <- UniFrac(physeq, weighted = FALSE)
weighted_unifrac <- UniFrac(physeq, weighted = TRUE, normalized = TRUE)

unweighted_matrix <- as.matrix(unweighted_unifrac)
weighted_matrix <- as.matrix(weighted_unifrac)
write.csv(
  data.frame(sample_id = rownames(unweighted_matrix), unweighted_matrix,
             check.names = FALSE),
  file.path(results_dir, "unifrac_unweighted_10k.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(sample_id = rownames(weighted_matrix), weighted_matrix,
             check.names = FALSE),
  file.path(results_dir, "unifrac_weighted_10k.csv"),
  row.names = FALSE
)

# 5. Ordinate distances with a Lingoes correction.
# UniFrac is not necessarily Euclidean, so use a Lingoes correction for PCoA.
unweighted_pcoa <- wcmdscale(
  unweighted_unifrac, k = 2, eig = TRUE, add = "lingoes"
)
weighted_pcoa <- wcmdscale(
  weighted_unifrac, k = 2, eig = TRUE, add = "lingoes"
)
unweighted_variance <- unweighted_pcoa$eig[1:2] /
  sum(unweighted_pcoa$eig[unweighted_pcoa$eig > 0])
weighted_variance <- weighted_pcoa$eig[1:2] /
  sum(weighted_pcoa$eig[weighted_pcoa$eig > 0])

unweighted_scores <- data.frame(
  sample_id = rownames(unweighted_pcoa$points),
  distance = "Unweighted UniFrac",
  pcoa1 = unweighted_pcoa$points[, 1],
  pcoa2 = unweighted_pcoa$points[, 2]
)
weighted_scores <- data.frame(
  sample_id = rownames(weighted_pcoa$points),
  distance = "Weighted UniFrac",
  pcoa1 = weighted_pcoa$points[, 1],
  pcoa2 = weighted_pcoa$points[, 2]
)
unifrac_scores <- rbind(unweighted_scores, weighted_scores)
unifrac_scores$site_code <- metadata$site_code[
  match(unifrac_scores$sample_id, metadata$sample_id)
]
unifrac_scores$habitat <- metadata$habitat[
  match(unifrac_scores$sample_id, metadata$sample_id)
]
write.csv(
  unifrac_scores,
  file.path(results_dir, "unifrac_pcoa_scores_10k.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    samples = nrow(counts),
    asvs = ncol(counts),
    reads_per_sample = rarefaction_depth,
    rooting = "midpoint root applied after pruning to the 10K table",
    unweighted_pcoa1_percent = 100 * unweighted_variance[1],
    unweighted_pcoa2_percent = 100 * unweighted_variance[2],
    weighted_pcoa1_percent = 100 * weighted_variance[1],
    weighted_pcoa2_percent = 100 * weighted_variance[2]
  ),
  file.path(results_dir, "unifrac_summary.csv"),
  row.names = FALSE
)

# 6. Plot the two ordinations.
unifrac_scores$habitat <- factor(
  unifrac_scores$habitat,
  levels = habitat_order
)
unifrac_panels <- list()
for (distance_name in c("Unweighted UniFrac", "Weighted UniFrac")) {
  plot_data <- unifrac_scores[unifrac_scores$distance == distance_name, ]
  axis_variance <- unweighted_variance
  file_label <- "unweighted"
  if (distance_name == "Weighted UniFrac") {
    axis_variance <- weighted_variance
    file_label <- "weighted"
  }
  unifrac_panels[[distance_name]] <- ggplot(plot_data, aes(pcoa1, pcoa2, color = habitat)) +
    geom_point(size = 2.2, alpha = 0.8) +
    coord_equal() +
    scale_color_manual(values = habitat_colors, labels = habitat_labels) +
    labs(
      x = paste0("PCoA1 (", round(100 * axis_variance[1], 1), "%)"),
      y = paste0("PCoA2 (", round(100 * axis_variance[2], 1), "%)"),
      color = "Habitat"
    ) +
    plot_theme
  ggsave(
    file.path(figures_dir, paste0("2016_unifrac_", file_label, "_pcoa_by_habitat.pdf")),
    unifrac_panels[[distance_name]],
    width = 6,
    height = 4.8,
    bg = "white"
  )
}
unifrac_plot <- wrap_plots(unifrac_panels, nrow = 1, guides = "collect")
ggsave(
  file.path(figures_dir, "2016_unifrac_pcoa_by_habitat.pdf"),
  unifrac_plot,
  width = 11,
  height = 4.8,
  bg = "white"
)

message("UniFrac analysis complete.")
