# Basic 2016 microbial alpha, gamma, and beta diversity.
# Run 01_prepare_diversity.R first, then run this file from the repository root.

library(vegan)
library(iNEXT)
library(ggplot2)

data_dir <- "data/derived/microbial_diversity_2016"
results_dir <- "results/diversity_2016/tables"
figures_dir <- "figures"
rarefaction_depth <- 10000L

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  file.path(data_dir, "sample_metadata.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)

unrarefied_table <- read.csv(
  file.path(data_dir, "asv_counts_unrarefied.csv"),
  check.names = FALSE
)
rownames(unrarefied_table) <- unrarefied_table$asv_id
unrarefied_table$asv_id <- NULL
unrarefied_counts <- t(as.matrix(unrarefied_table))
storage.mode(unrarefied_counts) <- "integer"

counts_table <- read.csv(
  file.path(data_dir, "asv_counts_10k.csv"),
  check.names = FALSE
)
rownames(counts_table) <- counts_table$asv_id
counts_table$asv_id <- NULL
counts <- t(as.matrix(counts_table))
storage.mode(counts) <- "integer"

metadata <- metadata[match(rownames(counts), metadata$sample_id), ]
if (anyNA(metadata$sample_id) ||
    !identical(rownames(counts), rownames(unrarefied_counts)) ||
    !all(rowSums(counts) == rarefaction_depth)) {
  stop("The prepared count and metadata tables are not aligned.")
}

habitat_order <- c("planktonic", "hyporheic", "sediment", "soil")
habitat_labels <- c(
  planktonic = "Planktonic streamwater",
  hyporheic = "Hyporheic porewater",
  sediment = "Stream sediment",
  soil = "Terrestrial soil"
)
habitat_colors <- c(
  planktonic = "#2C7FB8",
  hyporheic = "#41B6C4",
  sediment = "#8C6D31",
  soil = "#4D9221"
)
metadata$habitat <- factor(metadata$habitat, levels = habitat_order)

# Alpha Hill numbers from the same rarefied table used by all later analyses.
alpha <- data.frame(
  sample_id = metadata$sample_id,
  site_code = metadata$site_code,
  habitat = as.character(metadata$habitat),
  original_reads = metadata$sequence_reads,
  q0 = specnumber(counts),
  q1 = exp(diversity(counts, index = "shannon")),
  q2 = diversity(counts, index = "invsimpson")
)
write.csv(
  alpha,
  file.path(results_dir, "alpha_hill_numbers.csv"),
  row.names = FALSE
)

# iNEXT gives the expected Hill numbers at exactly 10K reads without relying
# on the particular random rarefaction saved by file 1.
inext_input <- lapply(seq_len(nrow(unrarefied_counts)), function(i) {
  x <- unrarefied_counts[i, ]
  x[x > 0]
})
names(inext_input) <- rownames(unrarefied_counts)
inext_alpha <- estimateD(
  inext_input,
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "size",
  level = rarefaction_depth,
  nboot = 0
)
names(inext_alpha)[names(inext_alpha) == "Assemblage"] <- "sample_id"
inext_alpha$site_code <- metadata$site_code[
  match(inext_alpha$sample_id, metadata$sample_id)
]
inext_alpha$habitat <- as.character(metadata$habitat)[
  match(inext_alpha$sample_id, metadata$sample_id)
]
write.csv(
  inext_alpha,
  file.path(results_dir, "alpha_inext_hill_at_10000.csv"),
  row.names = FALSE
)

# Habitat-level gamma diversity from the 10K table.
habitat_rows <- split(seq_len(nrow(counts)), metadata$habitat)
gamma_input <- lapply(habitat_rows, function(i) {
  t(counts[i, , drop = FALSE] > 0)
})
common_coverage <- min(DataInfo(gamma_input, datatype = "incidence_raw")$SC)
gamma <- estimateD(
  gamma_input,
  q = c(0, 1, 2),
  datatype = "incidence_raw",
  base = "coverage",
  level = common_coverage,
  nboot = 0
)
names(gamma)[names(gamma) == "Assemblage"] <- "habitat"
write.csv(
  gamma,
  file.path(results_dir, "gamma_hill_numbers.csv"),
  row.names = FALSE
)

# Hellinger-transform the 10K table for PCA and habitat RDA.
hellinger <- decostand(counts, method = "hellinger")
pca <- rda(hellinger)
pca_sites <- scores(pca, display = "sites", choices = 1:2, scaling = 1)
pca_variance <- eigenvals(pca) / sum(eigenvals(pca))
pca_scores <- data.frame(
  sample_id = rownames(pca_sites),
  habitat = as.character(metadata$habitat),
  pc1 = pca_sites[, 1],
  pc2 = pca_sites[, 2]
)
write.csv(
  pca_scores,
  file.path(results_dir, "hellinger_pca_scores.csv"),
  row.names = FALSE
)

habitat_rda <- rda(hellinger ~ habitat, data = metadata)
rda_fit <- RsquareAdj(habitat_rda)
write.csv(
  data.frame(
    samples = nrow(metadata),
    asvs = ncol(counts),
    rarefaction_depth = rarefaction_depth,
    r_squared = unname(rda_fit$r.squared),
    adjusted_r_squared = unname(rda_fit$adj.r.squared)
  ),
  file.path(results_dir, "hellinger_rda_habitat_summary.csv"),
  row.names = FALSE
)

# Figures.
alpha_long <- data.frame(
  alpha[rep(seq_len(nrow(alpha)), 3), c("sample_id", "site_code", "habitat")],
  order = factor(rep(c("q0", "q1", "q2"), each = nrow(alpha))),
  hill_number = c(alpha$q0, alpha$q1, alpha$q2)
)
hill_labels <- c(q0 = "q = 0 (richness)", q1 = "q = 1", q2 = "q = 2")

alpha_plot <- ggplot(alpha_long, aes(habitat, hill_number, color = habitat)) +
  geom_boxplot(outlier.shape = NA, color = "grey40") +
  geom_jitter(width = 0.14, alpha = 0.7, size = 1.5) +
  facet_wrap(~order, scales = "free_y", labeller = as_labeller(hill_labels)) +
  scale_color_manual(values = habitat_colors) +
  scale_x_discrete(labels = habitat_labels) +
  labs(x = NULL, y = "Effective number of ASVs") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
ggsave(
  file.path(figures_dir, "2016_alpha_hill_by_habitat.pdf"),
  alpha_plot,
  width = 9,
  height = 4.5
)

inext_alpha$order <- factor(
  paste0("q", inext_alpha$Order.q),
  levels = c("q0", "q1", "q2")
)
inext_plot <- ggplot(inext_alpha, aes(habitat, qD, color = habitat)) +
  geom_boxplot(outlier.shape = NA, color = "grey40") +
  geom_jitter(width = 0.14, alpha = 0.7, size = 1.5) +
  facet_wrap(~order, scales = "free_y", labeller = as_labeller(hill_labels)) +
  scale_color_manual(values = habitat_colors) +
  scale_x_discrete(labels = habitat_labels) +
  labs(
    x = NULL,
    y = "Effective number of ASVs",
    caption = "iNEXT interpolation to exactly 10,000 reads per sample."
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
ggsave(
  file.path(figures_dir, "2016_inext_alpha_at_10000_by_habitat.pdf"),
  inext_plot,
  width = 9,
  height = 4.5
)

pca_plot <- ggplot(pca_scores, aes(pc1, pc2, color = habitat)) +
  geom_point(size = 2.2, alpha = 0.8) +
  scale_color_manual(values = habitat_colors, labels = habitat_labels) +
  labs(
    x = paste0("PC1 (", round(100 * pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_variance[2], 1), "%)"),
    color = "Habitat"
  ) +
  theme_bw()
ggsave(
  file.path(figures_dir, "2016_hellinger_pca_by_habitat.pdf"),
  pca_plot,
  width = 8,
  height = 6
)

message("Basic diversity analysis complete.")
