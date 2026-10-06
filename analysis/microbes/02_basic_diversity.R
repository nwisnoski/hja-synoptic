# Basic 2016 microbial alpha, gamma, and beta diversity.
# Run 01_prepare_diversity.R first; no command-line arguments are needed.

library(here)
library(vegan)
library(iNEXT)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme.
data_dir <- here("data", "derived", "microbial_diversity_2016")
results_dir <- here("results", "diversity_2016", "tables")
figures_dir <- here("figures")
rarefaction_depth <- 10000L

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

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
alpha_axis_labels <- c(
  planktonic = "Planktonic",
  hyporheic = "Hyporheic",
  sediment = "Sediment",
  soil = "Soil"
)
hill_orders <- c(0, 1, 2)
# Fix only figure jitter; diversity estimates use nboot = 0 and are deterministic.
jitter_seed <- 2016L
hill_labels <- c(q0 = "q = 0 (richness)", q1 = "q = 1", q2 = "q = 2")

# 2. Read and align the fixed-depth counts and sample metadata.
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

metadata$habitat <- factor(metadata$habitat, levels = habitat_order)

# 3. Calculate alpha Hill numbers from the same rarefied table used by all later analyses.
alpha <- data.frame(
  sample_id = metadata$sample_id,
  site_code = metadata$site_code,
  habitat = as.character(metadata$habitat),
  original_reads = metadata$target_reads,
  q0 = specnumber(counts),
  q1 = exp(diversity(counts, index = "shannon")),
  q2 = diversity(counts, index = "invsimpson")
)
write.csv(
  alpha,
  file.path(results_dir, "alpha_hill_numbers.csv"),
  row.names = FALSE
)

# 4. Calculate expected alpha diversity without random rarefaction.
# iNEXT gives the expected Hill numbers at exactly 10K reads without relying
# on the particular random rarefaction saved by file 1.
inext_input <- vector("list", nrow(unrarefied_counts))
names(inext_input) <- rownames(unrarefied_counts)
for (i in seq_len(nrow(unrarefied_counts))) {
  sample_counts <- unrarefied_counts[i, ]
  inext_input[[i]] <- sample_counts[sample_counts > 0]
}
inext_alpha <- estimateD(
  inext_input,
  q = hill_orders,
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

# 5. Calculate habitat-level gamma diversity from the 10K table.
habitat_rows <- split(seq_len(nrow(counts)), metadata$habitat)
gamma_input <- vector("list", length(habitat_rows))
names(gamma_input) <- names(habitat_rows)
for (habitat in names(habitat_rows)) {
  rows <- habitat_rows[[habitat]]
  gamma_input[[habitat]] <- t(counts[rows, , drop = FALSE] > 0)
}
common_coverage <- min(DataInfo(gamma_input, datatype = "incidence_raw")$SC)
gamma <- estimateD(
  gamma_input,
  q = hill_orders,
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

# 6. Hellinger-transform the 10K table for PCA and habitat RDA.
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

# 7. Plot alpha diversity and Hellinger PCA.
alpha_long <- data.frame(
  alpha[rep(seq_len(nrow(alpha)), 3), c("sample_id", "site_code", "habitat")],
  order = factor(rep(c("q0", "q1", "q2"), each = nrow(alpha))),
  hill_number = c(alpha$q0, alpha$q1, alpha$q2)
)

alpha_long$habitat <- factor(alpha_long$habitat, levels = habitat_order)

# Save each Hill order separately, then combine those same panels.
alpha_panels <- list()
for (hill_order in names(hill_labels)) {
  plot_data <- alpha_long[alpha_long$order == hill_order, ]
  alpha_panels[[hill_order]] <- ggplot(plot_data, aes(habitat, hill_number, color = habitat)) +
    geom_boxplot(outlier.shape = NA, color = "grey40") +
    geom_point(
      position = position_jitter(width = 0.14, seed = jitter_seed),
      alpha = 0.7,
      size = 1.5
    ) +
    scale_color_manual(values = habitat_colors) +
    scale_x_discrete(labels = alpha_axis_labels) +
    labs(x = NULL, y = paste("Effective number of ASVs", hill_labels[hill_order], sep = "\n")) +
    plot_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  ggsave(
    file.path(figures_dir, paste0("2016_alpha_hill_", hill_order, "_by_habitat.pdf")),
    alpha_panels[[hill_order]],
    width = 6,
    height = 4.8,
    bg = "white"
  )
}
alpha_plot <- wrap_plots(alpha_panels, nrow = 1)
ggsave(
  file.path(figures_dir, "2016_alpha_hill_by_habitat.pdf"),
  alpha_plot,
  width = 11,
  height = 4.8,
  bg = "white"
)

inext_alpha$order <- factor(
  paste0("q", inext_alpha$Order.q),
  levels = c("q0", "q1", "q2")
)
inext_alpha$habitat <- factor(inext_alpha$habitat, levels = habitat_order)
inext_panels <- list()
for (hill_order in names(hill_labels)) {
  plot_data <- inext_alpha[inext_alpha$order == hill_order, ]
  inext_panels[[hill_order]] <- ggplot(plot_data, aes(habitat, qD, color = habitat)) +
    geom_boxplot(outlier.shape = NA, color = "grey40") +
    geom_point(
      position = position_jitter(width = 0.14, seed = jitter_seed),
      alpha = 0.7,
      size = 1.5
    ) +
    scale_color_manual(values = habitat_colors) +
    scale_x_discrete(labels = alpha_axis_labels) +
    labs(x = NULL, y = paste("Effective number of ASVs", hill_labels[hill_order], sep = "\n")) +
    plot_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  ggsave(
    file.path(figures_dir, paste0("2016_inext_alpha_", hill_order, "_at_10000_by_habitat.pdf")),
    inext_panels[[hill_order]],
    width = 6,
    height = 4.8,
    bg = "white"
  )
}
inext_plot <- wrap_plots(inext_panels, nrow = 1)
ggsave(
  file.path(figures_dir, "2016_inext_alpha_at_10000_by_habitat.pdf"),
  inext_plot,
  width = 11,
  height = 4.8,
  bg = "white"
)

pca_scores$habitat <- factor(pca_scores$habitat, levels = habitat_order)
pca_plot <- ggplot(pca_scores, aes(pc1, pc2, color = habitat)) +
  geom_point(size = 2.2, alpha = 0.8) +
  coord_equal() +
  scale_color_manual(values = habitat_colors, labels = habitat_labels) +
  labs(
    x = paste0("PC1 (", round(100 * pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_variance[2], 1), "%)"),
    color = "Habitat"
  ) +
  plot_theme
ggsave(
  file.path(figures_dir, "2016_hellinger_pca_by_habitat.pdf"),
  pca_plot,
  width = 8,
  height = 6,
  bg = "white"
)

message("Basic diversity analysis complete.")
