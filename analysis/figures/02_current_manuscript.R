# Plot habitat diversity and weighted UniFrac for the reviewed manuscript cohort.
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(vegan)

# 1. Settings and plotting theme ------------------------------------------------
habitats <- c("soil", "hyporheic", "sediment")
habitat_labels <- c(soil = "Soil", hyporheic = "Hyporheic", sediment = "Sediment")
habitat_colors <- c(soil = "#D55E00", hyporheic = "#0072B2", sediment = "#009E73")
habitat_shapes <- c(soil = 17, hyporheic = 15, sediment = 16)
jitter_seed <- 20261010
results_dir <- here("results", "manuscript_figures_2016")
summary_file <- here("results", "diversity_2016", "soil_porewater_sediment", "sample_summary.csv")
distance_file <- here("results", "diversity_2016", "tables", "unifrac_weighted_10k.csv")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    strip.background = element_blank(), legend.position = "bottom")
theme_set(plot_theme)

# 2. Read the reviewed communities and saved phylogenetic distances -------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input_files <- c(summary_file, distance_file)
input_hashes <- unname(tools::md5sum(input_files))
samples <- read_csv(summary_file, show_col_types = FALSE)
stopifnot(nrow(samples) == 76, !anyDuplicated(samples$sample_id),
  all(samples$included_manuscript), !any(samples$excluded_suspect_library),
  sum(samples$habitat == "soil") == 15, sum(samples$habitat == "sediment") == 34,
  sum(samples$habitat == "hyporheic") == 27)
samples$habitat <- factor(samples$habitat, levels = habitats)
distance_table <- read_csv(distance_file, show_col_types = FALSE, name_repair = "minimal")
distance_matrix <- as.matrix(distance_table[, -1])
rownames(distance_matrix) <- distance_table$sample_id
stopifnot(all(samples$sample_id %in% rownames(distance_matrix)),
  all(samples$sample_id %in% colnames(distance_matrix)))
distance_matrix <- distance_matrix[samples$sample_id, samples$sample_id]
stopifnot(all(is.finite(distance_matrix)), all(distance_matrix >= 0),
  max(abs(distance_matrix - t(distance_matrix))) < 1e-10,
  max(abs(diag(distance_matrix))) < 1e-10)

# 3. Combine the three fixed-depth Hill diversity panels -------------------------
# Points are libraries, including regional soils rather than paired soil sites.
alpha_panels <- list()
alpha_labels <- c(q0 = "ASV richness (q = 0)", q1 = "Hill diversity (q = 1)",
  q2 = "Hill diversity (q = 2)")
for (metric in names(alpha_labels)) {
  alpha_panels[[metric]] <- ggplot(samples, aes(habitat, .data[[metric]], color = habitat)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.4) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = jitter_seed),
      size = 1.6, alpha = 0.7) +
    scale_color_manual(values = habitat_colors, labels = habitat_labels) +
    scale_x_discrete(labels = c(soil = "Soil\n(n = 15)", hyporheic = "Hyporheic\n(n = 27)",
      sediment = "Sediment\n(n = 34)")) +
    labs(x = NULL, y = alpha_labels[[metric]]) + theme(legend.position = "none")
}
alpha_plot <- wrap_plots(alpha_panels, nrow = 1)
ggsave(here("figures", "2016_manuscript_alpha_diversity.pdf"), alpha_plot,
  width = 10, height = 3.6, bg = "white")

# 4. Refit PCoA on exactly the reviewed 76-library distance submatrix ------------
# Reuse the original midpoint-rooted tree distances; do not rebuild or reroot
# the phylogeny. The original 100-library PCoA scores are not reused.
ordination <- wcmdscale(as.dist(distance_matrix), k = 2, eig = TRUE, add = "lingoes")
variance <- ordination$eig[1:2] / sum(ordination$eig[ordination$eig > 0])
scores <- samples %>% select(sample_id, site_code, habitat, network_group) %>%
  mutate(pcoa1 = ordination$points[sample_id, 1], pcoa2 = ordination$points[sample_id, 2])
write_csv(scores, file.path(results_dir, "weighted_unifrac_scores.csv"))
write_csv(tibble(samples = nrow(samples), distance = "Normalized weighted UniFrac",
  correction = "Lingoes", correction_constant = ordination$ac,
  axis1_percent = 100 * variance[1], axis2_percent = 100 * variance[2],
  distance_reference = "Original fixed-10K table; tree midpoint-rooted after pruning"),
  file.path(results_dir, "weighted_unifrac_summary.csv"))
ordination_plot <- ggplot(scores, aes(pcoa1, pcoa2, color = habitat, shape = habitat)) +
  geom_point(size = 2.4, alpha = 0.8) + coord_equal() +
  scale_color_manual(values = habitat_colors, labels = habitat_labels) +
  scale_shape_manual(values = habitat_shapes, labels = habitat_labels) +
  labs(x = sprintf("PCoA1 (%.1f%%)", 100 * variance[1]),
    y = sprintf("PCoA2 (%.1f%%)", 100 * variance[2]), color = NULL, shape = NULL)
ggsave(here("figures", "2016_manuscript_weighted_unifrac.pdf"), ordination_plot,
  width = 6.4, height = 5.2, bg = "white")

# 5. Record unchanged inputs and the plotting session ----------------------------
stopifnot(identical(input_hashes, unname(tools::md5sum(input_files))))
write_csv(tibble(path = input_files, md5 = input_hashes, unchanged = TRUE),
  file.path(results_dir, "input_checksums.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(results_dir, "session_info.txt"))
message("Current manuscript alpha diversity and weighted UniFrac figures complete.")
