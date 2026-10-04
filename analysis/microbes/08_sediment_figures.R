# Sediment-only descriptive figures from the fixed 10K microbial table.

library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(vegan)

# 1. Settings and plotting theme ------------------------------------------------
data_dir <- here("data", "derived", "microbial_diversity_2016")
results_dir <- here("results", "diversity_2016", "tables")
figures_dir <- here("figures")
target_depth <- 10000
figure_width <- 6
figure_height <- 4.8
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
    strip.background = element_blank())

# 2. Align the sediment subset --------------------------------------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
metadata <- read_csv(file.path(data_dir, "sample_metadata.csv"),
  col_types = cols(site_code = col_character()))
count_table <- read_csv(file.path(data_dir, "asv_counts_10k.csv"), show_col_types = FALSE)
stopifnot(!anyDuplicated(count_table$asv_id), !anyDuplicated(metadata$sample_id))
counts <- t(as.matrix(count_table[, -1]))
colnames(counts) <- count_table$asv_id
storage.mode(counts) <- "numeric"
stopifnot(all(is.finite(counts)), all(counts >= 0), all(rowSums(counts) == target_depth))
metadata <- metadata[match(rownames(counts), metadata$sample_id), ]
stopifnot(!anyNA(metadata$sample_id), identical(rownames(counts), metadata$sample_id))
sediment_rows <- metadata$habitat == "sediment"
sediment_counts <- counts[sediment_rows, , drop = FALSE]
sediment_counts <- sediment_counts[, colSums(sediment_counts) > 0, drop = FALSE]
sediment_metadata <- metadata[sediment_rows, ]
stopifnot(all(is.finite(sediment_metadata$site_drainage_area_ha)),
  all(sediment_metadata$site_drainage_area_ha > 0))

# 3. Ordinate composition and calculate Hill diversity --------------------------
hellinger <- decostand(sediment_counts, method = "hellinger")
pca <- prcomp(hellinger, center = TRUE, scale. = FALSE)
percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
sediment_summary <- sediment_metadata %>%
  transmute(sample_id, site_code, site_drainage_area_ha, site_stream_order,
    has_fticr_profile,
    q0 = rowSums(sediment_counts > 0),
    q1 = exp(diversity(sediment_counts, index = "shannon")),
    q2 = diversity(sediment_counts, index = "invsimpson"),
    pc1 = pca$x[, 1], pc2 = pca$x[, 2])

# Verify these summaries against the existing baseline rather than replace it.
baseline <- read_csv(file.path(results_dir, "alpha_hill_numbers.csv"),
  col_types = cols(site_code = col_character()))
baseline <- baseline[match(sediment_summary$sample_id, baseline$sample_id), ]
stopifnot(!anyNA(baseline$sample_id))
stopifnot(isTRUE(all.equal(as.matrix(sediment_summary[c("q0", "q1", "q2")]),
  as.matrix(baseline[c("q0", "q1", "q2")]), check.attributes = FALSE, tolerance = 1e-8)))
write_csv(sediment_summary, file.path(results_dir, "sediment_10k_figure_summary.csv"))
write_csv(tibble(sites = nrow(sediment_counts), asvs = ncol(sediment_counts),
  pc1_percent = percent[1], pc2_percent = percent[2]),
  file.path(results_dir, "sediment_10k_pca_summary.csv"))

# 4. Save individual panels -----------------------------------------------------
pca_plot <- ggplot(sediment_summary, aes(pc1, pc2, color = log10(site_drainage_area_ha))) +
  geom_point(size = 2.8) +
  scale_color_viridis_c(name = "Drainage area\n(log10 ha)") +
  coord_equal() +
  labs(x = sprintf("Hellinger PC1 (%.1f%%)", percent[1]),
    y = sprintf("Hellinger PC2 (%.1f%%)", percent[2])) +
  plot_theme
ggsave(file.path(figures_dir, "2016_sediment_microbes_hellinger_pca.pdf"),
  pca_plot, width = figure_width, height = figure_height, bg = "white")

alpha_plot <- ggplot(sediment_summary, aes(site_drainage_area_ha, q1)) +
  geom_point(size = 2.8, color = "#0072B2") +
  scale_x_log10() +
  labs(x = "Drainage area (ha)", y = "Hill diversity (q = 1), 10K reads") +
  plot_theme
ggsave(file.path(figures_dir, "2016_sediment_microbes_q1_drainage.pdf"),
  alpha_plot, width = figure_width, height = figure_height, bg = "white")
session_lines <- capture.output(sessionInfo())
writeLines(trimws(session_lines, which = "right"),
  file.path(results_dir, "session_info_sediment_figures.txt"))
message("Sediment figures: ", nrow(sediment_counts), " samples; ", ncol(sediment_counts), " ASVs.")
