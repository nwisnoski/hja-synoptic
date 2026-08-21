# Explore 2016 microbial diversity and composition along catchment gradients.
# Run this file from the hja-synoptic repository root.

library(vegan)
library(ggplot2)

data_dir <- "data/derived/microbial_diversity_2016"
results_dir <- "results/diversity_2016/tables"
figures_dir <- "figures"

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
if (anyNA(metadata$sample_id) || !all(rowSums(counts) == 10000)) {
  stop("The 10K count table and sample metadata are not aligned.")
}

metadata$log_drainage_area <- log10(metadata$site_drainage_area_ha)
metadata$distance_to_outlet_km <- metadata$site_distance_to_outlet_m / 1000
hellinger <- decostand(counts, method = "hellinger")

habitat_order <- c("planktonic", "hyporheic", "sediment")
habitat_labels <- c(
  planktonic = "Planktonic streamwater",
  hyporheic = "Hyporheic porewater",
  sediment = "Stream sediment"
)
habitat_colors <- c(
  planktonic = "#2C7FB8",
  hyporheic = "#41B6C4",
  sediment = "#8C6D31"
)

# Join Hill diversity to the network variables without fitting models yet.
alpha <- read.csv(file.path(results_dir, "alpha_hill_numbers.csv"))
alpha <- cbind(
  alpha,
  metadata[match(alpha$sample_id, metadata$sample_id), c(
    "site_stream_order", "log_drainage_area", "distance_to_outlet_km"
  )]
)
alpha <- alpha[alpha$habitat %in% habitat_order, ]
alpha$habitat <- factor(alpha$habitat, levels = habitat_order)
write.csv(
  alpha,
  file.path(results_dir, "alpha_network_exploration.csv"),
  row.names = FALSE
)

alpha_plot_data <- rbind(
  data.frame(
    alpha,
    gradient = "Drainage area (log10 ha)",
    position = alpha$log_drainage_area
  ),
  data.frame(
    alpha,
    gradient = "Distance to outlet (km)",
    position = alpha$distance_to_outlet_km
  )
)
alpha_plot <- ggplot(alpha_plot_data, aes(position, q1, color = habitat)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_grid(habitat ~ gradient, scales = "free_x") +
  scale_color_manual(values = habitat_colors) +
  labs(x = NULL, y = "Hill diversity, q = 1") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(
  file.path(figures_dir, "2016_alpha_q1_network_gradients.pdf"),
  alpha_plot,
  width = 8,
  height = 7
)

# Make a separate unconstrained Hellinger PCA for each habitat. This lets us
# look within habitats before deciding which environmental models are useful.
ordination_scores <- list()
for (habitat in habitat_order) {
  rows <- metadata$habitat == habitat
  community <- hellinger[rows, colSums(counts[rows, , drop = FALSE]) > 0]
  pca <- rda(community)
  sites <- scores(pca, display = "sites", choices = 1:2, scaling = 1)
  variance <- eigenvals(pca) / sum(eigenvals(pca))

  ordination_scores[[habitat]] <- data.frame(
    sample_id = metadata$sample_id[rows],
    site_code = metadata$site_code[rows],
    habitat = habitat,
    pc1 = sites[, 1],
    pc2 = sites[, 2],
    pc1_percent = 100 * variance[1],
    pc2_percent = 100 * variance[2],
    utm_x_m = metadata$site_utm_x_m[rows],
    utm_y_m = metadata$site_utm_y_m[rows],
    stream_order = metadata$site_stream_order[rows],
    log_drainage_area = metadata$log_drainage_area[rows],
    distance_to_outlet_km = metadata$distance_to_outlet_km[rows]
  )
}
ordination_scores <- do.call(rbind, ordination_scores)
ordination_scores$habitat <- factor(
  ordination_scores$habitat,
  levels = habitat_order
)
write.csv(
  ordination_scores,
  file.path(results_dir, "community_pca_scores_by_habitat.csv"),
  row.names = FALSE
)

community_plot_data <- rbind(
  data.frame(
    ordination_scores,
    gradient = "Drainage area (log10 ha)",
    position = ordination_scores$log_drainage_area
  ),
  data.frame(
    ordination_scores,
    gradient = "Distance to outlet (km)",
    position = ordination_scores$distance_to_outlet_km
  )
)
community_plot <- ggplot(
  community_plot_data,
  aes(position, pc1, color = habitat)
) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_grid(habitat ~ gradient, scales = "free") +
  scale_color_manual(values = habitat_colors) +
  labs(x = NULL, y = "Within-habitat Hellinger PC1 score") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(
  file.path(figures_dir, "2016_community_pc1_network_gradients.pdf"),
  community_plot,
  width = 8,
  height = 7
)

spatial_plot <- ggplot(
  ordination_scores,
  aes(utm_x_m, utm_y_m, color = pc1)
) +
  geom_point(size = 3) +
  facet_wrap(~habitat, labeller = as_labeller(habitat_labels)) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  coord_equal() +
  labs(x = "UTM easting (m)", y = "UTM northing (m)", color = "PC1") +
  theme_bw()
ggsave(
  file.path(figures_dir, "2016_community_pc1_spatial_patterns.pdf"),
  spatial_plot,
  width = 10,
  height = 4
)

# Plot sediment community position against the four measured EEA rates. These
# are descriptive views only; formal models can follow after we inspect them.
sediment_scores <- ordination_scores[ordination_scores$habitat == "sediment", ]
sediment_metadata <- metadata[match(sediment_scores$sample_id, metadata$sample_id), ]
eea_names <- c("NAG", "LAP", "GLU", "AP")
eea_values <- sediment_metadata[c(
  "sediment_nag_umol_g_hr", "sediment_lap_umol_g_hr",
  "sediment_glu_umol_g_hr", "sediment_ap_umol_g_hr"
)]

eea_plot_data <- data.frame(
  sample_id = rep(sediment_scores$sample_id, times = 4),
  pc1 = rep(sediment_scores$pc1, times = 4),
  enzyme = factor(rep(eea_names, each = nrow(sediment_scores)), levels = eea_names),
  activity = unlist(eea_values, use.names = FALSE)
)
eea_plot <- ggplot(eea_plot_data, aes(pc1, activity)) +
  geom_point(color = habitat_colors["sediment"], size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, color = "grey30") +
  facet_wrap(~enzyme, scales = "free_y") +
  labs(
    x = "Sediment Hellinger PC1 score",
    y = "EEA (umol g-1 hr-1)"
  ) +
  theme_bw()
ggsave(
  file.path(figures_dir, "2016_sediment_community_eea_exploration.pdf"),
  eea_plot,
  width = 8,
  height = 6
)

message("Exploratory network and EEA figures complete.")
