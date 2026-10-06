# Explore 2016 microbial diversity and composition along catchment gradients.
# Run directly in Positron from any working directory inside this project.

library(here)
library(vegan)
library(ggplot2)
library(sandwich)

# 1. Settings and plotting theme.
data_dir <- here("data", "derived", "microbial_diversity_2016")
results_dir <- here("results", "diversity_2016", "tables")
figures_dir <- here("figures")

habitat_order <- c("planktonic", "hyporheic", "sediment")
habitat_labels <- c(
  planktonic = "Planktonic streamwater",
  hyporheic = "Hyporheic porewater",
  sediment = "Stream sediment"
)
habitat_colors <- c(
  planktonic = "#0072B2",
  hyporheic = "#E69F00",
  sediment = "#009E73"
)

rarefaction_depth <- 10000L
eea_names <- c("NAG", "LAP", "GLU", "AP")

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

gradient_theme <- plot_theme +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    strip.placement = "outside",
    strip.text.x = element_text(size = 12)
  )

# This layout is reused for the three overviews and nine individual panels.
plot_drainage <- function(plot_data, line_data, response_label) {
  ggplot(plot_data, aes(log_drainage_area, value, color = habitat)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_segment(
      data = line_data,
      aes(
        x = log_drainage_min,
        xend = log_drainage_max,
        y = intercept + slope_per_log10_ha * log_drainage_min,
        yend = intercept + slope_per_log10_ha * log_drainage_max,
        color = habitat
      ),
      inherit.aes = FALSE,
      linewidth = 0.7
    ) +
    scale_color_manual(values = habitat_colors) +
    labs(x = "Drainage area (log10 ha)", y = response_label) +
    gradient_theme +
    theme(legend.position = "none")
}

response_order <- c("q1", "q2", "pc1")
response_labels <- c(
  q1 = "Hill diversity, q = 1",
  q2 = "Hill diversity, q = 2",
  pc1 = "Within-habitat Hellinger PC1 score"
)
response_files <- c(
  q1 = "2016_alpha_q1_network_gradients.pdf",
  q2 = "2016_alpha_q2_network_gradients.pdf",
  pc1 = "2016_community_pc1_network_gradients.pdf"
)
significance_level <- 0.05
p_adjustment <- "BH"
cluster_covariance_type <- "HC1"
confidence_level <- 0.95
regression_design <- expand.grid(
  response = response_order,
  habitat = habitat_order,
  stringsAsFactors = FALSE
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

metadata$log_drainage_area <- log10(metadata$site_drainage_area_ha)
metadata$distance_to_outlet_km <- metadata$site_distance_to_outlet_m / 1000
hellinger <- decostand(counts, method = "hellinger")

# 3. Join Hill diversity to the recorded network variables.
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

# 4. Make a separate unconstrained Hellinger PCA for each habitat. This lets us
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

# 5. Test drainage-area slopes, allowing residual correlation within segments.
# The response and predictor stay on the scales used in the original plots.
# Separate three-habitat test families are declared for q1, q2, and PC1.
model_columns <- c("sample_id", "habitat", "log_drainage_area")
drainage_data <- rbind(
  data.frame(alpha[model_columns], response = "q1", value = alpha$q1),
  data.frame(alpha[model_columns], response = "q2", value = alpha$q2),
  data.frame(ordination_scores[model_columns], response = "pc1", value = ordination_scores$pc1)
)
drainage_data$segment <- metadata$site_segment_number[
  match(drainage_data$sample_id, metadata$sample_id)
]
drainage_data$valid_response <- is.finite(drainage_data$value)
drainage_data$valid_drainage <- is.finite(drainage_data$log_drainage_area)
drainage_data$segment_resolved <- !is.na(drainage_data$segment)
drainage_data$included <- drainage_data$valid_response &
  drainage_data$valid_drainage & drainage_data$segment_resolved
write.csv(
  drainage_data,
  file.path(results_dir, "network_drainage_model_inclusion.csv"),
  row.names = FALSE
)

regression_summary <- data.frame(
  regression_design,
  n_samples = NA_integer_,
  n_excluded = NA_integer_,
  n_segments = NA_integer_,
  log_drainage_min = NA_real_,
  log_drainage_max = NA_real_,
  intercept = NA_real_,
  slope_per_log10_ha = NA_real_,
  slope_cluster_se = NA_real_,
  slope_ci_lower = NA_real_,
  slope_ci_upper = NA_real_,
  r_squared = NA_real_,
  t_df = NA_integer_,
  p_ols = NA_real_,
  p_raw = NA_real_,
  p_adjusted = NA_real_,
  status = "pending",
  line_shown = FALSE
)
regression_summary$habitat <- factor(regression_summary$habitat, levels = habitat_order)
for (i in seq_len(nrow(regression_design))) {
  response <- regression_design$response[i]
  habitat <- regression_design$habitat[i]
  rows <- drainage_data$response == response & drainage_data$habitat == habitat
  current_data <- drainage_data[rows & drainage_data$included, ]
  n_samples <- nrow(current_data)
  n_segments <- length(unique(current_data$segment))
  regression_summary$n_samples[i] <- n_samples
  regression_summary$n_excluded[i] <- sum(rows & !drainage_data$included)
  regression_summary$n_segments[i] <- n_segments
  if (n_samples < 3L || n_segments < 2L ||
      length(unique(current_data$log_drainage_area)) < 2L) {
    regression_summary$status[i] <- "insufficient_data"
    next
  }
  model <- lm(value ~ log_drainage_area, data = current_data)
  model_summary <- summary(model)
  clustered_covariance <- vcovCL(
    model,
    cluster = current_data$segment,
    type = cluster_covariance_type,
    cadjust = TRUE
  )
  slope <- unname(coef(model)["log_drainage_area"])
  slope_se <- sqrt(clustered_covariance["log_drainage_area", "log_drainage_area"])
  if (!is.finite(slope_se) || slope_se <= 0) {
    regression_summary$status[i] <- "invalid_slope_standard_error"
    next
  }
  t_df <- n_segments - 1L
  critical_t <- qt(1 - (1 - confidence_level) / 2, df = t_df)
  regression_summary$log_drainage_min[i] <- min(current_data$log_drainage_area)
  regression_summary$log_drainage_max[i] <- max(current_data$log_drainage_area)
  regression_summary$intercept[i] <- unname(coef(model)["(Intercept)"])
  regression_summary$slope_per_log10_ha[i] <- slope
  regression_summary$slope_cluster_se[i] <- slope_se
  regression_summary$slope_ci_lower[i] <- slope - critical_t * slope_se
  regression_summary$slope_ci_upper[i] <- slope + critical_t * slope_se
  regression_summary$r_squared[i] <- model_summary$r.squared
  regression_summary$t_df[i] <- t_df
  regression_summary$p_ols[i] <- model_summary$coefficients["log_drainage_area", "Pr(>|t|)"]
  regression_summary$p_raw[i] <- 2 * pt(-abs(slope / slope_se), df = t_df)
  regression_summary$status[i] <- "ok"
}
for (response in response_order) {
  rows <- regression_summary$response == response
  regression_summary$p_adjusted[rows] <- p.adjust(
    regression_summary$p_raw[rows],
    method = p_adjustment,
    n = length(habitat_order)
  )
}
regression_summary$line_shown <- regression_summary$status == "ok" &
  !is.na(regression_summary$p_adjusted) &
  regression_summary$p_adjusted < significance_level
regression_summary$significance_level <- significance_level
regression_summary$p_adjustment <- p_adjustment
regression_summary$covariance_type <- cluster_covariance_type
regression_summary$confidence_level <- confidence_level
write.csv(
  regression_summary,
  file.path(results_dir, "network_drainage_regressions.csv"),
  row.names = FALSE
)

# 6. Plot drainage area only, retaining all points and only supported lines.
# Each line is limited to that habitat's observed drainage-area range.
# All fitted slopes and intervals remain in the table, including hidden lines.
drainage_plots <- list()
for (response in response_order) {
  plot_data <- drainage_data[
    drainage_data$response == response & drainage_data$included,
  ]
  line_data <- regression_summary[
    regression_summary$response == response & regression_summary$line_shown,
  ]
  facet_scales <- "fixed"
  if (response == "pc1") {
    facet_scales <- "free_y"
  }
  drainage_plots[[response]] <- plot_drainage(
    plot_data,
    line_data,
    response_labels[response]
  ) +
    facet_grid(habitat ~ ., scales = facet_scales)
  ggsave(
    file.path(figures_dir, response_files[response]),
    drainage_plots[[response]],
    width = 6,
    height = 7,
    bg = "white"
  )
}

# 7. Map community position with the established shared PC1 color scale.
spatial_plot <- ggplot(
  ordination_scores,
  aes(utm_x_m, utm_y_m, color = pc1)
) +
  geom_point(size = 3) +
  facet_wrap(~habitat, labeller = as_labeller(habitat_labels)) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  coord_equal() +
  labs(x = "UTM easting (m)", y = "UTM northing (m)", color = "PC1") +
  plot_theme
ggsave(
  file.path(figures_dir, "2016_community_pc1_spatial_patterns.pdf"),
  spatial_plot,
  width = 10,
  height = 4,
  bg = "white"
)

# 8. Plot sediment community position against the four measured EEA rates.
# Show observations only until corresponding slope tests are specified.
sediment_scores <- ordination_scores[ordination_scores$habitat == "sediment", ]
sediment_metadata <- metadata[match(sediment_scores$sample_id, metadata$sample_id), ]
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
  facet_wrap(~enzyme, scales = "free_y") +
  labs(
    x = "Sediment Hellinger PC1 score",
    y = "EEA (umol g-1 hr-1)"
  ) +
  plot_theme
ggsave(
  file.path(figures_dir, "2016_sediment_community_eea_exploration.pdf"),
  eea_plot,
  width = 8,
  height = 6,
  bg = "white"
)

# 9. Save individual drainage, spatial, and enzyme panels.
# Each panel receives only its own points and significance-filtered line.
for (response in response_order) {
  for (habitat in habitat_order) {
    plot_rows <- drainage_data$response == response &
      drainage_data$habitat == habitat & drainage_data$included
    line_rows <- regression_summary$response == response &
      regression_summary$habitat == habitat & regression_summary$line_shown
    panel <- plot_drainage(
      drainage_data[plot_rows, ],
      regression_summary[line_rows, ],
      response_labels[response]
    )
    file_prefix <- paste0("2016_alpha_", response)
    if (response == "pc1") {
      file_prefix <- "2016_community_pc1"
    }
    ggsave(
      file.path(figures_dir, paste0(file_prefix, "_", habitat, "_drainage.pdf")),
      panel,
      width = 6,
      height = 4.8,
      bg = "white"
    )
  }
}
# A shared color range keeps PC1 colors comparable across the spatial panels.
pc1_range <- range(ordination_scores$pc1)
for (habitat in habitat_order) {
  habitat_rows <- ordination_scores$habitat == habitat
  spatial_panel <- (spatial_plot + ordination_scores[habitat_rows, ]) +
    facet_null() +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = pc1_range)
  ggsave(
    file.path(figures_dir, paste0("2016_community_pc1_", habitat, "_spatial.pdf")),
    spatial_panel,
    width = 6,
    height = 4.8,
    bg = "white"
  )
}
for (enzyme in eea_names) {
  enzyme_rows <- eea_plot_data$enzyme == enzyme
  eea_panel <- (eea_plot + eea_plot_data[enzyme_rows, ]) +
    facet_null() +
    labs(y = paste0(enzyme, " (umol g-1 hr-1)"))
  ggsave(
    file.path(figures_dir, paste0("2016_sediment_community_", tolower(enzyme), ".pdf")),
    eea_panel,
    width = 6,
    height = 4.8,
    bg = "white"
  )
}

message("Exploratory network and EEA figures complete.")
