# Distance decay of sediment microbial and molecular composition on the 33 matched sites.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
library(igraph)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
dispersion_dir <- here("results", "fticr_2016", "matched_sediment_dispersion")
results_dir <- here("results", "fticr_2016", "matched_sediment_distance_decay")
figures_dir <- here("figures", "supplement", "sediment_comparison")
n_permutations <- 9999
permutation_seed <- 20261116
significance_threshold <- 0.05
# Tolerance for matching graph distances to the centerline's own outlet distances.
distance_tolerance_m <- 0.01
# Figures are built at printed size: each panel has a fixed size and the canvas
# is measured around it.
panel_width <- 2.6
panel_height <- 2.0
base_size <- 9
point_size <- 0.8
point_alpha <- 0.5
# Views of the dissimilarity between two sites. The primary views are the
# standard ones for each data type. Binary Jaccard dissimilarity divides into
# turnover, where features or ASVs are replaced, and nestedness, where the
# poorer site holds a subset of the richer one (Baselga 2010).
view_design <- tribble(
  ~metric, ~data_type, ~role, ~label,
  "microbial_hellinger", "microbial", "primary", "Microbial communities (Hellinger)",
  "molecular_jaccard", "molecular", "primary", "FT-ICR features (Jaccard)",
  "microbial_turnover", "microbial", "turnover", "Microbial communities (Jaccard turnover)",
  "molecular_turnover", "molecular", "turnover", "FT-ICR features (Jaccard turnover)",
  "microbial_nestedness", "microbial", "nestedness", "Microbial communities (Jaccard nestedness)",
  "molecular_nestedness", "molecular", "nestedness", "FT-ICR features (Jaccard nestedness)")
# Complete-data matrices saved by 04_sediment_convergence.R. Nestedness is
# binary Jaccard minus its turnover component.
saved_metrics <- c("microbial_hellinger", "molecular_jaccard", "microbial_jaccard",
  "microbial_turnover", "molecular_turnover")
# Two spatial separations: along the channel and in a straight line. Differences
# in network position (drainage area) are tested in 04_sediment_convergence.R.
predictor_design <- tribble(
  ~predictor, ~label,
  "watercourse_km", "Watercourse distance (km)",
  "straight_line_km", "Straight-line distance (km)")
connection_levels <- c("Flow-connected", "Not flow-connected")
connection_colors <- c("Flow-connected" = "#000000", "Not flow-connected" = "#999999")
plot_theme <- theme_bw(base_size = base_size) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), legend.position = "bottom",
    panel.widths = grid::unit(panel_width, "in"), panel.heights = grid::unit(panel_height, "in"),
    plot.tag = element_text(size = 11, face = "bold"))
# Save a figure whose panels have fixed sizes on a canvas measured around them.
# Measuring text needs an open graphics device, so a null device is opened.
save_at_panel_size <- function(plot, path) {
  pdf(NULL)
  if (inherits(plot, "patchwork")) {
    grob <- patchworkGrob(plot)
  } else {
    grob <- ggplotGrob(plot)
  }
  width <- grid::convertWidth(sum(grob$widths), "in", valueOnly = TRUE)
  height <- grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE)
  dev.off()
  ggsave(path, plot, width = width, height = height, bg = "white")
}

# 2. Read the matched sites and their saved dissimilarities -----------------------
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
site <- read_csv(file.path(dispersion_dir, "matched_site_manifest.csv"),
  col_types = cols(site_code = col_character()))
stopifnot(nrow(site) == 33, !anyDuplicated(site$site_code))
dissimilarities <- list()
for (metric in saved_metrics) {
  saved <- read_csv(file.path(dispersion_dir, paste0(metric, "_distance_matrix.csv")),
    col_types = cols(site_code = col_character(), .default = col_double()))
  stopifnot(identical(saved$site_code, site$site_code), identical(names(saved)[-1], site$site_code))
  current <- as.matrix(saved[, -1])
  rownames(current) <- site$site_code
  stopifnot(all(is.finite(current)), max(abs(current - t(current))) < 1e-12, all(diag(current) == 0))
  dissimilarities[[metric]] <- current
}
dissimilarities$microbial_nestedness <- dissimilarities$microbial_jaccard - dissimilarities$microbial_turnover
dissimilarities$molecular_nestedness <- dissimilarities$molecular_jaccard - dissimilarities$molecular_turnover
stopifnot(all(view_design$metric %in% names(dissimilarities)),
  min(dissimilarities$microbial_nestedness) > -1e-12, min(dissimilarities$molecular_nestedness) > -1e-12)

# 3. Build the stream network from the surveyed centerline ------------------------
# Each centerline row is one point; D is its channel distance to the basin outlet.
# The raw header has a trailing comma but data rows have seven fields.
geometry <- read_csv(here("data", "StreamCenterline.csv"), skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
geometry$node <- seq_len(nrow(geometry))
# Rows run down each segment, so consecutive rows of a segment are joined by the
# change in their outlet distance.
same_segment <- geometry$segment[-1] == geometry$segment[-nrow(geometry)]
channel_edges <- tibble(from = geometry$node[-nrow(geometry)][same_segment],
  to = geometry$node[-1][same_segment])
channel_edges$weight <- abs(geometry$outlet_distance[channel_edges$from] -
  geometry$outlet_distance[channel_edges$to])
# Segments meet at confluences, where their end points have identical coordinates.
shared_points <- split(geometry$node, paste(geometry$x, geometry$y))
shared_points <- shared_points[lengths(shared_points) > 1]
confluence_rows <- list()
for (i in seq_along(shared_points)) {
  nodes <- shared_points[[i]]
  confluence_rows[[i]] <- tibble(from = nodes[1], to = nodes[-1], weight = 0)
}
confluence_edges <- bind_rows(confluence_rows)
stream_graph <- graph_from_data_frame(bind_rows(channel_edges, confluence_edges),
  directed = FALSE, vertices = tibble(name = geometry$node))
# Check the graph against the source: along the connected network, the graph
# distance from every point to the outlet must equal the centerline's own D.
outlet_node <- which.min(geometry$outlet_distance)
graph_outlet_distance <- as.numeric(distances(stream_graph, v = outlet_node,
  weights = E(stream_graph)$weight))
membership <- components(stream_graph)$membership
connected <- membership == membership[outlet_node]
graph_checks <- tibble(
  centerline_points = nrow(geometry),
  segments = n_distinct(geometry$segment),
  confluences = length(shared_points),
  components = max(membership),
  segments_not_reaching_outlet = paste(sort(unique(geometry$segment[!connected])), collapse = ";"),
  points_not_reaching_outlet = sum(!connected),
  max_outlet_distance_error_m = max(abs(graph_outlet_distance[connected] -
    geometry$outlet_distance[connected])))
stopifnot(graph_checks$max_outlet_distance_error_m < distance_tolerance_m)

# 4. Measure the separation of every pair of sites ---------------------------------
# Each site is placed at the nearest centerline point of its own segment.
site$centerline_node <- NA_integer_
site$centerline_match_error_m <- NA_real_
for (i in seq_len(nrow(site))) {
  candidates <- which(geometry$segment == site$site_segment_number[i])
  stopifnot(length(candidates) > 0)
  offset <- sqrt((geometry$x[candidates] - site$site_utm_x_m[i])^2 +
    (geometry$y[candidates] - site$site_utm_y_m[i])^2)
  site$centerline_node[i] <- candidates[which.min(offset)]
  site$centerline_match_error_m[i] <- min(offset)
}
site_outlet_distance <- geometry$outlet_distance[site$centerline_node]
stopifnot(all(connected[site$centerline_node]), all(site$centerline_match_error_m < 1),
  max(abs(site_outlet_distance - site$basin_outlet_stream_distance_m)) < distance_tolerance_m)
watercourse_m <- distances(stream_graph, v = site$centerline_node, to = site$centerline_node,
  weights = E(stream_graph)$weight)
# Two sites are flow-connected when water passes from one to the other, so the
# channel distance between them equals the difference of their outlet distances.
outlet_difference_m <- abs(outer(site_outlet_distance, site_outlet_distance, "-"))
flow_connected <- abs(watercourse_m - outlet_difference_m) < distance_tolerance_m
separations <- list(
  watercourse_km = watercourse_m / 1000,
  straight_line_km = as.matrix(dist(cbind(site$site_utm_x_m, site$site_utm_y_m))) / 1000)
stopifnot(all(watercourse_m + distance_tolerance_m >= outlet_difference_m),
  all(separations$watercourse_km + 1e-6 >= separations$straight_line_km))
# One row per unordered pair of sites, with every separation and dissimilarity.
lower <- which(lower.tri(watercourse_m), arr.ind = TRUE)
pairs <- tibble(site_a = site$site_code[lower[, "col"]], site_b = site$site_code[lower[, "row"]],
  group_a = site$network_group[lower[, "col"]], group_b = site$network_group[lower[, "row"]],
  flow_connected = flow_connected[lower])
for (predictor in predictor_design$predictor) {
  pairs[[predictor]] <- separations[[predictor]][lower]
}
for (metric in view_design$metric) {
  pairs[[metric]] <- dissimilarities[[metric]][lower]
}
stopifnot(nrow(pairs) == choose(nrow(site), 2))
write_csv(pairs, file.path(results_dir, "site_pair_separations.csv"))
write_csv(graph_checks %>% mutate(flow_connected_pairs = sum(pairs$flow_connected),
  site_pairs = nrow(pairs), max_site_match_error_m = max(site$centerline_match_error_m)),
  file.path(results_dir, "stream_graph_checks.csv"))

# 5. Test distance decay for each view and separation ------------------------------
# Mantel tests permute whole sites, so every pair involving a site moves together.
# All tests use the same permutations and are two-sided, because turnover and
# nestedness can fall as well as rise with separation. The slope is the
# least-squares line through the site pairs, which is what the Pearson Mantel
# statistic tests.
set.seed(permutation_seed)
permutations <- shuffleSet(nrow(site), nset = n_permutations)
mantel_rows <- list()
for (view in seq_len(nrow(view_design))) {
  metric <- view_design$metric[view]
  for (predictor in predictor_design$predictor) {
    fit <- mantel(as.dist(separations[[predictor]]), as.dist(dissimilarities[[metric]]),
      method = "pearson", permutations = permutations)
    line <- lm(pairs[[metric]] ~ pairs[[predictor]])
    mantel_rows[[length(mantel_rows) + 1]] <- tibble(metric, data_type = view_design$data_type[view],
      role = view_design$role[view], predictor, sites = nrow(site), site_pairs = nrow(pairs),
      mantel_r = fit$statistic, intercept = coef(line)[[1]], slope = coef(line)[[2]],
      p_value = (sum(abs(fit$perm) >= abs(fit$statistic)) + 1) / (length(fit$perm) + 1))
  }
}
mantel_tests <- bind_rows(mantel_rows) %>% group_by(role) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
write_csv(mantel_tests, file.path(results_dir, "mantel_tests.csv"))

# 6. Test decay separately for flow-connected and unconnected pairs ----------------
# Water carries material between flow-connected sites, so composition may change
# with channel distance along a flow path and not between separate branches.
# The correlation with watercourse distance is computed within each set of
# pairs. Its null distribution permutes whole sites in the dissimilarity matrix
# while the channel distances and flow connections stay in place. The last test
# compares the two slopes.
connected <- pairs$flow_connected
distance <- pairs$watercourse_km
flow_rows <- list()
for (view in seq_len(nrow(view_design))) {
  metric <- view_design$metric[view]
  current <- dissimilarities[[metric]]
  observed <- current[lower]
  connected_line <- lm(observed[connected] ~ distance[connected])
  unconnected_line <- lm(observed[!connected] ~ distance[!connected])
  observed_r <- c(cor(observed[connected], distance[connected]), cor(observed[!connected], distance[!connected]))
  observed_difference <- coef(connected_line)[[2]] - coef(unconnected_line)[[2]]
  null_r <- matrix(NA_real_, n_permutations, 2)
  null_difference <- rep(NA_real_, n_permutations)
  for (k in seq_len(n_permutations)) {
    order <- permutations[k, ]
    shuffled <- current[order, order][lower]
    null_r[k, 1] <- cor(shuffled[connected], distance[connected])
    null_r[k, 2] <- cor(shuffled[!connected], distance[!connected])
    null_difference[k] <- cov(shuffled[connected], distance[connected]) / var(distance[connected]) -
      cov(shuffled[!connected], distance[!connected]) / var(distance[!connected])
  }
  flow_rows[[view]] <- tibble(metric, data_type = view_design$data_type[view], role = view_design$role[view],
    pair_set = c(connection_levels, "Slope difference"),
    site_pairs = c(sum(connected), sum(!connected), length(connected)),
    pearson_r = c(observed_r, NA_real_),
    intercept = c(coef(connected_line)[[1]], coef(unconnected_line)[[1]], NA_real_),
    slope_per_km = c(coef(connected_line)[[2]], coef(unconnected_line)[[2]], observed_difference),
    p_value = c((sum(abs(null_r[, 1]) >= abs(observed_r[1])) + 1) / (n_permutations + 1),
      (sum(abs(null_r[, 2]) >= abs(observed_r[2])) + 1) / (n_permutations + 1),
      (sum(abs(null_difference) >= abs(observed_difference)) + 1) / (n_permutations + 1)))
}
# Slopes within pair sets and slope differences are separate BH families per role.
flow_tests <- bind_rows(flow_rows) %>% group_by(role, pair_set == "Slope difference") %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  select(-`pair_set == "Slope difference"`)
write_csv(flow_tests, file.path(results_dir, "flow_connection_tests.csv"))

# 7. Compare microbial and molecular decay on the same site pairs ------------------
# The two data types use different dissimilarity scales, so each is standardized
# across the site pairs. Their difference for each pair is then related to
# separation. A positive association means molecular dissimilarity rises faster
# with separation than microbial dissimilarity. The test is two-sided.
contrast_rows <- list()
for (role in unique(view_design$role)) {
  microbial_metric <- view_design$metric[view_design$role == role & view_design$data_type == "microbial"]
  molecular_metric <- view_design$metric[view_design$role == role & view_design$data_type == "molecular"]
  microbial <- dissimilarities[[microbial_metric]]
  molecular <- dissimilarities[[molecular_metric]]
  difference <- (molecular - mean(molecular[lower])) / sd(molecular[lower]) -
    (microbial - mean(microbial[lower])) / sd(microbial[lower])
  diag(difference) <- 0
  for (predictor in predictor_design$predictor) {
    fit <- mantel(as.dist(separations[[predictor]]), as.dist(difference),
      method = "pearson", permutations = permutations)
    contrast_rows[[length(contrast_rows) + 1]] <- tibble(role, molecular_metric, microbial_metric,
      predictor, mantel_r = fit$statistic,
      p_value = (sum(abs(fit$perm) >= abs(fit$statistic)) + 1) / (length(fit$perm) + 1))
  }
}
contrast_tests <- bind_rows(contrast_rows) %>% group_by(role) %>%
  mutate(p_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
write_csv(contrast_tests, file.path(results_dir, "paired_decay_contrasts.csv"))

# 8. Plot every site pair against watercourse distance -----------------------------
# A line is drawn for a set of pairs only where its adjusted test passes. The
# other separations are reported in the Mantel table.
pairs$connection <- factor(if_else(pairs$flow_connected, connection_levels[1], connection_levels[2]),
  levels = connection_levels)
panels <- list()
for (view in seq_len(nrow(view_design))) {
  metric <- view_design$metric[view]
  panel <- tibble(separation = pairs$watercourse_km, dissimilarity = pairs[[metric]],
    connection = pairs$connection, metric)
  # Unconnected pairs are drawn first so the fewer flow-connected pairs stay visible.
  p <- ggplot(panel %>% arrange(desc(connection)), aes(separation, dissimilarity, color = connection)) +
    geom_point(size = point_size, alpha = point_alpha) +
    facet_wrap(~metric, labeller = as_labeller(setNames(view_design$label, view_design$metric))) +
    scale_color_manual(values = connection_colors, drop = FALSE) +
    labs(x = "Watercourse distance (km)", y = "Dissimilarity between sites", color = NULL) +
    plot_theme
  for (pair_set in connection_levels) {
    tested <- flow_tests %>% filter(.data$metric == .env$metric, .data$pair_set == .env$pair_set)
    if (tested$p_bh <= significance_threshold) {
      p <- p + geom_abline(intercept = tested$intercept, slope = tested$slope_per_km,
        color = connection_colors[[pair_set]], linewidth = 0.7)
    }
  }
  panels[[metric]] <- p
  save_at_panel_size(p, file.path(figures_dir, paste0("2016_matched_sediment_distance_decay_", metric, ".pdf")))
}
for (role in unique(view_design$role)) {
  metrics <- view_design$metric[view_design$role == role]
  combined <- wrap_plots(panels[metrics], ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
  save_at_panel_size(combined, file.path(figures_dir, paste0("2016_matched_sediment_distance_decay_", role, ".pdf")))
}

writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
print(graph_checks)
print(mantel_tests, n = Inf)
print(flow_tests, n = Inf)
print(contrast_tests, n = Inf)
