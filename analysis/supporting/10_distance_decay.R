# Distance decay of sediment microbial and molecular composition on the 33 matched sites.
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
dispersion_dir <- here("results", "fticr_2016", "matched_sediment_dispersion")
results_dir <- here("results", "fticr_2016", "matched_sediment_distance_decay")
figures_dir <- here("figures", "supplement", "sediment_comparison")
n_permutations <- 9999
permutation_seed <- 20261116
significance_threshold <- 0.05
# Figures are built at printed size: each panel has a fixed size and the canvas
# is measured around it.
panel_width <- 2.6
panel_height <- 2.0
base_size <- 9
point_size <- 0.8
point_alpha <- 0.5
# The standard dissimilarity for each data type. Its replacement and
# richness-difference components are tested in 04_sediment_convergence.R.
view_design <- tribble(
  ~metric, ~data_type, ~role, ~label,
  "microbial_hellinger", "microbial", "primary", "Microbial communities (Hellinger)",
  "molecular_jaccard", "molecular", "primary", "FT-ICR features (Jaccard)")
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

# 2. Read the matched sites, their dissimilarities, and their separations ----------
# 04_sediment_convergence.R saves the complete-data matrices and one row per pair
# of sites with its channel and straight-line distances and flow connection.
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
site <- read_csv(file.path(dispersion_dir, "matched_site_manifest.csv"),
  col_types = cols(site_code = col_character()))
stopifnot(nrow(site) == 33, !anyDuplicated(site$site_code))
dissimilarities <- list()
for (metric in view_design$metric) {
  saved <- read_csv(file.path(dispersion_dir, paste0(metric, "_distance_matrix.csv")),
    col_types = cols(site_code = col_character(), .default = col_double()))
  stopifnot(identical(saved$site_code, site$site_code), identical(names(saved)[-1], site$site_code))
  current <- as.matrix(saved[, -1])
  rownames(current) <- site$site_code
  stopifnot(all(is.finite(current)), max(abs(current - t(current))) < 1e-12, all(diag(current) == 0))
  dissimilarities[[metric]] <- current
}
pairs <- read_csv(file.path(dispersion_dir, "site_pair_components.csv"),
  col_types = cols(site_a = col_character(), site_b = col_character(), downstream_site = col_character())) %>%
  select(site_a, site_b, group_a, group_b, flow_connected, watercourse_km, straight_line_km)
# Rebuild full separation matrices from the pair rows, in the site order of the
# dissimilarity matrices.
lower <- cbind(row = match(pairs$site_b, site$site_code), col = match(pairs$site_a, site$site_code))
stopifnot(nrow(pairs) == choose(nrow(site), 2), !anyNA(lower), all(lower[, "row"] > lower[, "col"]))
separations <- list()
for (predictor in predictor_design$predictor) {
  current <- matrix(0, nrow(site), nrow(site), dimnames = list(site$site_code, site$site_code))
  current[lower] <- pairs[[predictor]]
  separations[[predictor]] <- current + t(current)
}
for (metric in view_design$metric) {
  pairs[[metric]] <- dissimilarities[[metric]][lower]
}
write_csv(pairs, file.path(results_dir, "site_pair_separations.csv"))

# 5. Test distance decay for each view and separation ------------------------------
# Mantel tests permute whole sites, so every pair involving a site moves together.
# All tests use the same permutations and are two-sided. The slope is the
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
combined <- wrap_plots(panels[view_design$metric], ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
save_at_panel_size(combined, file.path(figures_dir, "2016_matched_sediment_distance_decay_primary.pdf"))

writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(results_dir, "session_info.txt"))
print(mantel_tests, n = Inf)
print(flow_tests, n = Inf)
print(contrast_tests, n = Inf)
