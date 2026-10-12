# Descriptive catchment maps and network-position summaries of completed iCAMP.
# Pairwise assignments share samples; resampling ranges describe sensitivity,
# not confidence intervals or independent tests of network-position effects.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1. Settings and plotting theme ------------------------------------------------
input_dir <- here("results", "icamp_2016", "sediment_castor")
output_dir <- file.path(input_dir, "network_figures")
figures_dir <- here("figures", "supplement")
group_levels <- c("Headwater", "Intermediate", "Mainstem")
headwater_orders <- c(1, 2)
intermediate_orders <- c(3, 4)
mainstem_order <- 5
n_balanced_draws <- 1000L
segments_per_group <- 5L
base_seed <- 20261012L
centerline_stride <- 10L
map_width <- 6.8
map_height <- 5.6
comparison_width <- 5.4
comparison_height <- 4.5
process_design <- tribble(
  ~process, ~code, ~label,
  "Heterogeneous.Selection", "HeS", "Heterogeneous selection",
  "Homogeneous.Selection", "HoS", "Homogeneous selection",
  "Dispersal.Limitation", "DL", "Dispersal limitation",
  "Homogenizing.Dispersal", "HD", "Homogenizing dispersal",
  "Drift.and.Others", "DR", "Drift and other processes"
)
group_colors <- c(Headwater = "#009E73", Intermediate = "#E69F00", Mainstem = "#0072B2")
process_colors <- c(HeS = "#D55E00", HoS = "#009E73", DL = "#0072B2",
  HD = "#CC79A7", DR = "#999999")
plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 10),
    plot.margin = margin(8, 20, 8, 8)
  )

# 2. Check sample identities and completed pairwise assignments ------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c(
  file.path(input_dir, "pairwise_process_importance.csv"),
  file.path(input_dir, "sediment_samples.csv"),
  file.path(input_dir, "icamp_settings.csv"),
  here("data", "derived", "soil_porewater_sediment_2016", "sample_metadata.csv"),
  here("data", "StreamCenterline.csv")
)
checks <- tibble(
  path = sub(paste0(here(), "/"), "", input_paths, fixed = TRUE),
  initial_md5 = unname(tools::md5sum(input_paths))
)
pairwise <- read_csv(input_paths[1], show_col_types = FALSE)
samples <- read_csv(input_paths[2], col_types = cols(site_code = col_character()))
settings <- read_csv(input_paths[3], show_col_types = FALSE)
metadata <- read_csv(input_paths[4], col_types = cols(site_code = col_character())) %>%
  filter(habitat == "sediment")
stopifnot(
  nrow(samples) == 34L,
  !anyDuplicated(samples$sample_id),
  !anyDuplicated(samples$site_code),
  setequal(samples$sample_id, metadata$sample_id),
  settings$value[settings$setting == "randomizations"] == "1000",
  nrow(pairwise) == choose(nrow(samples), 2),
  all(pairwise$sample1 %in% samples$sample_id),
  all(pairwise$sample2 %in% samples$sample_id),
  all(pairwise$sample1 != pairwise$sample2)
)
metadata <- metadata[match(samples$sample_id, metadata$sample_id), ]
stopifnot(
  identical(samples$site_code, metadata$site_code),
  all(samples$site_stream_order == metadata$site_stream_order),
  max(abs(samples$site_utm_x_m - metadata$site_utm_x_m)) < 1e-8,
  max(abs(samples$site_utm_y_m - metadata$site_utm_y_m)) < 1e-8
)
site <- samples %>%
  mutate(
    site_segment_number = metadata$site_segment_number,
    network_group = case_when(
      site_stream_order %in% headwater_orders ~ "Headwater",
      site_stream_order %in% intermediate_orders ~ "Intermediate",
      site_stream_order == mainstem_order ~ "Mainstem"
    ),
    group_index = match(network_group, group_levels)
  )
stopifnot(!anyNA(site), !anyDuplicated(paste(pmin(pairwise$sample1, pairwise$sample2),
  pmax(pairwise$sample1, pairwise$sample2), sep = "|")))
fractions <- as.matrix(pairwise[, process_design$process])
stopifnot(all(is.finite(fractions)), all(fractions >= 0 & fractions <= 1),
  max(abs(rowSums(fractions) - 1)) < 1e-10)
sizes <- site %>% group_by(network_group) %>%
  summarize(sites = n(), segments = n_distinct(site_segment_number), .groups = "drop")
stopifnot(all(sizes$segments >= segments_per_group))
segment_groups <- site %>% distinct(site_segment_number, network_group)
stopifnot(!anyDuplicated(segment_groups$site_segment_number))

# 3. Describe each site relative to other sampled stream segments -----------------
# Average over sites within each partner segment, then weight partner segments
# equally. Excluding the focal segment prevents nearby replicates from defining
# a site's mapped value. These remain pairwise summaries, not local process rates.
first <- match(pairwise$sample1, site$sample_id)
second <- match(pairwise$sample2, site$sample_id)
pairwise <- pairwise %>%
  mutate(
    segment1 = site$site_segment_number[first],
    segment2 = site$site_segment_number[second],
    group1 = site$group_index[first],
    group2 = site$group_index[second],
    same_segment = segment1 == segment2,
    segment_low = pmin(segment1, segment2),
    segment_high = pmax(segment1, segment2),
    group_low = pmin(group1, group2),
    group_high = pmax(group1, group2)
  )
long_pairs <- pairwise %>% pivot_longer(all_of(process_design$process),
  names_to = "process", values_to = "importance")
incident <- bind_rows(
  long_pairs %>% transmute(sample_id = sample1, partner_id = sample2,
    partner_segment = segment2, same_segment, process, importance),
  long_pairs %>% transmute(sample_id = sample2, partner_id = sample1,
    partner_segment = segment1, same_segment, process, importance)
)
site_profiles <- incident %>% filter(!same_segment) %>%
  group_by(sample_id, partner_segment, process) %>%
  summarize(segment_mean = mean(importance), partner_sites = n(), .groups = "drop") %>%
  group_by(sample_id, process) %>%
  summarize(segment_balanced_mean = mean(segment_mean),
    partner_segments = n(), partner_sites = sum(partner_sites), .groups = "drop") %>%
  left_join(incident %>% group_by(sample_id, process) %>%
    summarize(all_partner_mean = mean(importance), .groups = "drop"),
    by = c("sample_id", "process")) %>%
  left_join(site, by = "sample_id")
stopifnot(nrow(site_profiles) == nrow(site) * nrow(process_design),
  all(site_profiles$partner_segments == nrow(segment_groups) - 1L))
site_sums <- site_profiles %>% group_by(sample_id) %>%
  summarize(total = sum(segment_balanced_mean), .groups = "drop")
stopifnot(max(abs(site_sums$total - 1)) < 1e-10)
write_csv(site_profiles, file.path(output_dir, "site_process_importance.csv"))

# 4. Give each distinct segment pair equal weight in each network comparison -----
# Within-segment pairs are retained in the audit but excluded from comparisons.
# Six cells distinguish within-group turnover from turnover between positions.
segment_pairs <- long_pairs %>% filter(!same_segment) %>%
  group_by(group_low, group_high, segment_low, segment_high, process) %>%
  summarize(importance = mean(importance), sample_pairs = n(), .groups = "drop")
group_profiles <- segment_pairs %>% group_by(group_low, group_high, process) %>%
  summarize(segment_balanced_mean = mean(importance),
    segment_pairs = n(), sample_pairs = sum(sample_pairs), .groups = "drop") %>%
  left_join(long_pairs %>% filter(!same_segment) %>%
    group_by(group_low, group_high, process) %>%
    summarize(sample_pair_mean = mean(importance), .groups = "drop"),
    by = c("group_low", "group_high", "process")) %>%
  mutate(group1 = group_levels[group_low], group2 = group_levels[group_high])
write_csv(group_profiles, file.path(output_dir, "network_pair_process_importance.csv"))
pair_audit <- pairwise %>%
  transmute(sample1, sample2, site1 = site$site_code[first], site2 = site$site_code[second],
    segment1, segment2, group1 = group_levels[group1], group2 = group_levels[group2],
    same_segment, included_network_comparison = !same_segment)
write_csv(pair_audit, file.path(output_dir, "pair_structure_audit.csv"))

# 5. Balance network groups by selecting one site in each of five segments -------
# The same 15 sites supply every process and every group comparison in a draw.
# This samples the completed assignments; it does not rerun the iCAMP null model.
draw_design <- tibble(run = seq_len(n_balanced_draws), seed = base_seed + run)
draw_rows <- vector("list", n_balanced_draws)
selection_rows <- vector("list", n_balanced_draws)
start_time <- proc.time()[["elapsed"]]
for (run in draw_design$run) {
  set.seed(draw_design$seed[run])
  selected <- integer(0)
  for (group in seq_along(group_levels)) {
    candidates <- which(site$group_index == group)
    segments <- unique(site$site_segment_number[candidates])
    selected_segments <- segments[sample.int(length(segments), segments_per_group)]
    for (segment in selected_segments) {
      indices <- candidates[site$site_segment_number[candidates] == segment]
      selected <- c(selected, indices[sample.int(length(indices), 1)])
    }
  }
  selected_ids <- site$sample_id[selected]
  selection_rows[[run]] <- site[selected, ] %>%
    transmute(run, seed = draw_design$seed[run], sample_id, site_code,
      network_group, segment = site_segment_number)
  draw_rows[[run]] <- long_pairs %>%
    filter(sample1 %in% selected_ids, sample2 %in% selected_ids) %>%
    group_by(group_low, group_high, process) %>%
    summarize(importance = mean(importance), sample_pairs = n(), .groups = "drop") %>%
    mutate(run, seed = draw_design$seed[run])
  if (run == 1L) {
    message("Representative balanced draw: ",
      round(proc.time()[["elapsed"]] - start_time, 3), " s")
  }
}
draws <- bind_rows(draw_rows)
selections <- bind_rows(selection_rows)
selection_audit <- selections %>% group_by(run, network_group) %>%
  summarize(sites = n(), segments = n_distinct(segment), .groups = "drop")
stopifnot(nrow(draws) == n_balanced_draws * 6L * nrow(process_design),
  all(selection_audit$sites == segments_per_group),
  all(selection_audit$segments == segments_per_group),
  all(draws$sample_pairs[draws$group_low == draws$group_high] == choose(segments_per_group, 2)),
  all(draws$sample_pairs[draws$group_low != draws$group_high] == segments_per_group^2))
balanced_summary <- draws %>% group_by(group_low, group_high, process) %>%
  summarize(draws = n(), median = median(importance),
    lower = quantile(importance, 0.025), upper = quantile(importance, 0.975),
    .groups = "drop") %>%
  mutate(group1 = group_levels[group_low], group2 = group_levels[group_high])
write_csv(draws, file.path(output_dir, "balanced_segment_draws.csv"))
write_csv(selections, file.path(output_dir, "balanced_site_selections.csv"))
# Keep the draw-level contrasts so the sensitivity of differences is inspectable.
contrasts <- tibble(group_a = c(1L, 1L, 2L), group_b = c(2L, 3L, 3L))
contrast_rows <- vector("list", nrow(contrasts))
within_draws <- draws %>% filter(group_low == group_high)
for (i in seq_len(nrow(contrasts))) {
  a <- within_draws %>% filter(group_low == contrasts$group_a[i])
  b <- within_draws %>% filter(group_low == contrasts$group_b[i])
  stopifnot(identical(a$run, b$run), identical(a$process, b$process))
  contrast_rows[[i]] <- tibble(run = a$run, seed = a$seed, process = a$process,
    contrast = paste(group_levels[contrasts$group_a[i]], group_levels[contrasts$group_b[i]], sep = " - "),
    difference = a$importance - b$importance)
}
write_csv(bind_rows(contrast_rows), file.path(output_dir, "balanced_within_group_contrasts.csv"))

# 6. Map segment-balanced pairwise importance at each sampled location ------------
geometry <- read_csv(input_paths[5], skip = 1,
  col_names = c("x", "y", "elevation", "outlet_distance", "upstream_area", "order", "segment"),
  col_types = cols(.default = col_double()))
stopifnot(nrow(problems(geometry)) == 0L,
  all(is.finite(geometry$x)), all(is.finite(geometry$y)))
map_lines <- geometry %>% group_by(segment) %>%
  filter(row_number() %% centerline_stride == 1L | row_number() == n()) %>% ungroup()
maps <- vector("list", nrow(process_design))
comparisons <- vector("list", nrow(process_design))
size_reference <- sizes[match(group_levels, sizes$network_group), ]
group_labels <- paste0(group_levels, "\n", size_reference$sites, " sites\n",
  size_reference$segments, " segments")
for (i in seq_len(nrow(process_design))) {
  current_process <- process_design$process[i]
  mapped <- site_profiles %>% filter(process == current_process)
  map <- ggplot(map_lines, aes(x, y, group = segment)) +
    geom_path(color = "#BBBBBB", linewidth = 0.3) +
    geom_point(data = mapped,
      aes(site_utm_x_m, site_utm_y_m, fill = 100 * segment_balanced_mean),
      inherit.aes = FALSE, shape = 21, size = 3.1, color = "black", stroke = 0.35) +
    # Each process has its own color range to reveal spatial differences;
    # its numerical legend remains in percentages for direct interpretation.
    scale_fill_viridis_c(option = "C",
      labels = function(x) sprintf("%.1f", x), name = "Mean pairwise importance (%)") +
    coord_equal() + labs(x = "UTM easting (m)", y = "UTM northing (m)") +
    facet_wrap(~process, labeller = as_labeller(setNames(process_design$label, process_design$process))) +
    plot_theme
  maps[[i]] <- map
  ggsave(file.path(figures_dir, paste0("2016_sediment_icamp_map_", process_design$code[i], ".pdf")),
    map, width = map_width, height = map_height, bg = "white")

  # Each range is the middle 95% of balanced draws, with their median marked.
  # Hollow diamonds retain the mean using every observed segment pair.
  panel <- balanced_summary %>%
    filter(process == current_process, group_low == group_high) %>%
    mutate(network_group = factor(group1, levels = group_levels))
  full <- group_profiles %>%
    filter(process == current_process, group_low == group_high) %>%
    mutate(network_group = factor(group1, levels = group_levels))
  comparison <- ggplot(panel,
    aes(group_low - 0.06, 100 * median, color = network_group)) +
    geom_linerange(aes(ymin = 100 * lower, ymax = 100 * upper), linewidth = 1) +
    geom_point(aes(shape = "Balanced draws: median"), size = 2.8) +
    geom_point(data = full,
      aes(x = group_low + 0.06, y = 100 * segment_balanced_mean,
        shape = "All segment pairs: mean"),
      fill = "white", size = 2.8, stroke = 0.8) +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_shape_manual(values = c("Balanced draws: median" = 16,
      "All segment pairs: mean" = 23), name = NULL) +
    scale_x_continuous(breaks = seq_along(group_levels), labels = group_labels,
      limits = c(0.5, 3.5)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.04, 0.08))) +
    labs(x = NULL, y = "Within-group pairwise importance (%)") +
    facet_wrap(~process, labeller = as_labeller(setNames(process_design$label, process_design$process))) +
    plot_theme + theme(aspect.ratio = 0.8, axis.text.x = element_text(size = 9))
  comparisons[[i]] <- comparison
  ggsave(file.path(figures_dir, paste0("2016_sediment_icamp_network_", process_design$code[i], ".pdf")),
    comparison, width = comparison_width, height = comparison_height, bg = "white")
}
# The sixth map identifies the network positions used in the comparisons.
position_map <- ggplot(map_lines, aes(x, y, group = segment)) +
  geom_path(color = "#BBBBBB", linewidth = 0.3) +
  geom_point(data = site %>% mutate(process = "Network position"),
    aes(site_utm_x_m, site_utm_y_m, fill = factor(network_group, levels = group_levels)),
    inherit.aes = FALSE, shape = 21, size = 3.1, color = "black", stroke = 0.35) +
  scale_fill_manual(values = group_colors, name = NULL) +
  coord_equal() + labs(x = "UTM easting (m)", y = "UTM northing (m)") +
  facet_wrap(~process) + plot_theme
ggsave(file.path(figures_dir, "2016_sediment_icamp_catchment_maps.pdf"),
  wrap_plots(c(maps, list(position_map)), ncol = 3),
  width = 15.6, height = 10.4, bg = "white")
ggsave(file.path(figures_dir, "2016_sediment_icamp_network_comparison.pdf"),
  wrap_plots(comparisons, ncol = 3, guides = "collect") & theme(legend.position = "bottom"),
  width = 15.6, height = 9, bg = "white")

# A compact composition view keeps all processes on the common 0-100% scale.
composition <- group_profiles %>% filter(group_low == group_high) %>%
  mutate(network_group = factor(group1, levels = group_levels),
    code = factor(process_design$code[match(process, process_design$process)],
      levels = process_design$code))
composition_plot <- ggplot(composition,
  aes(network_group, 100 * segment_balanced_mean, fill = code)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = process_colors,
    labels = setNames(process_design$label, process_design$code), name = NULL) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Within-group pairwise importance (%)") + plot_theme +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(size = 9))
ggsave(file.path(figures_dir, "2016_sediment_icamp_network_composition.pdf"),
  composition_plot, width = 7.2, height = 5.5, bg = "white")

# 7. Distinguish within-position and between-position comparisons -----------------
heatmap <- group_profiles %>%
  mutate(process = factor(process, levels = process_design$process),
    group1 = factor(group1, levels = group_levels),
    group2 = factor(group2, levels = rev(group_levels)),
    label = sprintf("%.1f", 100 * segment_balanced_mean))
pair_plot <- ggplot(heatmap, aes(group1, group2, fill = 100 * segment_balanced_mean)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label, color = segment_balanced_mean < 0.5), size = 3.6) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  scale_fill_viridis_c(option = "C", limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100), name = "Mean pairwise importance (%)") +
  facet_wrap(~process, ncol = 3,
    labeller = as_labeller(setNames(process_design$label, process_design$process))) +
  coord_equal() + labs(x = NULL, y = NULL) + plot_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(figures_dir, "2016_sediment_icamp_network_pair_comparison.pdf"),
  pair_plot, width = 13, height = 8.5, bg = "white")

# 8. Show every within-position sample pair with equal weight ---------------------
# Include same-segment pairs here, marking them separately. The comparison is
# 171 headwater, 45 intermediate, and 10 mainstem pairs, each using all sites.
# Individual pair rows share sites and are not independent replicates.
within_pairs <- long_pairs %>% filter(group1 == group2) %>%
  mutate(network_group = factor(group_levels[group1], levels = group_levels),
    pair_type = if_else(same_segment, "Same segment", "Different segments"))
within_summary <- within_pairs %>% group_by(network_group, process) %>%
  summarize(sites = n_distinct(c(sample1, sample2)), sample_pairs = n(),
    same_segment_pairs = sum(same_segment),
    mean = mean(importance), median = median(importance),
    lower_quartile = quantile(importance, 0.25),
    upper_quartile = quantile(importance, 0.75),
    minimum = min(importance), maximum = max(importance), .groups = "drop")
within_sizes <- within_summary %>% distinct(network_group, sites, sample_pairs)
stopifnot(identical(within_sizes$sample_pairs, c(171L, 45L, 10L)),
  all(within_sizes$sample_pairs == choose(within_sizes$sites, 2)))
within_labels <- paste0(group_levels, "\n", within_sizes$sites, " sites\n",
  within_sizes$sample_pairs, " pairs")
write_csv(within_summary, file.path(output_dir, "all_within_group_pair_summary.csv"))
all_pair_panels <- vector("list", nrow(process_design))
for (i in seq_len(nrow(process_design))) {
  current_process <- process_design$process[i]
  panel <- within_pairs %>% filter(process == current_process)
  summary <- within_summary %>% filter(process == current_process)
  pair_distribution <- ggplot(panel,
    aes(network_group, 100 * importance, color = network_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.6) +
    geom_point(aes(shape = pair_type),
      position = position_jitter(width = 0.12, height = 0, seed = base_seed),
      alpha = 0.5, size = 1.4) +
    geom_point(data = summary, aes(y = 100 * mean, shape = "Group mean"),
      fill = "white", color = "black", size = 2.8, stroke = 0.8) +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_shape_manual(values = c("Different segments" = 16,
      "Same segment" = 4, "Group mean" = 23), name = NULL) +
    guides(shape = guide_legend(override.aes = list(color = "black", alpha = 1))) +
    scale_x_discrete(labels = within_labels) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.06))) +
    labs(x = NULL, y = "Pairwise importance (%)") +
    facet_wrap(~process,
      labeller = as_labeller(setNames(process_design$label, process_design$process))) +
    plot_theme + theme(aspect.ratio = 0.8, axis.text.x = element_text(size = 9))
  all_pair_panels[[i]] <- pair_distribution
  ggsave(file.path(figures_dir, paste0("2016_sediment_icamp_all_within_pairs_",
    process_design$code[i], ".pdf")), pair_distribution,
    width = comparison_width, height = comparison_height, bg = "white")
}
ggsave(file.path(figures_dir, "2016_sediment_icamp_all_within_pairs.pdf"),
  wrap_plots(all_pair_panels, ncol = 3, guides = "collect") & theme(legend.position = "bottom"),
  width = 15.6, height = 9, bg = "white")
all_pair_composition <- within_summary %>%
  mutate(code = factor(process_design$code[match(process, process_design$process)],
    levels = process_design$code))
all_pair_composition_plot <- ggplot(all_pair_composition,
  aes(network_group, 100 * mean, fill = code)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = process_colors,
    labels = setNames(process_design$label, process_design$code), name = NULL) +
  scale_x_discrete(labels = within_labels) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean within-group pairwise importance (%)") + plot_theme +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(size = 9))
ggsave(file.path(figures_dir, "2016_sediment_icamp_all_within_pair_composition.pdf"),
  all_pair_composition_plot, width = 7.2, height = 5.5, bg = "white")

# 9. Preserve completed inputs and record the plotting environment ----------------
checks$final_md5 <- unname(tools::md5sum(input_paths))
stopifnot(identical(checks$initial_md5, checks$final_md5))
write_csv(checks, file.path(output_dir, "input_checksums.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
  file.path(output_dir, "session_info.txt"))
print(sizes)
print(group_profiles %>% filter(group_low == group_high) %>%
  select(group1, process, segment_balanced_mean, segment_pairs, sample_pairs))
message("Saved descriptive iCAMP catchment maps and network-position comparisons.")
if (interactive()) {
  print(wrap_plots(c(maps, list(position_map)), ncol = 3))
  print(composition_plot)
  print(wrap_plots(comparisons, ncol = 3, guides = "collect"))
  print(pair_plot)
  print(wrap_plots(all_pair_panels, ncol = 3, guides = "collect"))
  print(all_pair_composition_plot)
}
