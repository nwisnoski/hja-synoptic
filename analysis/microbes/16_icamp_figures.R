# Describe the completed catchment iCAMP run without repeating randomization.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

# 1. Settings and plotting theme -----------------------------------------------
input_dir <- here("results", "icamp_2016", "catchment_castor")
output_dir <- file.path(input_dir, "diagnostics")
figures_dir <- here("figures")
tolerance <- 1e-10
score_breaks <- seq(-1, 1, by = 0.05)
process_columns <- c("Heterogeneous.Selection", "Homogeneous.Selection",
  "Dispersal.Limitation", "Homogenizing.Dispersal", "Drift.and.Others")
process_labels <- c("Heterogeneous selection", "Homogeneous selection",
  "Dispersal limitation", "Homogenizing dispersal", "Drift/other")
process_colors <- c("#D55E00", "#E69F00", "#0072B2", "#56B4E9", "#999999")
habitat_order <- c("planktonic / planktonic", "hyporheic / hyporheic",
  "sediment / sediment", "soil / soil", "hyporheic / planktonic",
  "planktonic / sediment", "hyporheic / sediment", "planktonic / soil",
  "hyporheic / soil", "sediment / soil")
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
    legend.position = "bottom")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Read and check the saved result and exported pairs -------------------------
result <- readRDS(file.path(input_dir, "icamp_result.rds"))
pairs <- read_csv(file.path(input_dir, "pairwise_process_importance.csv"),
  show_col_types = FALSE)
samples <- read_csv(file.path(input_dir, "catchment_samples.csv"),
  show_col_types = FALSE)
community <- result$detail$comm
phylo_scores <- as.matrix(result$detail$SigbMPDi[, -(1:2)])
bray_scores <- as.matrix(result$detail$SigBCa[, -(1:2)])
weights <- as.matrix(result$detail$bin.weight[, -(1:2)])
cutoff <- result$detail$setting$conf.cut
stopifnot(result$detail$setting$sig.index == "Confidence")
stopifnot(nrow(pairs) == choose(nrow(community), 2))
stopifnot(isTRUE(all.equal(as.data.frame(pairs[names(result$CbMPDiCBraya)]),
  result$CbMPDiCBraya, check.attributes = FALSE, tolerance = tolerance)))
stopifnot(all(is.finite(phylo_scores)), all(is.finite(bray_scores)))
stopifnot(all(is.finite(weights)), all(weights >= 0))
stopifnot(max(abs(rowSums(weights) - 1)) < tolerance)
stopifnot(max(abs(rowSums(pairs[process_columns]) - 1)) < tolerance)
for (component in c("SigbMPDi", "SigBCa", "bin.weight")) {
  stopifnot(identical(unname(as.matrix(result$detail[[component]][1:2])),
    unname(as.matrix(result$CbMPDiCBraya[1:2]))))
}
stopifnot(identical(sub(".*[.]", "", colnames(phylo_scores)), colnames(weights)))
stopifnot(identical(sub(".*[.]", "", colnames(bray_scores)), colnames(weights)))

# 3. Mean abundance-weighted process fractions by habitat pair ------------------
habitat_summary <- pairs %>%
  group_by(habitat_pair, comparison) %>%
  summarize(n_pairs = n(), across(all_of(process_columns), mean), .groups = "drop")
write_csv(habitat_summary, file.path(output_dir, "habitat_process_summary.csv"))
process_plot_data <- habitat_summary %>%
  pivot_longer(all_of(process_columns), names_to = "process", values_to = "fraction") %>%
  mutate(habitat_pair = factor(habitat_pair, levels = rev(habitat_order)),
    process = factor(process, levels = process_columns, labels = process_labels))
process_plot <- ggplot(process_plot_data, aes(100 * fraction, habitat_pair, fill = process)) +
  geom_col(width = 0.75, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = setNames(process_colors, process_labels)) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  labs(x = "Mean abundance-weighted process fraction (%)", y = NULL, fill = NULL) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  plot_theme
ggsave(file.path(figures_dir, "2016_icamp_catchment_process_fractions.pdf"),
  process_plot, width = 8, height = 6.4, bg = "white")

# 4. Weighted distributions of the saved bin-level null-comparison scores ------
# Summarize matrices directly; avoid expanding millions of bin-pair rows.
score_rows <- list()
row_index <- 1
for (metric in c("phylogenetic", "compositional")) {
  scores <- phylo_scores
  if (metric == "compositional") scores <- bray_scores
  for (habitat_pair in habitat_order) {
    keep <- pairs$habitat_pair == habitat_pair
    values <- as.vector(scores[keep, , drop = FALSE])
    bin_weights <- as.vector(weights[keep, , drop = FALSE])
    ordinary <- abs(values) <= 1
    interval <- cut(values[ordinary], breaks = score_breaks,
      include.lowest = TRUE, labels = FALSE)
    mass <- rowsum(bin_weights[ordinary], group = interval, reorder = FALSE)
    fractions <- numeric(length(score_breaks) - 1)
    fractions[as.integer(rownames(mass))] <- as.vector(mass) / sum(keep)
    # +/-1.1 are special-case flags, not ordinary confidence probabilities.
    fractions <- c(sum(bin_weights[values < -1]) / sum(keep), fractions,
      sum(bin_weights[values > 1]) / sum(keep))
    score_rows[[row_index]] <- tibble(metric, habitat_pair,
      score_center = c(-1.1, head(score_breaks, -1) + 0.025, 1.1),
      special_case = c(TRUE, rep(FALSE, length(score_breaks) - 1), TRUE),
      fraction = fractions)
    stopifnot(abs(sum(fractions) - 1) < tolerance)
    row_index <- row_index + 1
  }
}
score_summary <- bind_rows(score_rows)
write_csv(score_summary, file.path(output_dir, "bin_score_distributions.csv"))
for (plot_metric in c("phylogenetic", "compositional")) {
  score_plot_data <- score_summary %>%
    filter(.data$metric == .env$plot_metric) %>%
    mutate(habitat_pair = factor(habitat_pair, levels = rev(habitat_order)))
  axis_breaks <- c(-1.1, -0.5, 0, 0.5, 1.1)
  axis_labels <- c("Special -", "-0.5", "0", "0.5", "Special +")
  axis_title <- "Phylogenetic null-comparison confidence"
  if (plot_metric == "compositional") {
    score_plot_data <- filter(score_plot_data, !special_case)
    axis_breaks <- c(-1, -0.5, 0, 0.5, 1)
    axis_labels <- as.character(axis_breaks)
    axis_title <- "Compositional null-comparison confidence"
  }
  score_plot <- ggplot(score_plot_data, aes(score_center, habitat_pair, fill = 100 * fraction)) +
    geom_tile(width = 0.05, height = 0.8) +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = "dashed", linewidth = 0.4) +
    scale_fill_viridis_c(name = "Mean bin weight (%)", transform = "sqrt") +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels) +
    guides(fill = guide_colorbar(barwidth = grid::unit(5, "cm"), title.position = "top")) +
    labs(x = axis_title, y = NULL) +
    plot_theme
  ggsave(file.path(figures_dir, paste0("2016_icamp_catchment_", plot_metric, "_scores.pdf")),
    score_plot, width = 8, height = 5.8, bg = "white")
}

# 5. Diagnose ASV turnover, turnover among bins, and absent-bin comparisons -----
bin_ids <- result$detail$taxabin$sp.bin[colnames(community), "bin.id.new"]
bin_counts <- t(rowsum(t(community), group = bin_ids, reorder = FALSE))
bin_counts <- bin_counts[, sub("bin", "", colnames(weights)), drop = FALSE]
stopifnot(all(rowSums(bin_counts) == rowSums(community)))
first <- match(pairs$sample1, rownames(community))
second <- match(pairs$sample2, rownames(community))
absent <- bin_counts[first, ] == 0 | bin_counts[second, ] == 0
asv_bray <- as.matrix(vegdist(community, method = "bray"))
bin_bray <- as.matrix(vegdist(bin_counts, method = "bray"))
turnover <- pairs %>%
  select(sample1, sample2, habitat_pair, comparison) %>%
  mutate(asv_bray = asv_bray[cbind(first, second)],
    bin_bray = bin_bray[cbind(first, second)],
    one_sided_bin_weight = rowSums(weights * absent),
    one_sided_zero_phylo_weight = rowSums(weights * absent * (phylo_scores == 0)),
    one_sided_special_phylo_weight = rowSums(weights * absent * (abs(phylo_scores) > 1)),
    special_phylo_weight = rowSums(weights * (abs(phylo_scores) > 1)))
write_csv(turnover, file.path(output_dir, "pairwise_turnover_diagnostics.csv"))
turnover_summary <- turnover %>%
  group_by(habitat_pair) %>%
  summarize(n_pairs = n(), across(asv_bray:special_phylo_weight, mean), .groups = "drop")
write_csv(turnover_summary, file.path(output_dir, "habitat_turnover_diagnostics.csv"))
turnover_plot_data <- turnover_summary %>%
  pivot_longer(c(asv_bray, bin_bray), names_to = "resolution", values_to = "bray") %>%
  mutate(habitat_pair = factor(habitat_pair, levels = rev(habitat_order)),
    resolution = factor(resolution, levels = c("asv_bray", "bin_bray"),
      labels = c("ASVs", "Phylogenetic bins")))
turnover_plot <- ggplot(turnover_plot_data, aes(bray, habitat_pair, color = resolution)) +
  geom_line(aes(group = habitat_pair), color = "#BBBBBB", orientation = "y") +
  geom_point(size = 2.7) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Mean observed Bray-Curtis dissimilarity", y = NULL, color = NULL) +
  plot_theme
ggsave(file.path(figures_dir, "2016_icamp_catchment_turnover_resolution.pdf"),
  turnover_plot, width = 8, height = 5.8, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info_figures.txt"))
message("Saved four catchment iCAMP figures and diagnostic tables; no null models rerun.")
