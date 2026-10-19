# Descriptive catalogue of Danczak mass-difference matches and assignment checks.
# Uses completed outputs from supporting 08; checks the curated motifs directly.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Settings and readable descriptions -----------------------------------------
results_dir <- here("results", "fticr_2016", "transformations")
figures_dir <- here("figures", "supplement", "transformations")
mass_tolerance <- 0.000010
motif_descriptions <- tribble(
  ~motif, ~compatible_change,
  "H2", "Hydrogen gain/loss",
  "O", "Oxygen gain/loss",
  "H2O", "Water gain/loss",
  "CO2", "Carboxyl-group equivalent",
  "CH2", "CH2 homologous-series difference",
  "NH2", "Amino-group equivalent",
  "HPO3", "Phosphate-group equivalent",
  "SO3", "SO3 gain/loss",
  "C6H10O5", "Hexose-residue equivalent",
  "C2H2O", "Acetyl-residue equivalent",
  "Glycine residue", "Glycine-residue equivalent",
  "Valine residue", "Valine-residue equivalent"
)
# Net elemental differences from lower to higher measured m/z do not imply flux.
motif_design <- tribble(
  ~motif, ~reference_label, ~C, ~H, ~O, ~N, ~S, ~P, ~Na,
  "H2", "hydrogenation/dehydrogenation", 0, 2, 0, 0, 0, 0, 0,
  "O", "hydroxylation (-H)", 0, 0, 1, 0, 0, 0, 0,
  "H2O", "condensation/dehydration", 0, 2, 1, 0, 0, 0, 0,
  "CO2", "Carboxylation", 1, 0, 2, 0, 0, 0, 0,
  "CH2", "methanol (-H2O)", 1, 2, 0, 0, 0, 0, 0,
  "NH2", "primary amine", 0, 2, 0, 1, 0, 0, 0,
  "HPO3", "phosphate", 0, 1, 3, 0, 0, 1, 0,
  "SO3", "sulfate (-H2O)", 0, 0, 3, 0, 1, 0, 0,
  "C6H10O5", "C6H10O5", 6, 10, 5, 0, 0, 0, 0,
  "C2H2O", "acetylation (-H2O)", 2, 2, 1, 0, 0, 0, 0,
  "Glycine residue", "Glycine", 2, 3, 1, 1, 0, 0, 0,
  "Valine residue", "Valine", 5, 9, 1, 1, 0, 0, 0
)
support_colors <- c(
  "Compatible; single formula candidates" = "#0072B2",
  "Compatible; multiple formula candidates" = "#E69F00",
  "Mass match; incompatible assigned formulas" = "#D9D9D9"
)
plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_blank())

# 2. Align completed pair assignments with the sediment inventory and formulas ---
input_paths <- c(
  file.path(results_dir, "transformation_profiles.csv"),
  file.path(results_dir, "peak_pair_assignments.csv.gz"),
  here("data", "derived", "fticr_2016", "sediment_44_presence_absence.csv"),
  here("results", "fticr_2016", "tables", "primary_peak_properties.csv"),
  here("reference", "danczak_2020", "Transformation_Database_07-2020.csv")
)
input_md5 <- unname(tools::md5sum(input_paths))
stopifnot(!anyNA(input_md5))
profiles <- read_csv(input_paths[1], show_col_types = FALSE)
assignments <- read_csv(input_paths[2], show_col_types = FALSE)
pa_table <- read_csv(input_paths[3], col_types = cols(site_code = col_character()))
peaks <- read_csv(input_paths[4], show_col_types = FALSE)
reference <- read_csv(input_paths[5], show_col_types = FALSE)
counts <- as.matrix(profiles[, -(1:3)])
pa <- as.matrix(pa_table[, -1])
rownames(pa) <- pa_table$site_code
stopifnot(ncol(counts) == 44, all(colSums(counts) > 0),
  identical(colnames(counts), rownames(pa)),
  identical(peaks$peak_id, colnames(pa)))

# 3. Check the curated formula differences without environmental models --------
pairs <- assignments %>% filter(published_assignment)
element_columns <- c(C = "carbon_count", H = "hydrogen_count", O = "oxygen_count",
  N = "nitrogen_count", S = "sulfur_count", P = "phosphorus_count", Na = "sodium_count")
lower <- match(pairs$lower_peak_id, peaks$peak_id)
upper <- match(pairs$upper_peak_id, peaks$peak_id)
stopifnot(!anyNA(lower), !anyNA(upper))
delta <- as.matrix(peaks[upper, element_columns]) -
  as.matrix(peaks[lower, element_columns])
colnames(delta) <- names(element_columns)
co_detected <- pa[, lower, drop = FALSE] > 0 & pa[, upper, drop = FALSE] > 0
stopifnot(all(colSums(co_detected) == pairs$observed_sites))
motif_design$mass <- reference$Mass[match(motif_design$reference_label, reference$Name)]
stopifnot(!anyNA(motif_design$mass))
motif_match <- matrix(FALSE, nrow(pairs), nrow(motif_design),
  dimnames = list(pairs$pair_id, motif_design$motif))
motif_audit <- vector("list", nrow(motif_design))
for (j in seq_len(nrow(motif_design))) {
  mass_match <- pairs$observed_difference + mass_tolerance >= motif_design$mass[j] &
    pairs$observed_difference - mass_tolerance <= motif_design$mass[j]
  formula_match <- rowSums(sweep(delta, 2,
    as.numeric(motif_design[j, names(element_columns)]), "!=") != 0) == 0
  motif_match[, j] <- mass_match & formula_match
  motif_audit[[j]] <- motif_design[j, ] %>% mutate(
    mass_matched_pairs = sum(mass_match), compatible_pairs = sum(mass_match & formula_match),
    compatible_single_candidate_pairs = sum(mass_match & formula_match &
      peaks$candidate_count[lower] == 1 & peaks$candidate_count[upper] == 1))
}
motifs <- bind_rows(motif_audit)
write_csv(motifs, file.path(results_dir, "motif_formula_audit.csv"))
motif_profiles <- co_detected %*% (1L * motif_match)
site_motifs <- bind_cols(tibble(site_code = rownames(pa)), as_tibble(motif_profiles)) %>%
  pivot_longer(-site_code, names_to = "motif", values_to = "pairs") %>%
  mutate(fraction_all_matches = pairs / rep(colSums(counts), each = nrow(motif_design)))
write_csv(site_motifs, file.path(results_dir, "site_formula_compatible_motifs.csv"))

# 4. Rank every database label, retaining overwritten alternative annotations ----
# A repeated peak pair is counted once across sites in the distinct-pair columns.
# Sample-pair counts and mean fractions describe occurrence, not reaction evidence.
candidate_summary <- assignments %>% group_by(database_row) %>% summarize(
  candidate_distinct_pairs = n(),
  pairs_with_alternative_labels = sum(candidate_assignments > 1),
  median_absolute_mass_error = median(abs(mass_error)), .groups = "drop")
published_summary <- assignments %>% filter(published_assignment) %>%
  count(database_row, name = "published_distinct_pairs")
catalogue <- profiles %>% select(database_row, transformation, mass) %>%
  mutate(
    sites_with_published_matches = rowSums(counts > 0),
    sample_pair_matches = rowSums(counts),
    mean_sample_pair_fraction = rowMeans(sweep(counts, 2, colSums(counts), "/"))
  ) %>%
  left_join(published_summary, by = "database_row") %>%
  left_join(candidate_summary, by = "database_row") %>%
  mutate(across(c(published_distinct_pairs, candidate_distinct_pairs,
    pairs_with_alternative_labels), ~coalesce(.x, 0L))) %>%
  arrange(desc(mean_sample_pair_fraction), database_row) %>%
  mutate(rank = row_number(), .before = 1)
write_csv(catalogue, file.path(results_dir, "transformation_label_catalogue.csv"))

# 5. Summarize the interpretable subset by frequency and formula-assignment checks
coverage <- site_motifs %>% group_by(motif) %>% summarize(
  sites_with_compatible_pairs = sum(pairs > 0),
  compatible_sample_pair_matches = sum(pairs),
  mean_sample_pair_fraction = mean(fraction_all_matches), .groups = "drop")
motif_summary <- motifs %>%
  left_join(motif_descriptions, by = "motif") %>%
  left_join(coverage, by = "motif") %>%
  transmute(motif, compatible_change, reference_label, mass,
    mass_matched_distinct_pairs = mass_matched_pairs,
    formula_compatible_distinct_pairs = compatible_pairs,
    compatible_single_candidate_pairs,
    compatible_multiple_candidate_pairs = compatible_pairs - compatible_single_candidate_pairs,
    mass_only_incompatible_pairs = mass_matched_pairs - compatible_pairs,
    sites_with_compatible_pairs, compatible_sample_pair_matches, mean_sample_pair_fraction) %>%
  arrange(desc(mass_matched_distinct_pairs))
stopifnot(!anyNA(motif_summary),
  all(motif_summary$compatible_multiple_candidate_pairs >= 0),
  all(motif_summary$mass_only_incompatible_pairs >= 0),
  sum(catalogue$published_distinct_pairs) == n_distinct(assignments$pair_id),
  sum(catalogue$sample_pair_matches) == sum(counts))
write_csv(motif_summary, file.path(results_dir, "interpretable_motif_summary.csv"))

# 6. Plot occurrence and assignment checks without assigning reaction confidence -
support_labels <- setNames(
  paste0(motif_summary$compatible_change, " (",
    recode(motif_summary$motif, "Glycine residue" = "C2H3NO",
      "Valine residue" = "C5H9NO"), ")"),
  motif_summary$motif)
support_distribution <- motif_summary %>%
  select(motif, compatible_single_candidate_pairs, compatible_multiple_candidate_pairs,
    mass_only_incompatible_pairs) %>%
  pivot_longer(-motif, names_to = "assignment_check", values_to = "distinct_pairs") %>%
  mutate(
    assignment_check = factor(assignment_check,
      levels = c("compatible_single_candidate_pairs", "compatible_multiple_candidate_pairs",
        "mass_only_incompatible_pairs"), labels = names(support_colors)),
    motif = factor(motif, levels = rev(motif_summary$motif))
  )
support_plot <- ggplot(support_distribution,
  aes(distinct_pairs, motif, fill = assignment_check)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  scale_y_discrete(labels = support_labels) +
  scale_fill_manual(values = support_colors, name = NULL) +
  labs(x = "Distinct mass-matched peak pairs across 44 profiles", y = NULL) +
  plot_theme + theme(legend.position = "bottom",
    legend.text = element_text(size = 10), legend.key.height = grid::unit(0.35, "cm")) +
  guides(fill = guide_legend(ncol = 1))
ggsave(file.path(figures_dir, "2016_transformation_assignment_summary.pdf"),
  support_plot, width = 11.5, height = 6.5, bg = "white")

rank_plot <- ggplot(filter(catalogue, sample_pair_matches > 0),
  aes(rank, mean_sample_pair_fraction)) +
  geom_point(size = 1.4, color = "#0072B2", alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_percent()) +
  labs(x = "Published transformation-label rank", y = "Mean fraction of matched pairs") +
  plot_theme
ggsave(file.path(figures_dir, "2016_transformation_label_rank_distribution.pdf"),
  rank_plot, width = 6, height = 4.8, bg = "white")

# 7. Preserve source identity ---------------------------------------------------
integrity <- tibble(path = substring(input_paths, nchar(here()) + 2L),
  md5_before = input_md5, md5_after = unname(tools::md5sum(input_paths))) %>%
  mutate(unchanged = md5_before == md5_after)
write_csv(integrity, file.path(results_dir, "summary_input_integrity.csv"))
stopifnot(all(integrity$unchanged))
message("Descriptive summary: ", sum(catalogue$sample_pair_matches > 0),
  " published labels; ", nrow(motif_summary), " interpretable motifs with assignment checks.")
