#!/usr/bin/env Rscript

# Basic taxonomic alpha, gamma, and beta diversity for the 2016 ASV table.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "helpers.R"))
source(file.path(script_dir, "config.R"))
args <- parse_microbe_args(commandArgs(trailingOnly = TRUE), script_dir)
project_root <- normalizePath(args$project_root, mustWork = TRUE)
output_root <- if (is.na(args$output_root)) {
  under_project_root(microbe_config$output, project_root)
} else {
  under_project_root(args$output_root, project_root)
}
tables_dir <- file.path(output_root, "tables")
derived_dir <- file.path(output_root, "derived")
figures_dir <- under_project_root(microbe_config$figures, project_root)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c("vegan", "iNEXT", "breakaway", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))

input <- load_diversity_input(output_root)
metadata <- input$sample_metadata
raw_counts <- input$asv_counts[metadata$sample_id, , drop = FALSE]

# Raw counts are retained only to describe sequencing depth and detection.
raw_summary <- data.frame(
  sample_id = metadata$sample_id,
  site_code = metadata$site_code,
  habitat = as.character(metadata$habitat),
  reads = rowSums(raw_counts),
  observed_asvs = vegan::specnumber(raw_counts),
  singleton_asvs = rowSums(raw_counts == 1),
  goods_coverage = 1 - rowSums(raw_counts == 1) / rowSums(raw_counts),
  meets_10000_reads = rowSums(raw_counts) >= microbe_config$rarefaction_depth
)
write_analysis_csv(raw_summary,
                   file.path(tables_dir, "raw_sequence_summary.csv"))

# Keep breakaway on the raw table as a richness sensitivity analysis only.
breakaway_results <- lapply(seq_len(nrow(raw_counts)), function(i) {
  frequencies <- as.data.frame(table(raw_counts[i, raw_counts[i, ] > 0]))
  names(frequencies) <- c("index", "frequency")
  frequencies$index <- as.integer(as.character(frequencies$index))
  fit <- try(breakaway::breakaway(frequencies), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(sample_id = metadata$sample_id[i], estimate = NA,
                      standard_error = NA, lower = NA, upper = NA,
                      reasonable = FALSE, warning = as.character(fit)))
  }
  data.frame(
    sample_id = metadata$sample_id[i], estimate = fit$estimate,
    standard_error = fit$error, lower = fit$interval[1], upper = fit$interval[2],
    model = if (is.null(fit$model)) NA_character_ else as.character(fit$model),
    reasonable = fit$reasonable, warning = paste(fit$warnings, collapse = "; ")
  )
})
breakaway_results <- merge(
  do.call(rbind, breakaway_results),
  transform(metadata[c("sample_id", "site_code", "habitat", "nonchim")],
            habitat = as.character(habitat)),
  by = "sample_id", all.x = TRUE, sort = FALSE
)
write_analysis_csv(breakaway_results,
                   file.path(tables_dir, "alpha_breakaway_richness.csv"))

# The main ecological branch starts here. Samples below 10K are excluded, and
# every downstream analysis uses this one saved 10K count table.
keep <- rowSums(raw_counts) >= microbe_config$rarefaction_depth
metadata <- metadata[keep, , drop = FALSE]
original_counts <- raw_counts[keep, , drop = FALSE]
original_reads <- rowSums(original_counts)
set.seed(microbe_config$random_seed)
counts <- vegan::rrarefy(
  original_counts,
  sample = microbe_config$rarefaction_depth
)

saveRDS(
  list(
    asv_counts = counts,
    sample_metadata = metadata,
    rarefaction_depth = microbe_config$rarefaction_depth,
    random_seed = microbe_config$random_seed
  ),
  file.path(derived_dir, "diversity_10k.rds")
)

# Expected Hill numbers at exactly 10K reads, calculated without a random draw.
inext_input <- lapply(seq_len(nrow(original_counts)), function(i) {
  x <- original_counts[i, ]
  x[x > 0]
})
names(inext_input) <- metadata$sample_id
set.seed(microbe_config$random_seed)
inext_alpha <- iNEXT::estimateD(
  inext_input, q = c(0, 1, 2), datatype = "abundance",
  base = "size", level = microbe_config$rarefaction_depth,
  nboot = microbe_config$hill_bootstraps
)
names(inext_alpha)[names(inext_alpha) == "Assemblage"] <- "sample_id"
inext_alpha <- merge(
  inext_alpha,
  transform(metadata[c("sample_id", "site_code", "habitat", "nonchim")],
            habitat = as.character(habitat)),
  by = "sample_id", all.x = TRUE, sort = FALSE
)
write_analysis_csv(inext_alpha,
                   file.path(tables_dir, "alpha_inext_hill_at_10000.csv"))

# Alpha Hill numbers: richness, exponential Shannon, and inverse Simpson.
alpha <- data.frame(
  sample_id = metadata$sample_id,
  site_code = metadata$site_code,
  habitat = as.character(metadata$habitat),
  original_reads = original_reads,
  rarefaction_depth = microbe_config$rarefaction_depth,
  q0 = vegan::specnumber(counts),
  q1 = exp(vegan::diversity(counts, index = "shannon")),
  q2 = vegan::diversity(counts, index = "invsimpson")
)
write_analysis_csv(alpha, file.path(tables_dir, "alpha_hill_numbers.csv"))

alpha_long <- data.frame(
  alpha[rep(seq_len(nrow(alpha)), 3),
        c("sample_id", "site_code", "habitat", "original_reads")],
  order = factor(rep(c("q0", "q1", "q2"), each = nrow(alpha)),
                 levels = c("q0", "q1", "q2")),
  hill_number = c(alpha$q0, alpha$q1, alpha$q2)
)

# Gamma Hill numbers use sample incidence and a common observed coverage, so
# habitats with different numbers of samples are compared at the same coverage.
habitat_rows <- split(seq_len(nrow(counts)), as.character(metadata$habitat))
gamma_input <- lapply(habitat_rows, function(i) t(counts[i, , drop = FALSE] > 0))
gamma_info <- iNEXT::DataInfo(gamma_input, datatype = "incidence_raw")
common_coverage <- min(gamma_info$SC)
set.seed(microbe_config$random_seed)
gamma <- iNEXT::estimateD(
  gamma_input,
  q = c(0, 1, 2),
  datatype = "incidence_raw",
  base = "coverage",
  level = common_coverage,
  nboot = microbe_config$hill_bootstraps
)
names(gamma)[names(gamma) == "Assemblage"] <- "habitat"
write_analysis_csv(gamma, file.path(tables_dir, "gamma_hill_numbers.csv"))

# Hellinger PCA and a simple habitat RDA on the rarefied table.
beta_metadata <- metadata
beta_counts <- counts
beta_counts <- beta_counts[, colSums(beta_counts) > 0, drop = FALSE]
hellinger <- vegan::decostand(beta_counts, method = "hellinger")
pca <- vegan::rda(hellinger)
pca_sites <- vegan::scores(pca, display = "sites", choices = 1:2, scaling = 1)
pca_variance <- vegan::eigenvals(pca) / sum(vegan::eigenvals(pca))
pca_scores <- data.frame(
  sample_id = rownames(pca_sites),
  habitat = as.character(beta_metadata$habitat),
  pc1 = pca_sites[, 1],
  pc2 = pca_sites[, 2]
)
write_analysis_csv(pca_scores, file.path(tables_dir, "hellinger_pca_scores.csv"))

habitat_rda <- vegan::rda(hellinger ~ habitat, data = beta_metadata)
rda_fit <- vegan::RsquareAdj(habitat_rda)
write_analysis_csv(
  data.frame(
    samples = nrow(beta_metadata),
    asvs = ncol(beta_counts),
    rarefaction_depth = microbe_config$rarefaction_depth,
    r_squared = unname(rda_fit$r.squared),
    adjusted_r_squared = unname(rda_fit$adj.r.squared)
  ),
  file.path(tables_dir, "hellinger_rda_habitat_summary.csv")
)

labels <- c(q0 = "q = 0 (richness)", q1 = "q = 1", q2 = "q = 2")
alpha_plot <- ggplot2::ggplot(
  alpha_long,
  ggplot2::aes(habitat, hill_number, color = habitat)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, color = "grey40") +
  ggplot2::geom_jitter(width = 0.14, alpha = 0.7, size = 1.5) +
  ggplot2::facet_wrap(~order, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
  ggplot2::scale_color_manual(values = microbe_config$habitat_colors) +
  ggplot2::scale_x_discrete(labels = microbe_config$habitat_labels) +
  ggplot2::labs(x = NULL, y = "Effective number of ASVs") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
ggplot2::ggsave(file.path(figures_dir, "2016_alpha_hill_by_habitat.pdf"), alpha_plot,
                width = 9, height = 4.5)

inext_alpha$order <- factor(
  paste0("q", inext_alpha$Order.q), levels = c("q0", "q1", "q2")
)
inext_plot <- ggplot2::ggplot(
  inext_alpha,
  ggplot2::aes(habitat, qD, color = habitat)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, color = "grey40") +
  ggplot2::geom_jitter(width = 0.14, alpha = 0.7, size = 1.5) +
  ggplot2::facet_wrap(~order, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
  ggplot2::scale_color_manual(values = microbe_config$habitat_colors) +
  ggplot2::scale_x_discrete(labels = microbe_config$habitat_labels) +
  ggplot2::labs(x = NULL, y = "Effective number of ASVs",
                caption = "iNEXT interpolation to exactly 10,000 reads per sample.") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
ggplot2::ggsave(
  file.path(figures_dir, "2016_inext_alpha_at_10000_by_habitat.pdf"),
  inext_plot, width = 9, height = 4.5
)

pca_plot <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(pc1, pc2, color = habitat)
) +
  ggplot2::geom_point(size = 2.2, alpha = 0.8) +
  ggplot2::scale_color_manual(values = microbe_config$habitat_colors,
                              labels = microbe_config$habitat_labels) +
  ggplot2::labs(
    x = paste0("PC1 (", round(100 * pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_variance[2], 1), "%)"),
    color = "Habitat"
  ) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(figures_dir, "2016_hellinger_pca_by_habitat.pdf"), pca_plot,
                width = 8, height = 6)

capture_analysis_session(file.path(output_root, "audit", "session_info_diversity.txt"))
message("Basic diversity analysis complete.")
