# Independent numerical and reproducibility checks for the bipartite analysis.
library(here)
library(readr)
library(dplyr)
library(igraph)

# 1. Settings and analytical examples -------------------------------------------
results_dir <- here("results", "fticr_2016", "bipartite_networks")
pilot_dir <- here("results", "fticr_2016", "bipartite_network_pilot")
tolerance <- 1e-10
source(here("analysis", "fticr", "bipartite_helpers.R"))
blocks <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1), 4, 4)
set.seed(1)
toy <- bipartite_modules(blocks, c(0L, 0L, 2L, 3L), 100L, tolerance)
stopifnot(abs(toy$modularity - 0.5) < tolerance, toy$converged)
complete <- bipartite_modules(matrix(1L, 3, 4), c(0L, 0L, 2L, 3L), 100L, tolerance)
stopifnot(abs(complete$modularity) < tolerance)
empty <- bipartite_modules(matrix(0L, 3, 4), c(0L, 0L, 2L, 3L), 100L, tolerance)
stopifnot(is.na(empty$modularity), all(is.na(empty$row)), all(is.na(empty$column)))

# 2. Compare Pearson correlation with the analytical binary phi formula ----------
x_table <- read_csv(file.path(results_dir, "asv_pattern_incidence.csv"), show_col_types = FALSE,
  col_types = cols(site_code = col_character()))
y_table <- read_csv(file.path(results_dir, "molecular_pattern_incidence.csv"), show_col_types = FALSE,
  col_types = cols(site_code = col_character()))
stopifnot(identical(x_table$site_code, y_table$site_code))
x <- as.matrix(x_table[, -1])
y <- as.matrix(y_table[, -1])
n_sites <- nrow(x)
px <- colSums(x)
py <- colSums(y)
phi <- (n_sites * crossprod(x, y) - outer(px, py)) /
  sqrt(outer(px * (n_sites - px), py * (n_sites - py)))
phi_error <- max(abs(phi - cor(x, y)))
stopifnot(phi_error < tolerance)
settings <- read_csv(file.path(results_dir, "settings.csv"), show_col_types = FALSE)
observed <- read_csv(file.path(results_dir, "observed_network.csv"), show_col_types = FALSE)
stopifnot(observed$positive_edges == sum(phi >= settings$edge_threshold),
  observed$negative_edges == sum(phi <= -settings$edge_threshold))

# 3. Independently calculate the module objective and check both coordinate steps -
memberships <- read_csv(file.path(results_dir, "pattern_membership.csv"), show_col_types = FALSE)
row_labels <- memberships$module[memberships$node_type == "ASV"]
column_labels <- memberships$module[memberships$node_type == "Molecular"]
web <- 1L * (phi >= settings$edge_threshold)
edges <- sum(web)
manual_q <- 0
for (module_id in sort(unique(na.omit(c(row_labels, column_labels))))) {
  rows <- which(row_labels == module_id)
  columns <- which(column_labels == module_id)
  within_edges <- sum(web[rows, columns, drop = FALSE])
  row_degree <- sum(rowSums(web)[rows])
  column_degree <- sum(colSums(web)[columns])
  manual_q <- manual_q + within_edges / edges - row_degree * column_degree / edges^2
}
stopifnot(abs(manual_q - observed$modularity) < tolerance)
active_rows <- which(rowSums(web) > 0)
active_columns <- which(colSums(web) > 0)
active <- web[active_rows, active_columns, drop = FALSE]
b <- active - outer(rowSums(active), colSums(active)) / edges
row_scores <- t(rowsum(t(b), column_labels[active_columns]))
column_scores <- t(rowsum(b, row_labels[active_rows]))
row_gain <- sum(apply(row_scores, 1, max)) / edges - manual_q
column_gain <- sum(apply(column_scores, 1, max)) / edges - manual_q
stopifnot(abs(row_gain) < tolerance, abs(column_gain) < tolerance)

# 4. Verify all nulls and the full-size pilot agree at the same per-run seeds ------
nulls <- read_csv(file.path(results_dir, "null_runs.csv"), show_col_types = FALSE)
stopifnot(nrow(nulls) == 4L * settings$n_null, all(nulls$successful),
  all(nulls$margins_preserved[nulls$scheme %in% c("graph_fixed_degrees", "molecular_fixed_margins")]))
pilot_available <- file.exists(file.path(pilot_dir, "null_runs.csv"))
reproduced_pilot_runs <- 0L
if (pilot_available) {
  pilot <- read_csv(file.path(pilot_dir, "null_runs.csv"), show_col_types = FALSE)
  comparison <- inner_join(pilot, nulls, by = c("scheme", "run", "seed", "module_seed"),
    suffix = c("_pilot", "_full"))
  columns <- c("positive_edges", "negative_edges", "active_asv_patterns", "active_molecular_patterns",
    "modularity", "modules", "converged", "changed_fraction")
  for (column in columns) {
    stopifnot(isTRUE(all.equal(comparison[[paste0(column, "_pilot")]],
      comparison[[paste0(column, "_full")]], tolerance = tolerance)))
  }
  stopifnot(nrow(comparison) == nrow(pilot))
  reproduced_pilot_runs <- nrow(comparison)
}
write_csv(tibble(toy_block_modularity = toy$modularity, maximum_phi_error = phi_error,
  objective_error = abs(manual_q - observed$modularity), row_update_gain = row_gain,
  column_update_gain = column_gain, checked_null_runs = nrow(nulls),
  pilot_available, reproduced_pilot_runs, passed = TRUE),
  file.path(results_dir, "numerical_validation.csv"))
message("All bipartite numerical and reproducibility checks passed.")
