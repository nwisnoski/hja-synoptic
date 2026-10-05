# Repeated bipartite module fitting for observed, null, and resampled networks.
# BRIM alternates assignments between the two node sets (Barber 2007).
# It maximizes a heuristic lower bound on modularity, not a proven global optimum.
bipartite_modules <- function(web, start_groups, max_sweeps, tolerance) {
  row_active <- which(rowSums(web) > 0)
  column_active <- which(colSums(web) > 0)
  row_membership <- rep(NA_integer_, nrow(web))
  column_membership <- rep(NA_integer_, ncol(web))
  if (length(row_active) == 0L) {
    return(list(modularity = NA_real_, row = row_membership,
      column = column_membership, modules = 0L, converged = TRUE,
      starts = data.frame(start = integer(), modularity = numeric(),
        sweeps = integer(), converged = logical())))
  }
  active <- web[row_active, column_active, drop = FALSE]
  row_degree <- rowSums(active)
  column_degree <- colSums(active)
  edges <- sum(active)
  components <- igraph::components(igraph::graph_from_biadjacency_matrix(active))$membership
  best_q <- -Inf
  diagnostics <- vector("list", length(start_groups))
  for (start in seq_along(start_groups)) {
    if (start == 1L) {
      row_labels <- components[seq_len(nrow(active))]
    } else if (start == 2L) {
      row_labels <- seq_len(nrow(active))
    } else {
      # Keep disconnected components separate: merging them only lowers Q.
      row_labels <- integer(nrow(active))
      next_label <- 0L
      for (component in unique(components[seq_len(nrow(active))])) {
        positions <- which(components[seq_len(nrow(active))] == component)
        groups <- min(start_groups[start], length(positions))
        row_labels[positions] <- next_label + sample.int(groups, length(positions), replace = TRUE)
        next_label <- next_label + groups
      }
    }
    previous_q <- -Inf
    converged <- FALSE
    for (sweep in seq_len(max_sweeps)) {
      row_groups <- sort(unique(row_labels))
      column_scores <- t(rowsum(active, row_labels, reorder = TRUE)) -
        outer(column_degree, as.vector(rowsum(row_degree, row_labels, reorder = TRUE))) / edges
      column_labels <- row_groups[max.col(column_scores, ties.method = "first")]
      column_groups <- sort(unique(column_labels))
      row_scores <- t(rowsum(t(active), column_labels, reorder = TRUE)) -
        outer(row_degree, as.vector(rowsum(column_degree, column_labels, reorder = TRUE))) / edges
      row_labels <- column_groups[max.col(row_scores, ties.method = "first")]
      same_module <- outer(row_labels, column_labels, "==")
      q <- sum((active - outer(row_degree, column_degree) / edges)[same_module]) / edges
      if (is.finite(previous_q)) {
        stopifnot(q >= previous_q - tolerance)
        if (abs(q - previous_q) <= tolerance) {
          converged <- TRUE
          break
        }
      }
      previous_q <- q
    }
    diagnostics[[start]] <- data.frame(start, modularity = q, sweeps = sweep, converged)
    if (q > best_q + tolerance) {
      best_q <- q
      best_row <- row_labels
      best_column <- column_labels
    }
  }
  labels <- sort(unique(c(best_row, best_column)))
  row_membership[row_active] <- match(best_row, labels)
  column_membership[column_active] <- match(best_column, labels)
  return(list(modularity = best_q, row = row_membership, column = column_membership,
    modules = length(labels), converged = all(vapply(diagnostics, function(x) x$converged, TRUE)),
    starts = dplyr::bind_rows(diagnostics)))
}

# Module assignment is repeated thousands of times; centralize the common metrics.
bipartite_statistics <- function(associations, edge_threshold, start_groups,
    max_sweeps, tolerance) {
  web <- 1L * (associations >= edge_threshold)
  fit <- bipartite_modules(web, start_groups, max_sweeps, tolerance)
  statistics <- data.frame(positive_edges = sum(web),
    negative_edges = sum(associations <= -edge_threshold),
    active_asv_patterns = sum(rowSums(web) > 0),
    active_molecular_patterns = sum(colSums(web) > 0),
    modularity = fit$modularity, modules = fit$modules, converged = fit$converged)
  return(list(statistics = statistics, fit = fit))
}
