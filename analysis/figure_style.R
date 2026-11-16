# Shared figure style for every plotting script in this repository.
# Figures are built at printed size: each panel has a fixed physical size and
# the canvas is measured around it, so the text sizes below are the sizes that
# print. Three standard panels fit across a double-column journal figure
# (about 7.1 in). Source this file after loading ggplot2 and patchwork.

# 1. Text, panel, and mark sizes --------------------------------------------------
# Axis titles, legend titles, and facet labels use the base size; tick labels
# and legend entries are one point smaller. All text is black and plain.
figure_text_size <- 9
figure_small_text_size <- 8
figure_tag_size <- 10
# Panel sizes in inches. The standard panel is about 4:3. A large panel holds a
# single plot that overlays several groups. Map and ordination heights follow
# the data ranges so the axes keep equal scales.
figure_panel_width <- 2.1
figure_panel_height <- 1.6
figure_large_panel_width <- 2.6
figure_large_panel_height <- 2.0
figure_map_width <- 3.0
# Marks at printed scale. Text drawn inside a panel (geom_text) takes its size
# in millimeters, so the label size converts the small text size from points.
figure_point_size <- 1.4
figure_small_point_size <- 0.8
figure_line_width <- 0.5
figure_label_size <- figure_small_text_size / .pt
figure_raster_dpi <- 600

# 2. Theme ------------------------------------------------------------------------
figure_theme <- theme_bw(base_size = figure_text_size) +
  theme(
    text = element_text(color = "black"),
    axis.title = element_text(size = figure_text_size, color = "black"),
    axis.text = element_text(size = figure_small_text_size, color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = figure_text_size, color = "black"),
    legend.title = element_text(size = figure_text_size, face = "plain", color = "black"),
    legend.text = element_text(size = figure_small_text_size, color = "black"),
    legend.key.size = grid::unit(0.35, "cm"),
    legend.position = "bottom",
    plot.tag = element_text(size = figure_tag_size, face = "bold"),
    plot.margin = margin(4, 6, 4, 4)
  )

# Maps drawn in UTM meters are labeled in kilometers, as in Figures 1 and 4.
utm_km_scales <- list(
  scale_x_continuous(labels = function(v) v / 1000),
  scale_y_continuous(labels = function(v) v / 1000),
  labs(x = "Easting (km, UTM zone 10N)", y = "Northing (km)")
)

# 3. Fixed panel sizes and export ---------------------------------------------------
# A single panel.widths value sets the width of all panels together, so faceted
# plots give one width per column and one height per row.
fixed_panels <- function(width = figure_panel_width, height = figure_panel_height, n_col = 1, n_row = 1) {
  theme(panel.widths = grid::unit(rep(width, n_col), "in"),
    panel.heights = grid::unit(rep(height, n_row), "in"))
}

# Panels for maps and ordinations: the height follows the ratio of the data
# ranges so one unit on each axis has the same printed length.
fixed_equal_panels <- function(x, y, width = figure_map_width, n_col = 1, n_row = 1) {
  height <- width * diff(range(y, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
  fixed_panels(width, height, n_col, n_row)
}

# Save a figure on a canvas measured around its fixed panels. Measuring text
# needs an open graphics device, so a null device is opened. A legend below
# the panels is centered on the columns it spans and can be wider than them;
# the table gets an empty column on whichever side the legend would overflow.
# PNG output uses the raster resolution; PDF output is vector.
save_figure <- function(plot, path) {
  pdf(NULL)
  if (inherits(plot, "patchwork")) {
    grob <- patchworkGrob(plot)
  } else {
    grob <- ggplotGrob(plot)
  }
  widths <- grid::convertWidth(grob$widths, "in", valueOnly = TRUE)
  legend_rows <- which(grepl("^guide-box-(bottom|top)", grob$layout$name))
  for (i in legend_rows) {
    legend_width <- grid::convertWidth(grid::grobWidth(grob$grobs[[i]]), "in", valueOnly = TRUE)
    left <- grob$layout$l[i]
    right <- grob$layout$r[i]
    overflow <- (legend_width - sum(widths[left:right])) / 2
    left_pad <- max(0, overflow - sum(widths[seq_len(left - 1)]))
    right_pad <- max(0, overflow - sum(widths[-seq_len(right)]))
    grob <- gtable::gtable_add_cols(grob, grid::unit(left_pad, "in"), pos = 0)
    grob <- gtable::gtable_add_cols(grob, grid::unit(right_pad, "in"), pos = -1)
    widths <- c(left_pad, widths, right_pad)
  }
  width <- sum(widths)
  height <- grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE)
  dev.off()
  ggsave(path, grob, width = width, height = height, dpi = figure_raster_dpi, bg = "white")
}
