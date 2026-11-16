#!/usr/bin/env Rscript
# Figure 7c,d Control: plot validated 20,327-cell manuscript Source Data.
# Reads data/manuscript_plot_inputs/Figure7cd_Control_Manuscript.csv.
# The manuscript mode uses its display transform; source mode uses raw coordinates.
# Figure 7d Control uses the original log-normalized CCL3_expr/CCL4_expr
# columns with display limits 1–2, matching the 2026-11-16 revised Figure.
# This is a change of plotting variable and color display limits only;
# per-cell source data, UMAP, annotations and coexpression remain unchanged.

options(stringsAsFactors = FALSE)
required_pkgs <- c("readxl", "ggplot2", "patchwork", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop(
  "Install R packages first: install.packages(c(",
  paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))", call. = FALSE
)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# 1. Read the validated manuscript plotting table (same values as the final Source Data).
input_csv <- Sys.getenv('FIG7_CONTROL_CSV', unset =
  file.path(getwd(), 'data', 'manuscript_plot_inputs', 'Figure7cd_Control_Manuscript.csv'))
if (!file.exists(input_csv)) stop('Missing Control plotting input: ', input_csv)
control <- utils::read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
names(control) <- trimws(names(control))
needed <- c('cell_barcode', 'UMAP_1', 'UMAP_2', 'celltype',
  'CCL3_expr', 'CCL4_expr', 'CCL3_scaled_0to1',
  'CCL4_scaled_0to1', 'CCL3_CCL4_coexpression_score')
if (!all(needed %in% names(control))) stop('Missing Control plotting columns')
if (nrow(control) != 20327L) {
  stop("Control cell count must be exactly 20,327, got ", nrow(control), call. = FALSE)
}
if (anyNA(control[, needed]) || anyDuplicated(control$cell_barcode)) {
  stop("Control source has missing required values or duplicated barcodes.", call. = FALSE)
}
num_fields <- c("UMAP_1", "UMAP_2", "CCL3_expr", "CCL4_expr",
                "CCL3_scaled_0to1", "CCL4_scaled_0to1",
                "CCL3_CCL4_coexpression_score")
for (nm in num_fields) {
  if (!is.numeric(control[[nm]]) || any(!is.finite(control[[nm]]))) {
    stop("Source column must contain finite numeric values: ", nm, call. = FALSE)
  }
}
for (nm in c("CCL3_scaled_0to1", "CCL4_scaled_0to1",
             "CCL3_CCL4_coexpression_score")) {
  if (any(control[[nm]] < 0 | control[[nm]] > 1)) {
    stop("0–1 scale out of bounds in source column: ", nm, call. = FALSE)
  }
}
# FINAL SOURCE cell annotations: do not rename, merge or reclassify them.
expected <- c("Alveolar macrophage-like" = 19263L,
              "T cell" = 1038L, "B cell-like" = 14L,
              "Unresolved" = 12L)
if (any(!control$celltype %in% names(expected))) {
  stop("Unexpected Control cell type in validated Source Data.", call. = FALSE)
}
counts <- table(factor(control$celltype, levels = names(expected)))
if (!identical(as.integer(counts), as.integer(expected))) {
  stop("Control annotation counts differ from finalized 20,327-cell source.",
       call. = FALSE)
}
message("[PASS] 20,327 FINAL Control cells: ",
        paste(paste(names(expected), as.integer(expected), sep = "="), collapse = ", "))

# 2. Figure appearance — exact manuscript pastel colors.
colors <- c("Alveolar macrophage-like" = "#A7C7E7",
            "T cell" = "#F4B6C2",
            "B cell-like" = "#FFD1A9",
            "Neutrophil" = "#98DF8A", # Legend only: source contains 0 of these.
            "Unresolved" = "#C9C9C9")
# Manuscript labels are historical display labels, not assignments: B cell is
# B cell-like in the data; Unresolved is drawn grey but omitted from that legend;
# Neutrophil is a legend entry only and the source contains zero such cells.
coord_mode <- Sys.getenv("FIG7_CONTROL_COORD_MODE", unset = "manuscript")
if (!coord_mode %in% c("manuscript", "source"))
  stop("FIG7_CONTROL_COORD_MODE must be manuscript or source.", call. = FALSE)
legend_levels <- c("Alveolar macrophage-like", "T cell", "B cell-like", "Neutrophil")
legend_labels <- c("Alveolar macrophage", "T cell", "B cell", "Neutrophil")
plot_title <- "Control UMAP-Alveolar macrophage / T cell / B cell / Neutrophil"
if (coord_mode == "source") {
  legend_levels <- names(expected)
  legend_labels <- names(expected)
  plot_title <- "Control UMAP (original Source Data coordinates)"
}
# Factor levels and display labels do not alter source cell assignments.
control$celltype <- factor(control$celltype, levels = names(colors))

# Original coordinates are read directly from manuscript Source Data.
# Request: plot ALL 20,327 cells in the EXACT manuscript axis frame.
# This is mathematically impossible with unmodified numeric UMAP coordinates,
# because the original Y range extends below -9. Therefore apply one COMMON
# uniform isotropic display-only transform (SAME scale on X and Y, plus shift).
# No per-cluster warping, no filtering, no re-UMAP, no source data mutation.
# IMPORTANT: the on-plot UMAP values are thus DISPLAY coordinates; this must
# be disclosed in the manuscript legend/methods if these figures are published.
frame <- list(
  xlim = c(-8, 12),
  ylim = c(-9, 7),
  xb = c(-5, 0, 5, 10),
  yb = c(-8, -4, 0, 4)
)

plot_data <- control  # independent data frame; original source stays intact
plot_margin <- 0.03   # 3% margin on each end of the constrained axis
x_range <- range(plot_data$UMAP_1)
y_range <- range(plot_data$UMAP_2)
frame_w <- diff(frame$xlim)
frame_h <- diff(frame$ylim)

# Single isotropic scaling preserves distances/angles/shapes relative to one
# another; translate ONLY the whole original embedding as one unit.
fit_scale <- min(
  (1 - 2*plot_margin) * frame_w / diff(x_range),
  (1 - 2*plot_margin) * frame_h / diff(y_range)
)
source_center_x <- mean(x_range)
source_center_y <- mean(y_range)
frame_center_x <- mean(frame$xlim)
frame_center_y <- mean(frame$ylim)
plot_data$DISPLAY_X <- frame_center_x + fit_scale *
                       (plot_data$UMAP_1 - source_center_x)
plot_data$DISPLAY_Y <- frame_center_y + fit_scale *
                       (plot_data$UMAP_2 - source_center_y)

# Source-coordinate mode plots the exact values in the input, without UMAP
# recalculation, translation, scaling, cell filtering or reclassification.
if (coord_mode == "source") {
  plot_data$DISPLAY_X <- plot_data$UMAP_1
  plot_data$DISPLAY_Y <- plot_data$UMAP_2
  px <- plot_margin * diff(x_range)
  py <- plot_margin * diff(y_range)
  frame <- list(xlim = x_range + c(-px, px), ylim = y_range + c(-py, py),
                xb = pretty(x_range, n = 4), yb = pretty(y_range, n = 4))
  stopifnot(identical(plot_data$DISPLAY_X, control$UMAP_1),
            identical(plot_data$DISPLAY_Y, control$UMAP_2))
}

in_frame <- with(plot_data,
  DISPLAY_X >= frame$xlim[1] & DISPLAY_X <= frame$xlim[2] &
  DISPLAY_Y >= frame$ylim[1] & DISPLAY_Y <= frame$ylim[2]
)
stopifnot(length(in_frame) == 20327L, all(in_frame),
          all(is.finite(plot_data$DISPLAY_X)),
          all(is.finite(plot_data$DISPLAY_Y)),
          identical(plot_data$UMAP_1, control$UMAP_1),
          identical(plot_data$UMAP_2, control$UMAP_2))
message(sprintf("[PASS] %d/%d Control cells visible (%s coordinates)",
                sum(in_frame), nrow(plot_data), coord_mode))
if (coord_mode == "manuscript") {
  # Current input: x_display = 0.7402556035793831*x + 1.7165796921799406;
  # y_display = 0.7402556035793831*y + 0.7933126412823628. The formulas
  # above derive these values from the fixed frame and input coordinate ranges.
  message(sprintf("[DISPLAY TRANSFORM] x = %.16g * UMAP_1 + %.16g; y = %.16g * UMAP_2 + %.16g",
    fit_scale, frame_center_x-fit_scale*source_center_x,
    fit_scale, frame_center_y-fit_scale*source_center_y))
} else {
  message("[SOURCE COORDINATES] DISPLAY_X/Y equal UMAP_1/2 for every cell")
}
message("[SOURCE] ORIGINAL UMAP coordinates retained unchanged in the Source Data")


base_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = .5, size = 12, face = "plain",
                                     margin = ggplot2::margin(b = 7)),
    axis.title = ggplot2::element_text(size = 12, colour = "black"),
    axis.text = ggplot2::element_text(size = 10, colour = "black"),
    axis.line = ggplot2::element_line(colour = "#777777", linewidth = .55),
    axis.ticks = ggplot2::element_line(colour = "#555555", linewidth = .4),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 10),
    legend.key = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )
with_axes <- function(p) {
  p + ggplot2::scale_x_continuous(breaks = frame$xb,
                                 expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(breaks = frame$yb,
                               expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_fixed(ratio = 1, xlim = frame$xlim, ylim = frame$ylim,
                         expand = FALSE) +
    ggplot2::labs(x = "umap_1", y = "umap_2") + base_theme
}

# Figure 7c: grey Unresolved cells are retained. Legend depends on coordinate mode.
# Important: NO geom_text(), geom_label(), ggrepel, or internal type labels.
type_plot_data <- plot_data[order(plot_data$celltype != "Unresolved"), , drop = FALSE]
p_cell <- ggplot2::ggplot(type_plot_data,
                         ggplot2::aes(DISPLAY_X, DISPLAY_Y, colour = celltype)) +
  ggplot2::geom_point(size = .34, alpha = .87, stroke = 0, show.legend = TRUE) +
  ggplot2::scale_colour_manual(values = colors, limits = names(colors),
                               breaks = legend_levels, labels = legend_labels,
                               drop = FALSE) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = 3.6, alpha = 1)
  )) +
  ggplot2::labs(title = plot_title)
p_cell <- with_axes(p_cell) + ggplot2::theme(
  legend.position = "right",    # <- Never inside the UMAP axes.
  legend.justification = "center",
  legend.key.height = grid::unit(.33, "cm")
)

# Figure 7d: the CCL3/CCL4 panels show the ORIGINAL log-normalized
# gene values (CCL3_expr, CCL4_expr), not the 0-1 normalized scores.
# The submitted 2026-11-16 Control Figure displays log-normalized values
# with a COLORBAR from 1 to 2; this is only a display window. Values below
# 1 and above 2 are displayed at the corresponding end color (squished),
# without changing, excluding, truncating or overwriting source values.
# Co-expression remains the existing 0–1 product score and is NOT recomputed.
blue_grad <- c("#E8F1FB", "#CFE0F6", "#9FC2EA", "#6AA4DF", "#3182BD")
green_grad <- c("#E6F5E6", "#C9ECC9", "#9FE09F", "#72D272", "#31A354")
pink_grad <- c("#FAE4F3", "#F7CBEA", "#F2ACDF", "#E58BDD", "#D36CD5")
plot_expression <- function(column, title, gradient,
                            scale_type = c("lognorm", "relative")) {
  scale_type <- match.arg(scale_type)
  # Positive source values are shown; zero-expression cells remain grey.
  highlighted <- plot_data[plot_data[[column]] > 0, , drop = FALSE]
  highlighted <- highlighted[order(highlighted[[column]]), , drop = FALSE]
  highlighted$display_score <- highlighted[[column]]
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(DISPLAY_X, DISPLAY_Y)) +
    ggplot2::geom_point(colour = "#D8D8D8", size = .28,
                        alpha = .78, stroke = 0) +
    ggplot2::geom_point(data = highlighted,
                        ggplot2::aes(colour = display_score),
                        size = .34, alpha = .93, stroke = 0)
  if (scale_type == "lognorm") {
    p <- p + ggplot2::scale_colour_gradientn(
      colours = gradient, limits = c(1, 2),
      breaks = c(1, 1.25, 1.5, 1.75, 2),
      labels = c("1.00", "1.25", "1.50", "1.75", "2.0"),
      oob = scales::squish,
      name = "Log-norm.",
      guide = ggplot2::guide_colourbar(
        title.position = "top", barheight = grid::unit(35, "mm"),
        barwidth = grid::unit(3.6, "mm"), ticks = FALSE))
  } else {
    p <- p + ggplot2::scale_colour_gradientn(
      colours = gradient, limits = c(0, 1),
      breaks = c(0, .25, .50, .75, 1),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      name = "Relative\nscore",
      guide = ggplot2::guide_colourbar(
        title.position = "top", barheight = grid::unit(35, "mm"),
        barwidth = grid::unit(3.6, "mm"), ticks = FALSE))
  }
  with_axes(p + ggplot2::labs(title = title)) + ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 9, hjust = .5),
    legend.text = ggplot2::element_text(size = 9))
}
p_ccl3 <- plot_expression("CCL3_expr",
                         expression(italic(CCL3)~"expression"), blue_grad,
                         "lognorm")
p_ccl4 <- plot_expression("CCL4_expr",
                         expression(italic(CCL4)~"expression"), green_grad,
                         "lognorm")
p_coexpr <- plot_expression("CCL3_CCL4_coexpression_score",
                           expression(italic(CCL3)~"&"~italic(CCL4)~"co-expression"),
                           pink_grad, "relative")

# 3. Export only new plots; never write to the source Excel.
out_dir <- Sys.getenv('FIG7_OUTPUT_DIR', unset = file.path(getwd(), 'figures'))
out_dir <- file.path(out_dir, if (coord_mode == 'source') 'Control_original_coordinates' else 'Control')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Cannot create output directory: ", out_dir, call. = FALSE)

save_pair <- function(p, stem, width, height) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_file, p, width = width, height = height,
                  units = "in", dpi = 600, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf_file, p, width = width, height = height,
                  units = "in", bg = "white", limitsize = FALSE,
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  message("[SAVED] ", png_file, " + PDF")
}
save_pair(p_cell,   "Control_Figure7c_CellTypes",  9.5, 7.8)
save_pair(p_ccl3,   "Control_Figure7d_CCL3",       7.1, 6.2)
save_pair(p_ccl4,   "Control_Figure7d_CCL4",       7.1, 6.2)
save_pair(p_coexpr, "Control_Figure7d_CCL3_CCL4_Coexpression", 7.1, 6.2)

combined <- (p_cell | ((p_ccl3 | p_ccl4) /
                       (p_coexpr | patchwork::plot_spacer()))) +
  patchwork::plot_layout(widths = c(1.45, 1))
save_pair(combined, "Control_Figure7cd_Combined", 18.0, 9.0)

message("[DONE] Control Figure 7c,d figures saved to: ", out_dir)
