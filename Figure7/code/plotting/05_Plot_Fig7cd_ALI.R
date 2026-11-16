#!/usr/bin/env Rscript
# Figure 7c,d ALI: plot 23,918-cell manuscript Source Data without changing UMAP.
# Reads data/manuscript_plot_inputs/Figure7cd_ALI_Manuscript.csv.
# Point sizes / alpha match the Control plotting configuration.

options(stringsAsFactors = FALSE)

required_packages <- c("readxl", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Run install.packages(c('readxl','ggplot2','patchwork')) and retry.",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

# 1. Read the 23,918-cell manuscript plotting input without changing the values.
input_csv <- Sys.getenv('FIG7_ALI_CSV', unset =
  file.path(getwd(), 'data', 'manuscript_plot_inputs', 'Figure7cd_ALI_Manuscript.csv'))
if (!file.exists(input_csv)) stop('Missing ALI plotting input: ', input_csv)
ali <- utils::read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
names(ali) <- trimws(names(ali))
need <- c('cell_barcode', 'UMAP_1', 'UMAP_2', 'celltype',
          'CCL3_expr', 'CCL4_expr', 'CCL3_scaled_0to1',
          'CCL4_scaled_0to1', 'CCL3_CCL4_coexpression_score')
stopifnot(identical(names(ali), need))
stopifnot(nrow(ali) == 23918L,
          !anyNA(ali$cell_barcode),
          !anyDuplicated(ali$cell_barcode))

# Verify — do not change cell types or values.
expected <- c("Alveolar macrophage" = 11048L,
              "T cell" = 6767L,
              "B cell" = 2325L,
              "Neutrophil" = 3778L)
actual <- table(factor(ali$celltype, levels = names(expected)))
stopifnot(!anyNA(ali$celltype),
          all(ali$celltype %in% names(expected)),
          identical(as.integer(actual), as.integer(expected)))
num_fields <- setdiff(need, c("cell_barcode", "celltype"))
for (nm in num_fields) {
  stopifnot(is.numeric(ali[[nm]]), all(is.finite(ali[[nm]])))
}
for (nm in c("CCL3_scaled_0to1", "CCL4_scaled_0to1",
             "CCL3_CCL4_coexpression_score")) {
  stopifnot(all(ali[[nm]] >= 0 & ali[[nm]] <= 1))
}
message(sprintf("[PASS] Exactly %s source cells, all four existing labels.",
                format(nrow(ali), big.mark = ",")))
message(sprintf("[ORIGINAL UMAP] x=[%.6f, %.6f], y=[%.6f, %.6f]",
                min(ali$UMAP_1), max(ali$UMAP_1),
                min(ali$UMAP_2), max(ali$UMAP_2)))

# 3. Manuscript colors. No UMAP point labels; legends stay OUTSIDE.
cell_colors <- c("Alveolar macrophage" = "#A7C7E7",
                 "T cell" = "#F4B6C2",
                 "B cell" = "#FFD1A9",
                 "Neutrophil" = "#98DF8A")
blue  <- c("#E8F1FB", "#CFE0F6", "#9FC2EA", "#6AA4DF", "#3182BD")
green <- c("#E6F5E6", "#C9ECC9", "#9FE09F", "#72D272", "#31A354")
pink  <- c("#FAE4F3", "#F7CBEA", "#F2ACDF", "#E58BDD", "#D36CD5")

# Factor levels only set the legend order. Original cell classifications stay intact.
ali$celltype <- factor(ali$celltype, levels = names(cell_colors))
base_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        plot.margin = margin(8, 8, 8, 8))
# Approximate outer boundaries read from the supplied manuscript Figure 7 PDF.
# Exact axis tick values match the manuscript. The PDF shows no numeric
# coordinates at the plot edges, so outer boundaries are visual estimates.
# IMPORTANT: no source UMAP values are altered in any way.
frames <- list(
  cell = list(x = c(-15, 7), y = c(-9, 14),
              xb = c(-10, -5, 0, 5), yb = c(-5, 0, 5, 10)),
  expr = list(x = c(-16, 7), y = c(-11, 14),
              xb = c(-15, -10, -5, 0, 5),
              yb = c(-10, -5, 0, 5, 10))
)

# Check for clipping BEFORE plotting. Never silently hide source cells.
for (nm in names(frames)) {
  fr <- frames[[nm]]
  inside <- ali$UMAP_1 >= fr$x[1] & ali$UMAP_1 <= fr$x[2] &
            ali$UMAP_2 >= fr$y[1] & ali$UMAP_2 <= fr$y[2]
  stopifnot(all(inside))
  message(sprintf("[AXIS PASS] %s frame contains %d/%d ALI cells", nm,
                  sum(inside), nrow(ali)))
}

apply_axes <- function(p, panel = c("cell", "expr")) {
  panel <- match.arg(panel)
  fr <- frames[[panel]]
  p +
    scale_x_continuous(breaks = fr$xb, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = fr$yb, expand = expansion(mult = 0)) +
    coord_fixed(xlim = fr$x, ylim = fr$y, ratio = 1, expand = FALSE) +
    labs(x = "umap_1", y = "umap_2") +
    base_theme
}

# Figure 7c: Cell-type UMAP; ALL 23,918 approved cells shown.
p_celltypes <- apply_axes(
  ggplot(ali, aes(UMAP_1, UMAP_2, color = celltype)) +
    geom_point(size = 0.34, alpha = 0.87, stroke = 0) +
    scale_color_manual(values = cell_colors, limits = names(cell_colors),
                       breaks = names(cell_colors), drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
    labs(title = "ALI UMAP-Alveolar macrophage / T cell / B cell / Neutrophil"),
  panel = "cell"
)

# Figure 7d — original expression data and product score, read verbatim.
# Consistent with the original manuscript ALI panel:
# CCL3 and CCL4 use unmodified source CCL3_expr / CCL4_expr with colorbar 1..5;
# values >5 saturate ONLY the display color (not the source values).
# Zero-expression cells appear as faint gray points.
# Co-expression uses the professor-provided 0..1 score, NOT recalculated.
expression_plot <- function(field, title, colors, legend_type = c("raw", "score")) {
  legend_type <- match.arg(legend_type)
  pos <- ali[ali[[field]] > 0, , drop = FALSE]
  pos <- pos[order(pos[[field]]), , drop = FALSE]
  pos$plot_value <- pos[[field]] # only mapping to color; unchanged source value
  p <- ggplot(ali, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "#D5D5D5", size = 0.28, alpha = 0.78, stroke = 0) +
    geom_point(data = pos, aes(color = plot_value),
               size = 0.34, alpha = 0.93, stroke = 0)
  if (legend_type == "raw") {
    p <- p + scale_color_gradientn(colours = colors, limits = c(1, 5),
                breaks = 1:5, oob = scales::squish)
  } else {
    p <- p + scale_color_gradientn(colours = colors, limits = c(0, 1),
                breaks = c(0, 0.25, 0.5, 0.75, 1),
                labels = c("0.00", "0.25", "0.50", "0.75", "1.00"))
  }
  apply_axes(p + labs(title = title), panel = "expr")
}

p_ccl3 <- expression_plot("CCL3_expr", "CCL3 expression", blue, "raw")
p_ccl4 <- expression_plot("CCL4_expr", "CCL4 expression", green, "raw")
p_coexp <- expression_plot("CCL3_CCL4_coexpression_score",
                           "CCL3 & CCL4 co-expression", pink, "score")

# 4. Export — separate folder, verifies write permission before rendering.
out_dir <- file.path(Sys.getenv('FIG7_OUTPUT_DIR', unset = file.path(getwd(), 'figures')), 'ALI')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) {
  stop("Output directory could not be created: ", out_dir, call. = FALSE)
}
test_path <- file.path(out_dir, "_write_test.txt")
if (!isTRUE(file.create(test_path))) {
  stop("Output folder is not writable: ", out_dir, call. = FALSE)
}
unlink(test_path)

save_plot <- function(plot, stem, width, height) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(png_file, plot, width = width, height = height,
                  units = "in", dpi = 600, bg = "white", device = png_device)
  ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                  units = "in", bg = "white", device = "pdf")
  message("[SAVED] ", stem, " (PNG + PDF)")
}

save_plot(p_celltypes, "ALI_Figure7c_CellTypes", 9.2, 7.7)
save_plot(p_ccl3,     "ALI_Figure7d_CCL3", 7.2, 6.2)
save_plot(p_ccl4,     "ALI_Figure7d_CCL4", 7.2, 6.2)
save_plot(p_coexp,    "ALI_Figure7d_CCL3_CCL4_Coexpression", 7.2, 6.2)

# ALI d MANUSCRIPT ORDER (not Control d order):
#  top-left:  co-expression     top-right: intentionally blank
#  bottom-left: CCL3            bottom-right: CCL4
# This is the layout in the professor's Figure 7 PDF.
right_panel <- (p_coexp | patchwork::plot_spacer()) /
               (p_ccl3 | p_ccl4) +
  patchwork::plot_layout(heights = c(1, 1))
combined <- (p_celltypes | right_panel) +
  patchwork::plot_layout(widths = c(1.35, 2))
save_plot(combined, "ALI_Figure7cd_Manuscript_Layout", 17.8, 8.4)

message("[DONE] Manuscript-matched ALI AXES and layout; approved Source Data only")
message("[OUTPUT] ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
