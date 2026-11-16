#!/usr/bin/env Rscript
# ===========================================================================
# FIGURE 7f + 7g -- REPRODUCTION USING THE USER-CONFIRMED METHODS
# Source: Figure_7f_7g_plotting_data_GSE58682.xlsx / sheet Sample_Data
#
# FIGURE 7f:
#   CCL3: two-sided Welch two-sample t-test (unequal variances)
#   CCL4: two-sided Welch two-sample t-test (unequal variances)
#   CCL3-CCL4 JOINT PROFILE: two-response MANOVA with Pillai's trace
# FIGURE 7g:
#   Pearson (two-sided), one test for each displayed feature
#
# IMPORTANT: The third 7f boxplot DISPLAYS the mean of z(CCL3), z(CCL4).
# The joint two-response MANOVA is calculated separately and saved to CSV.
# No test is printed on the Figure 7f panel.
# No displayed p-values or significance stars are hardcoded.
# The 56 profiles come from 19 subjects. These profile-level tests do NOT
# adjust for repeated profiles within subject; see ANALYSIS_NOTES.txt.
# ===========================================================================

options(stringsAsFactors = FALSE)
required <- c("readxl", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing packages: install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))")
}
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggplot2))

# 1) INPUT and OUTPUT. Put the .R and Excel together in this Windows folder.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
base_dir <- Sys.getenv('FIG7_BASE_DIR', unset = getwd())
input_candidates <- file.path(base_dir, c(
  "Figure_7f_7g_plotting_data_GSE58682.xlsx",
  "Figure_7f_7g_plotting_data_GSE58682(1).xlsx"
))
existing_inputs <- input_candidates[file.exists(input_candidates)]
default_input <- if (length(existing_inputs)) existing_inputs[1] else input_candidates[1]
input_path <- Sys.getenv("FIG7_INPUT_XLSX", unset = default_input)
output_dir <- file.path(Sys.getenv('FIG7_OUTPUT_DIR', unset = file.path(getwd(),'figures')), 'Fig7fg')
if (!file.exists(input_path)) stop("Missing Excel input: ", input_path)
if (!("Sample_Data" %in% readxl::excel_sheets(input_path))) {
  stop("Excel input needs a sheet named Sample_Data")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2) READ and VALIDATE
D <- as.data.frame(readxl::read_excel(input_path, sheet = "Sample_Data"))
needed <- c("Sample_ID", "Exposure_group", "CCL3_log2", "CCL4_log2",
            "CCL3+CCL4_meanZ", "CCL3_z", "CCL4_z", "CCL3+CCL4_z",
            "AM_activation_score", "AM_activation_score_z")
if (!all(needed %in% names(D))) {
  stop("Missing columns: ", paste(setdiff(needed, names(D)), collapse = ", "))
}
if (nrow(D) != 56L || anyDuplicated(D$Sample_ID)) {
  stop("Expected 56 rows with distinct sample IDs")
}
for (nm in setdiff(needed, c("Sample_ID", "Exposure_group"))) {
  D[[nm]] <- suppressWarnings(as.numeric(D[[nm]]))
}
if (any(!is.finite(as.matrix(D[setdiff(needed, c("Sample_ID", "Exposure_group"))])))) {
  stop("Missing/non-finite measurements found")
}
D$group <- factor(D$Exposure_group, levels = c("Clean Air", "Ozone"))
if (anyNA(D$group) || sum(D$group == "Clean Air") != 18L ||
    sum(D$group == "Ozone") != 38L) {
  stop("Expected 18 Clean Air and 38 Ozone profiles")
}
D$subject_id <- sub("^GSM[0-9]+_([0-9]+)_(000|100|200)$", "\\1", D$Sample_ID)
if (any(D$subject_id == D$Sample_ID)) stop("Unexpected Sample_ID format")
n_subjects <- length(unique(D$subject_id))
if (n_subjects != 19L) warning("Expected 19 subjects; saw ", n_subjects)

# 3) RECOMPUTE NORMALIZATIONS FOR QC, but NEVER overwrite input gene values.
z <- function(x) as.numeric(scale(x))
z3 <- z(D$CCL3_log2)
z4 <- z(D$CCL4_log2)
composite_mean_z <- (z3 + z4) / 2
assert_equal <- function(name, actual, expected, tolerance = 1e-6) {
  if (max(abs(actual - expected)) > tolerance)
    stop("Excel normalization mismatch: ", name)
}
assert_equal("CCL3_z", D$CCL3_z, z3)
assert_equal("CCL4_z", D$CCL4_z, z4)
assert_equal("CCL3+CCL4_meanZ", D$`CCL3+CCL4_meanZ`, composite_mean_z)
assert_equal("CCL3+CCL4_z", D$`CCL3+CCL4_z`, z(composite_mean_z))
assert_equal("AM_activation_score_z", D$AM_activation_score_z,
             z(D$AM_activation_score))

# 4) FIGURE 7f: Welch FOR TWO INDIVIDUAL GENES ONLY
feature_levels <- c("CCL3", "CCL4", "CCL3+CCL4")
welch_for <- function(gene, y) {
  clean <- y[D$group == "Clean Air"]
  ozone <- y[D$group == "Ozone"]
  result <- stats::t.test(ozone, clean, alternative = "two.sided",
                          var.equal = FALSE)
  data.frame(feature = gene, test = "Welch two-sided unequal-variance t-test",
             n_clean = length(clean), n_ozone = length(ozone),
             mean_clean = mean(clean), mean_ozone = mean(ozone),
             t_statistic = unname(result$statistic),
             df = unname(result$parameter), p_raw = result$p.value)
}
welch_stats <- rbind(
  welch_for("CCL3", D$CCL3_log2),
  welch_for("CCL4", D$CCL4_log2)
)

# 5) FIGURE 7f: TWO-RESPONSE MANOVA / PILLAI
# NOTE: CCL3_log2 and CCL4_log2 are the TWO separate response variables.
manova_fit <- stats::manova(cbind(CCL3_log2, CCL4_log2) ~ group, data = D)
Pillai <- summary(manova_fit, test = "Pillai")$stats
joint <- data.frame(
  feature = "CCL3-CCL4 joint profile",
  test = "MANOVA (Pillai trace, two responses)",
  n_clean = sum(D$group == "Clean Air"),
  n_ozone = sum(D$group == "Ozone"),
  Pillai_trace = unname(Pillai[1, "Pillai"]),
  approx_F = unname(Pillai[1, "approx F"]),
  num_df = unname(Pillai[1, "num Df"]),
  den_df = unname(Pillai[1, "den Df"]),
  p_raw = unname(Pillai[1, "Pr(>F)"])
)

p_symbol <- function(p) ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
p_label <- function(p) ifelse(p < 0.001, sprintf("%.3g", p), sprintf("%.4f", p))
welch_stats$significance <- p_symbol(welch_stats$p_raw)
joint$significance <- p_symbol(joint$p_raw)

# 6) FIGURE 7f: three FACETS; the third visually shows mean-z.
value_list <- list(CCL3 = D$CCL3_log2,
                   CCL4 = D$CCL4_log2,
                   "CCL3+CCL4" = D$`CCL3+CCL4_meanZ`)
f_long <- do.call(rbind, lapply(feature_levels, function(f) {
  data.frame(sample = D$Sample_ID, subject_id = D$subject_id,
             group = D$group, feature = f, expression = value_list[[f]])
}))
f_long$feature <- factor(f_long$feature, levels = feature_levels)
row.names(f_long) <- NULL

# Figure 7f display: boxplots, individual points and calculated statistical marks.
# Statistics are calculated above and still exported to CSV below.

set.seed(7)
p7f <- ggplot(f_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(width = 0.46, outlier.shape = NA,
               color = "#777777", alpha = 0.80) +
  geom_point(aes(color = group),
             position = position_jitter(width = 0.10, height = 0, seed = 7),
             size = 1.8, alpha = 0.84, show.legend = FALSE) +
  facet_wrap(~ feature, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Clean Air" = "#F6AFAF", "Ozone" = "#AFC8EB")) +
  scale_color_manual(values = c("Clean Air" = "#E66F6F", "Ozone" = "#3877BB")) +
  scale_x_discrete(labels = c("Clean Air" = "Clean Air", "Ozone" = "Ozone\ninhalation")) +
  scale_y_continuous(expand = expansion(mult = c(0.045, 0.16))) +
  labs(x = NULL, y = "Expression (log2) / combined mean z-score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "#EDEDED", color = "#BBBBBB"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9))

# Figure 7f annotation: first two labels use Welch p, the joint label uses
# Pillai MANOVA of the underlying two-response CCL3/CCL4 profile.
mark <- data.frame(feature=factor(feature_levels,levels=feature_levels),
  label=c(as.character(welch_stats$significance), as.character(joint$significance)),
  y=unname(vapply(feature_levels,function(f) max(f_long$expression[f_long$feature==f],na.rm=TRUE),numeric(1))) +
    unname(vapply(feature_levels,function(f) diff(range(f_long$expression[f_long$feature==f]))*.10,numeric(1))))
mark$label_y <- mark$y + unname(vapply(feature_levels,function(f)
  max(diff(range(f_long$expression[f_long$feature==f]))*.07, .12),numeric(1)))
p7f <- p7f +
  ggplot2::geom_segment(data=mark, ggplot2::aes(x=1,xend=2,y=y,yend=y),
                        inherit.aes=FALSE, linewidth=.4, colour='black') +
  ggplot2::geom_text(data=mark,ggplot2::aes(x=1.5,y=label_y,label=label),
                     inherit.aes=FALSE,size=3.2,colour='black')

# 7) FIGURE 7g: Pearson (raw P values; NO automatic BH adjustment)
# Plot the exact normalized values saved in the manuscript workbook. The
# independent full-precision recomputation remains the established test input.
z_features <- list(CCL3 = D$CCL3_z, CCL4 = D$CCL4_z,
                   "CCL3+CCL4" = D$`CCL3+CCL4_z`)
z_features_for_tests <- list(CCL3 = z3, CCL4 = z4,
                             "CCL3+CCL4" = z(composite_mean_z))
g_long <- do.call(rbind, lapply(feature_levels, function(f) {
  data.frame(sample = D$Sample_ID, subject_id = D$subject_id,
             feature = f, score_z = D$AM_activation_score_z,
             expression_z = z_features[[f]])
}))
g_long$feature <- factor(g_long$feature, levels = feature_levels)
row.names(g_long) <- NULL

g_stats <- do.call(rbind, lapply(feature_levels, function(f) {
  a <- g_long[g_long$feature == f, ]
  ct <- stats::cor.test(a$score_z, z_features_for_tests[[f]],
                        method = "pearson", alternative = "two.sided")
  data.frame(feature = f, n_profiles = nrow(a),
             Pearson_R = unname(ct$estimate), p_raw = ct$p.value)
}))
g_stats$plot_label <- paste0("R = ", sprintf("%.2f", g_stats$Pearson_R),
                             ", P = ", p_label(g_stats$p_raw))
g_anno <- data.frame(feature = factor(g_stats$feature, levels = feature_levels),
                     label = g_stats$plot_label)

p7g <- ggplot(g_long, aes(x = score_z, y = expression_z)) +
  geom_point(aes(color = feature), size = 2.15, alpha = 0.86) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linetype = "dashed", color = "#333333", linewidth = 0.65) +
  geom_text(data = g_anno, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.12, vjust = 1.4, size = 3.0) +
  facet_grid(feature ~ ., scales = "free_y") +
  scale_color_manual(values = c("CCL3" = "#77AFE9", "CCL4" = "#A4D7A8",
                                "CCL3+CCL4" = "#F1A0B8")) +
  labs(
       x = "Macrophage activation score (z-normalized)",
       y = "Gene expression (z-normalized)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "#EDEDED", color = "#BBBBBB"),
        strip.text.y.right = element_text(angle = 270),
        panel.grid.minor = element_blank(),
        panel.spacing.y = grid::unit(0.16, "cm"))

# 8) EXPORT EXACTLY WHAT WAS CALCULATED
write.csv(f_long, file.path(output_dir, "Figure7f_plotted_points.csv"), row.names = FALSE)
write.csv(welch_stats, file.path(output_dir, "Figure7f_Welch_CCL3_CCL4.csv"), row.names = FALSE)
write.csv(joint, file.path(output_dir, "Figure7f_MANOVA_Pillai_joint_CCL3_CCL4.csv"), row.names = FALSE)
write.csv(g_long, file.path(output_dir, "Figure7g_plotted_points.csv"), row.names = FALSE)
write.csv(g_stats, file.path(output_dir, "Figure7g_Pearson_statistics.csv"), row.names = FALSE)

ggsave(file.path(output_dir, "Figure7f_Welch_CCL3_CCL4_MANOVA_joint.png"),
       plot = p7f, width = 9.0, height = 5.2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Figure7f_Welch_CCL3_CCL4_MANOVA_joint.pdf"),
       plot = p7f, width = 9.0, height = 5.2, bg = "white")
ggsave(file.path(output_dir, "Figure7g_Pearson.png"),
       plot = p7g, width = 9.0, height = 5.4, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Figure7g_Pearson.pdf"),
       plot = p7g, width = 9.0, height = 5.4, bg = "white")

writeLines(c(
  "Source workbook: Figure_7f_7g_plotting_data_GSE58682.xlsx / Sample_Data",
  "Figure 7f methods: Welch two-sided t-tests for CCL3, CCL4; MANOVA Pillai trace for the two-response joint CCL3-CCL4 profile.",
  "Figure 7f third panel displays composite mean(z(CCL3),z(CCL4)); joint MANOVA uses the two underlying variables.",
  paste0("CCL3 Welch P = ", format(welch_stats$p_raw[1], digits = 10)),
  paste0("CCL4 Welch P = ", format(welch_stats$p_raw[2], digits = 10)),
  paste0("Joint MANOVA Pillai = ", format(joint$Pillai_trace, digits = 10),
         "; F = ", format(joint$approx_F, digits = 10),
         "; df = ", joint$num_df, ",", joint$den_df,
         "; P = ", format(joint$p_raw, digits = 10)),
  "Figure 7g: two-sided Pearson correlation (raw P) per feature; original manuscript numbers are not forced.",
  paste0("Profiles = 56; subjects = ", n_subjects,
         ". All tests above ignore within-subject clustering. Treat significance as exploratory."),
  "Figure 7f displays Welch stars (CCL3,CCL4) and joint-Pillai MANOVA star (CCL3+CCL4); Figure 7g Pearson labels are included."
), file.path(output_dir, "ANALYSIS_NOTES.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))

# 9) SHOW the figures in RStudio's Plots pane as well.
if (interactive()) { print(p7f); print(p7g) }
cat("\n=== FIGURE 7f: CONFIRMED TEST TYPES, RECOMPUTED VALUES ===\n")
print(welch_stats, row.names = FALSE, digits = 10)
print(joint, row.names = FALSE, digits = 10)
cat("\n=== FIGURE 7g: PEARSON ===\n")
print(g_stats, row.names = FALSE, digits = 10)
cat("\nSaved figures and CSVs to:\n", normalizePath(output_dir, winslash = "/"), "\n")
cat("Figure 7f plot object: p7f; Figure 7g plot object: p7g\n")
