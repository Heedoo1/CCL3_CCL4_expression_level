#!/usr/bin/env Rscript

# ============================================================
# AM activation score vs CCL3 / CCL4
# Pearson correlation, z-normalized
#
# INPUT  : df_corr.csv
#          (Sample, AM_score, CCL3, CCL4)
# OUTPUT : 
#   1) Figure_AM_vs_CCL3_CCL4_Pearson.png
#   2) Stats_AM_vs_CCL3_CCL4_Pearson.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# ===============================
# 1. Оролтын файл
# ===============================
in_file <- "df_corr.csv"   # GitHub-д хамт оруулна
out_fig <- "Figure_AM_vs_CCL3_CCL4_Pearson.png"
out_stat <- "Stats_AM_vs_CCL3_CCL4_Pearson.csv"

df_corr <- read_csv(in_file, show_col_types = FALSE)

# ===============================
# 2. Z-normalization
# ===============================
z <- function(x) as.numeric(scale(x))

df_corr <- df_corr %>%
  mutate(
    AM_z   = z(AM_score),
    CCL3_z = z(CCL3),
    CCL4_z = z(CCL4),
    CCL3_CCL4_z = z((CCL3_z + CCL4_z) / 2)
  )

# ===============================
# 3. Long формат (3 панель)
# ===============================
df_long <- df_corr %>%
  select(Sample, AM_z, CCL3_z, CCL4_z, CCL3_CCL4_z) %>%
  pivot_longer(
    cols = c(CCL3_z, CCL4_z, CCL3_CCL4_z),
    names_to = "Marker",
    values_to = "Expr_z"
  ) %>%
  mutate(
    Marker = recode(
      Marker,
      "CCL3_z" = "CCL3",
      "CCL4_z" = "CCL4",
      "CCL3_CCL4_z" = "CCL3+CCL4"
    ),
    Marker = factor(Marker, levels = c("CCL3", "CCL4", "CCL3+CCL4"))
  ) %>%
  drop_na()

# ===============================
# 4. Pearson correlation
# ===============================
stat_tbl <- df_long %>%
  group_by(Marker) %>%
  summarise(
    n = n(),
    R = cor(AM_z, Expr_z, method = "pearson"),
    P = cor.test(AM_z, Expr_z, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("R = %.2f,  P = %.3g", R, P)
  )

write.csv(stat_tbl, out_stat, row.names = FALSE)

# Label байрлал
label_pos <- df_long %>%
  group_by(Marker) %>%
  summarise(
    x = min(AM_z) + 0.05 * diff(range(AM_z)),
    y = max(Expr_z) - 0.05 * diff(range(Expr_z)),
    .groups = "drop"
  ) %>%
  left_join(stat_tbl, by = "Marker")

# ===============================
# 5. Зураг
# ===============================
p <- ggplot(df_long, aes(x = AM_z, y = Expr_z)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.6) +
  facet_grid(Marker ~ ., scales = "free_y", switch = "y") +
  geom_text(
    data = label_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 4
  ) +
  labs(
    title = "AM activation score vs CCL3 / CCL4 (Pearson correlation, normalized)",
    x = "Macrophage activation score (z-normalized)",
    y = "Gene expression (z-normalized)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid.minor = element_blank()
  )

ggsave(out_fig, p, width = 8.5, height = 5.5, dpi = 300)
