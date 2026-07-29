install.packages("showtext")

library(tidyverse)
library(ggsci)
library(scales)
library(showtext)

# ============================================================
# 0. Font: Times New Roman
# ============================================================

# For macOS
font_add(
  family = "Times New Roman",
  regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
)

showtext_auto()

theme_bgc <- theme_bw(base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman", color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Number of BGC types
n_bgc_type <- length(unique(bgc_stats$BGC_type))

# NPG color palette
npg_cols <- colorRampPalette(
  ggsci::pal_npg("nrc")(10)
)(n_bgc_type)

names(npg_cols) <- sort(unique(bgc_stats$BGC_type))

cell_cols <- c(
  "#3C5488FF", "#00A087FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
  "#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF",
  "#F39B7FFF", "#91D1C2FF", "#8491B4FF", "#7E6148FF",
  "#B09C85FF", "#A73030FF", "#0072B5FF", "#E18727FF"
)

cell_cols <- colorRampPalette(cell_cols)(n_bgc_type)
names(cell_cols) <- sort(unique(bgc_stats$BGC_type))

plot_cols <- npg_cols
# plot_cols <- cell_cols

bgc_stats_pct <- bgc_stats %>%
  group_by(BAPS_level1) %>%
  mutate(
    group_total = sum(total_count, na.rm = TRUE),
    percentage = total_count / group_total
  ) %>%
  ungroup()

p1_pct <- ggplot(
  bgc_stats_pct,
  aes(
    x = factor(BAPS_level1),
    y = percentage,
    fill = BGC_type
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.15,
    width = 0.75
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = plot_cols) +
  theme_bgc +
  labs(
    x = "BAPS_level1",
    y = "BGC type percentage",
    fill = "BGC type",
    title = "Relative composition of BGC types by BAPS_level1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p1_pct

ggsave(
  "BGC_percentage_by_BAPS_level1_stacked_bar.pdf",
  p1_pct,
  width = 10,
  height = 6,
  device = cairo_pdf
)

ggsave(
  "BGC_percentage_by_BAPS_level1_stacked_bar.png",
  p1_pct,
  width = 10,
  height = 6,
  dpi = 600
)

p1_pct_label <- ggplot(
  bgc_stats_pct,
  aes(
    x = factor(BAPS_level1),
    y = percentage,
    fill = BGC_type
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.15,
    width = 0.75
  ) +
  geom_text(
    aes(label = ifelse(percentage >= 0.05, percent(percentage, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    family = "Times New Roman"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = plot_cols) +
  theme_bgc +
  labs(
    x = "BAPS_level1",
    y = "BGC type percentage",
    fill = "BGC type",
    title = "Relative composition of BGC types by BAPS_level1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p1_pct_label

ggsave(
  "BGC_percentage_by_BAPS_level1_stacked_bar_with_labels.pdf",
  p1_pct_label,
  width = 10,
  height = 6,
  device = cairo_pdf
)

baps_cols <- colorRampPalette(
  ggsci::pal_npg("nrc")(10)
)(length(unique(bgc_stats$BAPS_level1)))

names(baps_cols) <- sort(unique(as.character(bgc_stats$BAPS_level1)))

p2 <- ggplot(
  bgc_stats,
  aes(
    x = factor(BAPS_level1),
    y = mean_count,
    fill = factor(BAPS_level1)
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.15,
    width = 0.75
  ) +
  facet_wrap(~ BGC_type, scales = "free_y") +
  scale_fill_manual(values = baps_cols) +
  theme_bgc +
  labs(
    x = "BAPS_level1",
    y = "Mean BGC count per genome",
    fill = "BAPS_level1",
    title = "Mean BGC type count per genome by BAPS_level1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p2

ggsave(
  "BGC_mean_count_by_BAPS_level1_facet_bar.pdf",
  p2,
  width = 14,
  height = 10,
  device = cairo_pdf
)

ggsave(
  "BGC_mean_count_by_BAPS_level1_facet_bar.png",
  p2,
  width = 14,
  height = 10,
  dpi = 600
)

p3 <- ggplot(
  bgc_stats,
  aes(
    x = factor(BAPS_level1),
    y = BGC_type,
    fill = mean_count
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradientn(
    colors = ggsci::pal_npg("nrc")(10),
    name = "Mean count"
  ) +
  theme_bgc +
  labs(
    x = "BAPS_level1",
    y = "BGC type",
    title = "Mean BGC type counts by BAPS_level1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p3

ggsave(
  "BGC_mean_count_by_BAPS_level1_heatmap.pdf",
  p3,
  width = 8,
  height = 7,
  device = cairo_pdf
)

ggsave(
  "BGC_mean_count_by_BAPS_level1_heatmap.png",
  p3,
  width = 8,
  height = 7,
  dpi = 600
)