# ============================================================
# VF class distribution by BAPS_level1
# Corrected final script
# Input VF format: wide CSV with first column VFC_id
# NPG palette + Times New Roman
# Includes Adherence check and corrected relative stacked bar
# ============================================================

rm(list = ls())

# -------------------------------
# 0. Packages
# -------------------------------

packages <- c(
  "tidyverse",
  "ggsci",
  "showtext",
  "sysfonts",
  "patchwork"
)

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

if (length(to_install) > 0) {
  install.packages(to_install)
}

library(tidyverse)
library(ggsci)
library(showtext)
library(sysfonts)
library(patchwork)

# -------------------------------
# 1. File paths
# -------------------------------

vf_file   <- "data/pm_vf_vfc.csv"
meta_file <- "data/phylogeny/metadata.txt"

outdir <- "figs/VF_BAPS_outputs_final"
dir.create(outdir, showWarnings = FALSE)

# -------------------------------
# 2. Times New Roman font
# -------------------------------

font_family <- "Times New Roman"

font_paths <- c(
  "C:/Windows/Fonts/times.ttf",
  "C:/Windows/Fonts/timesbd.ttf",
  "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
  "/Library/Fonts/Times New Roman.ttf",
  "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"
)

font_path <- font_paths[file.exists(font_paths)][1]

if (!is.na(font_path)) {
  try(
    sysfonts::font_add(
      family = font_family,
      regular = font_path
    ),
    silent = TRUE
  )
} else {
  message("Times New Roman not found. Trying Tinos fallback.")
  try(
    sysfonts::font_add_google("Tinos", "Times New Roman"),
    silent = TRUE
  )
}

showtext::showtext_auto()

# ============================================================
# 3. Read VF wide CSV and convert to long format
# ============================================================

vf_wide <- readr::read_csv(
  vf_file,
  show_col_types = FALSE
)

# Rename first column safely
names(vf_wide)[1] <- "VF_class"

cat("\nVF wide table check:\n")
print(dim(vf_wide))
print(head(vf_wide[, 1:5]))

vf_long <- vf_wide %>%
  mutate(
    VF_class = as.character(VF_class),
    VF_class = stringr::str_replace_all(VF_class, '"', ""),
    VF_class = stringr::str_squish(VF_class)
  ) %>%
  pivot_longer(
    cols = -VF_class,
    names_to = "tip",
    values_to = "count"
  ) %>%
  mutate(
    tip = as.character(tip),
    tip = stringr::str_squish(tip),
    count = suppressWarnings(as.numeric(count)),
    count = replace_na(count, 0),
    VFC_code = stringr::str_match(VF_class, "\\((VFC[0-9]+)\\)")[, 2],
    VF_name = stringr::str_remove(VF_class, "\\s*\\(VFC[0-9]+\\)$"),
    VF_name = stringr::str_squish(VF_name)
  )

# Collapse duplicate genome/VF rows if present
vf_long <- vf_long %>%
  group_by(VF_class, VF_name, VFC_code, tip) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nVF long table check:\n")
print(dim(vf_long))
print(head(vf_long, 20))
print(unique(vf_long$VF_name))

cat("\nVF class count:\n")
print(vf_long %>% count(VF_class, VF_name, VFC_code))

if (nrow(vf_long) == 0) {
  stop("VF parsing failed: vf_long is empty.")
}

# ============================================================
# 4. Read metadata
# ============================================================

metadata <- readr::read_tsv(
  meta_file,
  show_col_types = FALSE,
  na = c("NA", "n/a", "")
) %>%
  mutate(
    tip = as.character(tip),
    tip = stringr::str_squish(tip),
    BAPS_level1 = as.factor(BAPS_level1)
  )

required_cols <- c("tip", "BAPS_level1")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "metadata.txt is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

cat("\nMetadata check:\n")
print(dim(metadata))
print(head(metadata))

# Order BAPS numerically
baps_order <- metadata %>%
  distinct(BAPS_level1) %>%
  mutate(
    BAPS_num = suppressWarnings(as.numeric(as.character(BAPS_level1)))
  ) %>%
  arrange(BAPS_num, BAPS_level1) %>%
  pull(BAPS_level1)

metadata <- metadata %>%
  mutate(
    BAPS_level1 = factor(BAPS_level1, levels = baps_order)
  )

# ============================================================
# 5. Matching check
# ============================================================

missing_in_metadata <- anti_join(
  vf_long %>% distinct(tip),
  metadata %>% distinct(tip),
  by = "tip"
)

missing_in_vf <- anti_join(
  metadata %>% distinct(tip),
  vf_long %>% distinct(tip),
  by = "tip"
)

cat("\nVF genomes missing from metadata:\n")
print(missing_in_metadata)

cat("\nMetadata genomes missing from VF file:\n")
print(missing_in_vf)

write_csv(
  missing_in_metadata,
  file.path(outdir, "VF_genomes_missing_from_metadata.csv")
)

write_csv(
  missing_in_vf,
  file.path(outdir, "metadata_genomes_missing_from_VF.csv")
)

# ============================================================
# 6. Complete genome x VF_class table
# ============================================================

vf_classes <- vf_long %>%
  distinct(VF_class, VF_name, VFC_code) %>%
  arrange(VF_name)

all_vfs <- vf_classes$VF_class

vf_complete <- metadata %>%
  select(tip, BAPS_level1, Source, Source2, Country, Contenient) %>%
  crossing(VF_class = all_vfs) %>%
  left_join(vf_classes, by = "VF_class") %>%
  left_join(
    vf_long,
    by = c("tip", "VF_class", "VF_name", "VFC_code")
  ) %>%
  mutate(
    count = replace_na(count, 0),
    VF_class = factor(VF_class, levels = all_vfs),
    VF_name = factor(VF_name, levels = vf_classes$VF_name),
    BAPS_level1 = factor(BAPS_level1, levels = baps_order)
  )

cat("\nComplete VF table check:\n")
print(
  vf_complete %>%
    summarise(
      n_genomes = n_distinct(tip),
      n_VF_classes = n_distinct(VF_class),
      n_rows = n()
    )
)

write_csv(
  vf_complete,
  file.path(outdir, "VF_complete_genome_by_class_table.csv")
)

# ============================================================
# 7. Summarise by BAPS_level1
# ============================================================

vf_summary <- vf_complete %>%
  group_by(BAPS_level1, VF_class, VF_name, VFC_code) %>%
  summarise(
    n_genomes = n_distinct(tip),
    total_VF = sum(count, na.rm = TRUE),
    mean_VF_per_genome = mean(count, na.rm = TRUE),
    genomes_with_VF = sum(count > 0, na.rm = TRUE),
    prevalence_pct = 100 * genomes_with_VF / n_genomes,
    .groups = "drop"
  )

vf_relative <- vf_summary %>%
  group_by(BAPS_level1) %>%
  mutate(
    total_VF_in_BAPS = sum(total_VF, na.rm = TRUE),
    percent_of_VFs = if_else(
      total_VF_in_BAPS > 0,
      100 * total_VF / total_VF_in_BAPS,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    percent_of_VFs = replace_na(percent_of_VFs, 0)
  )

write_csv(
  vf_summary,
  file.path(outdir, "VF_summary_by_BAPS_level1.csv")
)

write_csv(
  vf_relative,
  file.path(outdir, "VF_relative_percentage_by_BAPS_level1.csv")
)

# ============================================================
# 8. Adherence-specific check
# ============================================================

adherence_summary <- vf_summary %>%
  filter(
    stringr::str_detect(
      as.character(VF_name),
      regex("^Adherence$", ignore_case = TRUE)
    )
  ) %>%
  select(
    BAPS_level1,
    VF_name,
    n_genomes,
    total_VF,
    mean_VF_per_genome,
    genomes_with_VF,
    prevalence_pct
  )

cat("\nAdherence summary by BAPS_level1:\n")
print(adherence_summary)

write_csv(
  adherence_summary,
  file.path(outdir, "Adherence_summary_by_BAPS_level1_CHECK.csv")
)

adherence_relative <- vf_relative %>%
  filter(
    stringr::str_detect(
      as.character(VF_name),
      regex("^Adherence$", ignore_case = TRUE)
    )
  ) %>%
  select(
    BAPS_level1,
    VF_name,
    total_VF,
    total_VF_in_BAPS,
    percent_of_VFs
  )

cat("\nAdherence relative percentage by BAPS_level1:\n")
print(adherence_relative)

write_csv(
  adherence_relative,
  file.path(outdir, "Adherence_relative_percentage_by_BAPS_level1_CHECK.csv")
)

# Check sum of percentages per BAPS
percent_check <- vf_relative %>%
  group_by(BAPS_level1) %>%
  summarise(
    sum_percent = sum(percent_of_VFs, na.rm = TRUE),
    n_classes = n(),
    n_NA = sum(is.na(percent_of_VFs)),
    .groups = "drop"
  )

cat("\nPercentage sum check by BAPS_level1:\n")
print(percent_check)

write_csv(
  percent_check,
  file.path(outdir, "VF_percent_sum_check_by_BAPS_level1.csv")
)

# ============================================================
# 9. Colors and Times New Roman theme
# ============================================================

base_npg <- ggsci::pal_npg("nrc")(10)

if (length(all_vfs) <= 10) {
  vf_cols <- base_npg[seq_along(all_vfs)]
} else {
  vf_cols <- colorRampPalette(base_npg)(length(all_vfs))
}

names(vf_cols) <- all_vfs

vf_labels <- setNames(
  as.character(vf_classes$VF_name),
  as.character(vf_classes$VF_class)
)

write_csv(
  tibble(
    VF_class = names(vf_cols),
    VF_name = vf_labels[names(vf_cols)],
    color = unname(vf_cols)
  ),
  file.path(outdir, "VF_NPG_color_mapping.csv")
)

theme_times <- function(base_size = 13) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = "black"),
      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = base_size + 2,
        hjust = 0.5
      ),
      axis.title = element_text(
        family = font_family,
        face = "bold",
        size = base_size
      ),
      axis.text = element_text(
        family = font_family,
        color = "black",
        size = base_size - 2
      ),
      legend.title = element_text(
        family = font_family,
        face = "bold",
        size = base_size - 1
      ),
      legend.text = element_text(
        family = font_family,
        size = base_size - 3
      ),
      strip.text = element_text(
        family = font_family,
        face = "bold",
        color = "black",
        size = base_size - 2
      )
    )
}

# ============================================================
# 10. Plot A: Relative VF composition
#     IMPORTANT FIX:
#     Do NOT use limits = c(0, 100)
#     Use coord_cartesian() to avoid removing stacked segments.
# ============================================================

p_relative <- ggplot(
  vf_relative,
  aes(x = BAPS_level1, y = percent_of_VFs, fill = VF_class)
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.8
  ) +
  scale_fill_manual(
    values = vf_cols,
    labels = vf_labels,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  coord_cartesian(
    ylim = c(0, 100),
    clip = "off"
  ) +
  labs(
    title = "Relative composition of VF classes",
    x = "BAPS level 1 cluster",
    y = "VF class percentage",
    fill = "VF class"
  ) +
  theme_times(13) +
  theme(
    legend.position = "right"
  )

# ============================================================
# 11. Plot B: Total VF counts
# ============================================================

p_total <- ggplot(
  vf_summary,
  aes(x = BAPS_level1, y = total_VF, fill = VF_class)
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.8
  ) +
  scale_fill_manual(
    values = vf_cols,
    labels = vf_labels,
    drop = FALSE
  ) +
  labs(
    title = "Total VF counts",
    x = "BAPS level 1 cluster",
    y = "Total VF count",
    fill = "VF class"
  ) +
  theme_times(13) +
  theme(
    legend.position = "right"
  )

# ============================================================
# 12. Plot C: Mean VF count per genome heatmap
# ============================================================

p_mean <- ggplot(
  vf_summary,
  aes(x = BAPS_level1, y = VF_name, fill = mean_VF_per_genome)
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = round(mean_VF_per_genome, 2)),
    family = font_family,
    size = 3.2,
    color = "black"
  ) +
  scale_fill_gradientn(
    colors = base_npg
  ) +
  labs(
    title = "Mean VF count per genome",
    x = "BAPS level 1 cluster",
    y = "VF class",
    fill = "Mean VFs\ngenome⁻¹"
  ) +
  theme_times(13)

# ============================================================
# 13. Plot D: Prevalence of VF classes
# ============================================================

p_prev <- ggplot(
  vf_summary,
  aes(x = BAPS_level1, y = prevalence_pct, fill = VF_class)
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  facet_wrap(~ VF_name, scales = "free_y") +
  scale_fill_manual(
    values = vf_cols,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Prevalence of VF classes",
    x = "BAPS level 1 cluster",
    y = "Genomes carrying VF class",
    fill = "VF class"
  ) +
  theme_times(12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA)
  )

# Print plots
print(p_relative)
print(p_total)
print(p_mean)
print(p_prev)

# ============================================================
# 14. Save individual plots
# ============================================================

ggsave(
  file.path(outdir, "VF_relative_composition_by_BAPS_level1_FIXED_NPG_Times.pdf"),
  p_relative,
  width = 9,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "VF_relative_composition_by_BAPS_level1_FIXED_NPG_Times.png"),
  p_relative,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "VF_total_counts_by_BAPS_level1_NPG_Times.pdf"),
  p_total,
  width = 9,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "VF_total_counts_by_BAPS_level1_NPG_Times.png"),
  p_total,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "VF_mean_per_genome_heatmap_by_BAPS_level1_NPG_Times.pdf"),
  p_mean,
  width = 8,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "VF_mean_per_genome_heatmap_by_BAPS_level1_NPG_Times.png"),
  p_mean,
  width = 8,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "VF_prevalence_by_BAPS_level1_NPG_Times.pdf"),
  p_prev,
  width = 10,
  height = 7,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "VF_prevalence_by_BAPS_level1_NPG_Times.png"),
  p_prev,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# 15. Combined multi-panel figure
# ============================================================

combined_VF_plot <-
  (p_relative + p_total) /
  (p_mean + p_prev) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Distribution of virulence factor classes across BAPS level 1 clusters"
  ) &
  theme(
    text = element_text(family = font_family, color = "black"),
    plot.title = element_text(
      family = font_family,
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    plot.tag = element_text(
      family = font_family,
      face = "bold",
      size = 18
    )
  )

print(combined_VF_plot)

ggsave(
  file.path(outdir, "VF_combined_BAPS_level1_FIXED_NPG_Times.pdf"),
  combined_VF_plot,
  width = 17,
  height = 12,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "VF_combined_BAPS_level1_FIXED_NPG_Times.png"),
  combined_VF_plot,
  width = 17,
  height = 12,
  dpi = 300
)

# ============================================================
# 16. Finished
# ============================================================

cat("\nFinished successfully.\n")
cat("Outputs saved in: ", outdir, "\n")
cat("\nKey check file:\n")
cat(file.path(outdir, "Adherence_summary_by_BAPS_level1_CHECK.csv"), "\n")
