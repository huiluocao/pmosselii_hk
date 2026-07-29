# ============================================================
# ARG class distribution by BAPS_level1
# Robust block parser + NPG palette + Times New Roman
# Individual plots + combined figure
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

arg_file  <- "data/pm_arg_sum.csv"
meta_file <- "data/phylogeny/metadata.txt"

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
  sysfonts::font_add(
    family = font_family,
    regular = font_path
  )
} else {
  message("Times New Roman font file was not found automatically. Using system font name.")
}

showtext::showtext_auto()

# -------------------------------
# 3. Function to read ARG file
#    Handles your block format:
#    class: Aminoglycosides
#    SRR28587992: 1
#    Pse_Mo6: 0 class: Beta-Lactams
#
#    Also handles wide CSV if needed.
# -------------------------------

read_arg_any <- function(arg_file) {
  
  first_lines <- readLines(arg_file, n = 5, warn = FALSE)
  first_nonempty <- first_lines[str_trim(first_lines) != ""][1]
  
  # -------------------------------
  # Case A: wide CSV format
  # -------------------------------
  if (!is.na(first_nonempty) &&
      str_detect(first_nonempty, '^"?class"?\\s*,')) {
    
    message("Detected ARG file as WIDE CSV format.")
    
    arg_wide <- readr::read_csv(arg_file, show_col_types = FALSE)
    names(arg_wide)[1] <- "ARG_class"
    
    arg_long <- arg_wide %>%
      pivot_longer(
        cols = -ARG_class,
        names_to = "tip",
        values_to = "count"
      ) %>%
      mutate(
        ARG_class = as.character(ARG_class),
        tip = as.character(tip),
        count = suppressWarnings(as.numeric(count))
      )
    
    return(arg_long)
  }
  
  # -------------------------------
  # Case B: block/text format
  # -------------------------------
  message("Detected ARG file as BLOCK/TEXT format.")
  
  arg_txt <- readr::read_file(arg_file)
  
  # Normalize line endings
  arg_txt_clean <- arg_txt %>%
    str_replace_all("\\r\\n?", "\n")
  
  # Remove accidental header text if copied into file
  arg_txt_clean <- arg_txt_clean %>%
    str_replace_all("Content of pm_arg_sum\\.csv:\\s*", "")
  
  # Critical fix:
  # Pse_Mo6: 0 class: Beta-Lactams
  # becomes
  # Pse_Mo6: 0
  # class: Beta-Lactams
  arg_txt_clean <- arg_txt_clean %>%
    str_replace_all("\\s+class:\\s*", "\nclass: ")
  
  # If sample-count entries are separated by spaces rather than newlines,
  # force new line before the next genome: count entry.
  arg_txt_clean <- arg_txt_clean %>%
    str_replace_all(
      "(?<=\\d)\\s+(?=[A-Za-z0-9_.-]+:\\s*-?\\d+)",
      "\n"
    )
  
  arg_lines <- str_split(arg_txt_clean, "\n")[[1]] %>%
    str_trim() %>%
    discard(~ .x == "")
  
  current_class <- NA_character_
  arg_list <- list()
  j <- 1
  
  for (ln in arg_lines) {
    
    # ARG class line
    if (str_detect(ln, "^class:\\s*")) {
      current_class <- ln %>%
        str_remove("^class:\\s*") %>%
        str_trim()
      next
    }
    
    # Genome-count line
    if (!is.na(current_class) &&
        str_detect(ln, "^[A-Za-z0-9_.-]+:\\s*-?[0-9]+\\s*$")) {
      
      m <- str_match(ln, "^([A-Za-z0-9_.-]+):\\s*(-?[0-9]+)\\s*$")
      
      arg_list[[j]] <- tibble(
        ARG_class = current_class,
        tip       = str_trim(m[1, 2]),
        count     = as.numeric(m[1, 3])
      )
      
      j <- j + 1
    }
  }
  
  arg_long <- bind_rows(arg_list)
  
  return(arg_long)
}

# -------------------------------
# 4. Read ARG file
# -------------------------------

arg_long <- read_arg_any(arg_file)

cat("\nARG parsing check:\n")
print(dim(arg_long))
print(head(arg_long, 20))
print(unique(arg_long$ARG_class))

if (nrow(arg_long) == 0) {
  stop("ARG parsing failed: arg_long is empty.")
}

# Remove duplicated genome/class entries if present
arg_long <- arg_long %>%
  group_by(ARG_class, tip) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nARG class count:\n")
print(arg_long %>% count(ARG_class))

# -------------------------------
# 5. Read metadata
# -------------------------------

metadata <- readr::read_tsv(
  meta_file,
  show_col_types = FALSE,
  na = c("NA", "n/a", "")
)

cat("\nMetadata columns:\n")
print(colnames(metadata))

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

metadata <- metadata %>%
  mutate(
    tip = as.character(tip),
    BAPS_level1 = as.factor(BAPS_level1)
  )

# Order BAPS clusters numerically
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

cat("\nMetadata genome number:\n")
print(n_distinct(metadata$tip))

# -------------------------------
# 6. Check matching genome names
# -------------------------------

missing_in_metadata <- anti_join(
  arg_long %>% distinct(tip),
  metadata %>% distinct(tip),
  by = "tip"
)

missing_in_arg <- anti_join(
  metadata %>% distinct(tip),
  arg_long %>% distinct(tip),
  by = "tip"
)

cat("\nARG genomes missing from metadata:\n")
print(missing_in_metadata)

cat("\nMetadata genomes missing from ARG file:\n")
print(missing_in_arg)

write_csv(missing_in_metadata, "ARG_genomes_missing_from_metadata.csv")
write_csv(missing_in_arg, "metadata_genomes_missing_from_ARG.csv")

# -------------------------------
# 7. Complete genome x ARG_class table
# -------------------------------

all_classes <- sort(unique(arg_long$ARG_class))

cat("\nARG classes detected:\n")
print(all_classes)

if (length(all_classes) == 0) {
  stop("No ARG classes detected.")
}

arg_complete <- metadata %>%
  select(tip, BAPS_level1) %>%
  crossing(ARG_class = all_classes) %>%
  left_join(arg_long, by = c("tip", "ARG_class")) %>%
  mutate(
    count = replace_na(count, 0),
    ARG_class = factor(ARG_class, levels = all_classes),
    BAPS_level1 = factor(BAPS_level1, levels = baps_order)
  )

cat("\nComplete table check:\n")
print(
  arg_complete %>%
    summarise(
      n_genomes = n_distinct(tip),
      n_ARG_classes = n_distinct(ARG_class),
      n_rows = n()
    )
)

# Expected if 107 genomes and 7 classes:
# n_rows = 107 x 7 = 749

# -------------------------------
# 8. Summarise by BAPS_level1
# -------------------------------

arg_summary <- arg_complete %>%
  group_by(BAPS_level1, ARG_class) %>%
  summarise(
    n_genomes = n_distinct(tip),
    total_ARG = sum(count, na.rm = TRUE),
    mean_ARG_per_genome = mean(count, na.rm = TRUE),
    genomes_with_ARG = sum(count > 0, na.rm = TRUE),
    prevalence_pct = 100 * genomes_with_ARG / n_genomes,
    .groups = "drop"
  )

arg_relative <- arg_summary %>%
  group_by(BAPS_level1) %>%
  mutate(
    total_ARG_in_BAPS = sum(total_ARG, na.rm = TRUE),
    percent_of_ARGs = if_else(
      total_ARG_in_BAPS > 0,
      100 * total_ARG / total_ARG_in_BAPS,
      0
    )
  ) %>%
  ungroup()

write_csv(arg_summary, "ARG_summary_by_BAPS_level1.csv")
write_csv(arg_relative, "ARG_relative_percentage_by_BAPS_level1.csv")

cat("\nARG summary by BAPS_level1:\n")
print(arg_summary)

# -------------------------------
# 9. NPG color palette
# -------------------------------

base_npg <- ggsci::pal_npg("nrc")(10)

if (length(all_classes) <= 10) {
  npg_cols <- base_npg[seq_along(all_classes)]
} else {
  npg_cols <- colorRampPalette(base_npg)(length(all_classes))
}

names(npg_cols) <- all_classes

cat("\nNPG colors used:\n")
print(npg_cols)

# -------------------------------
# 10. Common Times New Roman theme
# -------------------------------

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
# Plot A: Relative ARG composition
# ============================================================

p_relative <- ggplot(
  arg_relative,
  aes(x = BAPS_level1, y = percent_of_ARGs, fill = ARG_class)
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.8
  ) +
  scale_fill_manual(values = npg_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Relative composition of ARG classes",
    x = "BAPS level 1 cluster",
    y = "ARG class percentage",
    fill = "ARG class"
  ) +
  theme_times(13) +
  theme(
    legend.position = "right"
  )

# ============================================================
# Plot B: Total ARG counts
# ============================================================

p_total <- ggplot(
  arg_summary,
  aes(x = BAPS_level1, y = total_ARG, fill = ARG_class)
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.8
  ) +
  scale_fill_manual(values = npg_cols, drop = FALSE) +
  labs(
    title = "Total ARG counts",
    x = "BAPS level 1 cluster",
    y = "Total ARG count",
    fill = "ARG class"
  ) +
  theme_times(13) +
  theme(
    legend.position = "right"
  )

# ============================================================
# Plot C: Mean ARG count per genome
# ============================================================

p_mean <- ggplot(
  arg_summary,
  aes(x = BAPS_level1, y = ARG_class, fill = mean_ARG_per_genome)
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = round(mean_ARG_per_genome, 2)),
    family = font_family,
    size = 3.5,
    color = "black"
  ) +
  scale_fill_gradientn(
    colors = base_npg
  ) +
  labs(
    title = "Mean ARG count per genome",
    x = "BAPS level 1 cluster",
    y = "ARG class",
    fill = "Mean ARGs\ngenome⁻¹"
  ) +
  theme_times(13)

# ============================================================
# Plot D: ARG prevalence
# ============================================================

p_prev <- ggplot(
  arg_summary,
  aes(x = BAPS_level1, y = prevalence_pct, fill = ARG_class)
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  facet_wrap(~ ARG_class, scales = "free_y") +
  scale_fill_manual(values = npg_cols, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Prevalence of ARG classes",
    x = "BAPS level 1 cluster",
    y = "Genomes carrying ARG class",
    fill = "ARG class"
  ) +
  theme_times(12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA)
  )

# -------------------------------
# 11. Print plots
# -------------------------------

print(p_relative)
print(p_total)
print(p_mean)
print(p_prev)

# -------------------------------
# 12. Save individual plots
# -------------------------------

ggsave(
  "ARG_relative_composition_by_BAPS_level1_NPG_Times.pdf",
  p_relative,
  width = 8,
  height = 5
)

ggsave(
  "ARG_relative_composition_by_BAPS_level1_NPG_Times.png",
  p_relative,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "ARG_total_counts_by_BAPS_level1_NPG_Times.pdf",
  p_total,
  width = 8,
  height = 5
)

ggsave(
  "ARG_total_counts_by_BAPS_level1_NPG_Times.png",
  p_total,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "ARG_mean_per_genome_heatmap_by_BAPS_level1_NPG_Times.pdf",
  p_mean,
  width = 7,
  height = 5
)

ggsave(
  "ARG_mean_per_genome_heatmap_by_BAPS_level1_NPG_Times.png",
  p_mean,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "ARG_prevalence_by_BAPS_level1_NPG_Times.pdf",
  p_prev,
  width = 9,
  height = 6
)

ggsave(
  "ARG_prevalence_by_BAPS_level1_NPG_Times.png",
  p_prev,
  width = 9,
  height = 6,
  dpi = 300
)

# -------------------------------
# 13. Combined multi-panel figure
# -------------------------------

combined_ARG_plot <-
  (p_relative + p_total) /
  (p_mean + p_prev) +
  plot_annotation(
    tag_levels = "A",
    title = "Distribution of ARG classes across BAPS level 1 clusters"
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

print(combined_ARG_plot)

ggsave(
  "ARG_combined_BAPS_level1_NPG_Times.pdf",
  combined_ARG_plot,
  width = 16,
  height = 12
)

ggsave(
  "ARG_combined_BAPS_level1_NPG_Times.png",
  combined_ARG_plot,
  width = 16,
  height = 12,
  dpi = 300
)

cat("\nFinished successfully.\n")
cat("Output files:\n")
cat("1. ARG_summary_by_BAPS_level1.csv\n")
cat("2. ARG_relative_percentage_by_BAPS_level1.csv\n")
cat("3. ARG_relative_composition_by_BAPS_level1_NPG_Times.pdf/png\n")
cat("4. ARG_total_counts_by_BAPS_level1_NPG_Times.pdf/png\n")
cat("5. ARG_mean_per_genome_heatmap_by_BAPS_level1_NPG_Times.pdf/png\n")
cat("6. ARG_prevalence_by_BAPS_level1_NPG_Times.pdf/png\n")
cat("7. ARG_combined_BAPS_level1_NPG_Times.pdf/png\n")