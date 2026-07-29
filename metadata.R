# Read the metadata file
# Assuming the file is tab-separated
metadata <- read.delim(
  "data/phylogeny/metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Load required libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(knitr)
library(gridExtra)
library(grid)  # For textGrob and gpar

# Set global theme to use Times New Roman
theme_set(theme_minimal(base_family = "Times New Roman"))

# Read the metadata file
#metadata <- read.delim("metadata.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Clean up column names
colnames(metadata) <- trimws(colnames(metadata))

# Fix typo in column name
colnames(metadata)[which(colnames(metadata) == "Contenient")] <- "Continent"

# Display basic information
cat("========================================\n")
cat("DATASET OVERVIEW\n")
cat("========================================\n")
cat("Total records:", nrow(metadata), "\n")
cat("Columns:", paste(colnames(metadata), collapse = ", "), "\n\n")
cat("First 5 rows:\n")
print(head(metadata, 5))

# =============================================
# 1. SOURCE STATISTICS
# =============================================
cat("\n========================================\n")
cat("SOURCE STATISTICS\n")
cat("========================================\n")

# Clean Source2 - replace 'n/a' with NA
metadata$Source2_clean <- ifelse(metadata$Source2 == "n/a", NA, metadata$Source2)

source_stats <- metadata %>%
  filter(!is.na(Source2_clean)) %>%
  count(Source2_clean, name = "Count") %>%
  mutate(Percentage = round(Count / nrow(metadata) * 100, 1)) %>%
  arrange(desc(Count))

cat("\nSource Distribution:\n")
print(source_stats)

# =============================================
# 2. CONTINENT STATISTICS
# =============================================
cat("\n========================================\n")
cat("CONTINENT STATISTICS\n")
cat("========================================\n")

continent_stats <- metadata %>%
  count(Continent, name = "Count") %>%
  mutate(Percentage = round(Count / nrow(metadata) * 100, 1)) %>%
  arrange(desc(Count))

cat("\nContinent Distribution:\n")
print(continent_stats)

# =============================================
# 3. COUNTRY STATISTICS
# =============================================
cat("\n========================================\n")
cat("COUNTRY STATISTICS (Top 10)\n")
cat("========================================\n")

country_stats <- metadata %>%
  count(Country, name = "Count") %>%
  mutate(Percentage = round(Count / nrow(metadata) * 100, 1)) %>%
  arrange(desc(Count)) %>%
  head(10)

cat("\nTop 10 Countries:\n")
print(country_stats)

# =============================================
# 4. SOURCE TYPE STATISTICS
# =============================================
cat("\n========================================\n")
cat("SOURCE TYPE STATISTICS\n")
cat("========================================\n")

source_type_stats <- metadata %>%
  count(Source, name = "Count") %>%
  mutate(Percentage = round(Count / nrow(metadata) * 100, 1)) %>%
  arrange(desc(Count))

cat("\nSource Type Distribution:\n")
print(source_type_stats)

# =============================================
# 5. DATABASE STATISTICS
# =============================================
cat("\n========================================\n")
cat("DATABASE STATISTICS\n")
cat("========================================\n")

db_stats <- metadata %>%
  count(database, name = "Count") %>%
  mutate(Percentage = round(Count / nrow(metadata) * 100, 1)) %>%
  arrange(desc(Count))

cat("\nDatabase Distribution:\n")
print(db_stats)

# =============================================
# 6. YEAR GROUP STATISTICS
# =============================================
cat("\n========================================\n")
cat("YEAR GROUP STATISTICS\n")
cat("========================================\n")

year_stats <- metadata %>%
  count(year_group, name = "Count") %>%
  mutate(Percentage = round(Count / nrow(metadata) * 100, 1)) %>%
  arrange(year_group)

cat("\nYear Group Distribution:\n")
print(year_stats)

# =============================================
# 7. CROSS-TABULATION: Source by Continent
# =============================================
cat("\n========================================\n")
cat("SOURCE BY CONTINENT (Top 5 sources per continent)\n")
cat("========================================\n")

source_continent_cross <- metadata %>%
  filter(!is.na(Source2_clean)) %>%
  count(Continent, Source2_clean) %>%
  group_by(Continent) %>%
  mutate(Percentage = round(n / sum(n) * 100, 1)) %>%
  arrange(Continent, desc(Percentage)) %>%
  group_by(Continent) %>%
  slice_head(n = 5) %>%
  ungroup()

print(source_continent_cross)

# =============================================
# 8. SUMMARY STATISTICS
# =============================================
cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n")

summary_stats <- data.frame(
  Metric = c(
    "Total records",
    "Unique tree values",
    "Unique poppunk values",
    "Unique BAPS_level1 values",
    "Unique BAPS_level2 values",
    "Unique continents",
    "Unique countries",
    "Unique sources",
    "NCBI records",
    "This study records"
  ),
  Value = c(
    nrow(metadata),
    length(unique(metadata$tree)),
    length(unique(metadata$poppunk)),
    length(unique(metadata$BAPS_level1)),
    length(unique(metadata$BAPS_level2)),
    length(unique(metadata$Continent)),
    length(unique(metadata$Country)),
    length(unique(metadata$Source2_clean[!is.na(metadata$Source2_clean)])),
    sum(metadata$database == "NCBI"),
    sum(metadata$database == "This study")
  )
)

print(summary_stats)

# =============================================
# VISUALIZATION - WITH TIMES NEW ROMAN FONT
# =============================================

# Create color palettes with enough colors
get_color_palette <- function(n) {
  if (n <= 8) {
    return(brewer.pal(n, "Set3"))
  } else {
    return(colorRampPalette(brewer.pal(8, "Set3"))(n))
  }
}

# Define common theme for all plots
common_theme <- theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 14, face = "bold", family = "Times New Roman", hjust = 0.5),
    axis.title = element_text(size = 12, family = "Times New Roman"),
    axis.text = element_text(size = 10, family = "Times New Roman"),
    legend.text = element_text(size = 10, family = "Times New Roman"),
    legend.title = element_text(size = 11, family = "Times New Roman", face = "bold"),
    strip.text = element_text(size = 11, family = "Times New Roman", face = "bold")
  )

# 1. Source Distribution (Top 10)
p1 <- source_stats %>%
  head(10) %>%
  ggplot(aes(x = reorder(Source2_clean, Count), y = Count, fill = Source2_clean)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, fontface = "bold", family = "Times New Roman") +
  coord_flip() +
  scale_fill_manual(values = get_color_palette(nrow(head(source_stats, 10)))) +
  labs(title = "Top 10 Sources",
       x = "Source",
       y = "Count") +
  common_theme +
  theme(legend.position = "none")

# 2. Continent Distribution (Pie Chart)
p2 <- ggplot(continent_stats, aes(x = "", y = Count, fill = Continent)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold", family = "Times New Roman") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Continent Distribution") +
  theme_void(base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 14, face = "bold", family = "Times New Roman", hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 11, family = "Times New Roman", face = "bold"))

# 3. Top 10 Countries
p3 <- country_stats %>%
  ggplot(aes(x = reorder(Country, Count), y = Count, fill = Country)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, fontface = "bold", family = "Times New Roman") +
  coord_flip() +
  scale_fill_manual(values = get_color_palette(nrow(country_stats))) +
  labs(title = "Top 10 Countries",
       x = "Country",
       y = "Count") +
  common_theme +
  theme(legend.position = "none")

# 4. Source Type Distribution
p4 <- source_type_stats %>%
  ggplot(aes(x = reorder(Source, Count), y = Count, fill = Source)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 4, fontface = "bold", family = "Times New Roman") +
  coord_flip() +
  scale_fill_manual(values = get_color_palette(nrow(source_type_stats))) +
  labs(title = "Source Type Distribution",
       x = "Source Type",
       y = "Count") +
  common_theme +
  theme(legend.position = "none")

# 5. Database Distribution
p5 <- db_stats %>%
  ggplot(aes(x = "", y = Count, fill = database)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold", family = "Times New Roman") +
  scale_fill_manual(values = c("#66b3ff", "#ff9999")) +
  labs(title = "Database Distribution") +
  theme_void(base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 14, face = "bold", family = "Times New Roman", hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 11, family = "Times New Roman", face = "bold"))

# 6. Year Group Distribution
p6 <- year_stats %>%
  filter(!is.na(year_group)) %>%
  ggplot(aes(x = year_group, y = Count, fill = year_group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold", family = "Times New Roman") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Year Group Distribution",
       x = "Year Group",
       y = "Count") +
  common_theme +
  theme(legend.position = "none")

# 7. Source by Year Group (Top 5 sources)
top5_sources <- source_stats$Source2_clean[1:5]
p7 <- metadata %>%
  filter(Source2_clean %in% top5_sources) %>%
  count(year_group, Source2_clean) %>%
  ggplot(aes(x = year_group, y = n, fill = Source2_clean)) +
  geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.8) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold", family = "Times New Roman") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Top 5 Sources by Year Group",
       x = "Year Group",
       y = "Count",
       fill = "Source") +
  common_theme +
  theme(legend.position = "bottom")

# 8. Source Type by Continent
p8 <- metadata %>%
  count(Continent, Source) %>%
  ggplot(aes(x = Continent, y = n, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5, fontface = "bold", family = "Times New Roman") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Source Type by Continent",
       x = "Continent",
       y = "Count",
       fill = "Source Type") +
  common_theme +
  theme(legend.position = "bottom")

# Arrange all plots using grid.arrange with Times New Roman title
combined_plot <- grid.arrange(
  p1, p3, p4, p6,
  p2, p5, p7, p8,
  ncol = 2,
  nrow = 4,
  top = textGrob("Metadata Statistics for Poppunk/BAPS Analysis", 
                 gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Times New Roman"))
)

# Save the combined plot
ggsave("metadata_statistics_combined.png", combined_plot, width = 16, height = 20, dpi = 300)

# Save individual plots as a multi-page PDF
pdf("metadata_statistics.pdf", width = 10, height = 8, family = "Times New Roman")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)

dev.off()

# =============================================
# EXPORT STATISTICS TO CSV
# =============================================

write.csv(source_stats, "source_statistics.csv", row.names = FALSE)
write.csv(continent_stats, "continent_statistics.csv", row.names = FALSE)
write.csv(country_stats, "country_statistics.csv", row.names = FALSE)
write.csv(source_type_stats, "source_type_statistics.csv", row.names = FALSE)
write.csv(db_stats, "database_statistics.csv", row.names = FALSE)
write.csv(year_stats, "year_statistics.csv", row.names = FALSE)
write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)
write.csv(source_continent_cross, "source_continent_cross.csv", row.names = FALSE)

cat("\n========================================\n")
cat("STATISTICS EXPORTED TO CSV FILES:\n")
cat("- summary_statistics.csv\n")
cat("- source_statistics.csv\n")
cat("- continent_statistics.csv\n")
cat("- country_statistics.csv\n")
cat("- source_type_statistics.csv\n")
cat("- database_statistics.csv\n")
cat("- year_statistics.csv\n")
cat("- source_continent_cross.csv\n")
cat("\nPLOTS EXPORTED:\n")
cat("- metadata_statistics_combined.png (8-panel combined plot)\n")
cat("- metadata_statistics.pdf (individual plots)\n")
cat("========================================\n")

# =============================================
# ADDITIONAL CUSTOM PLOTS WITH TIMES NEW ROMAN
# =============================================

# Heatmap of Source vs Continent
heatmap_data <- metadata %>%
  filter(!is.na(Source2_clean)) %>%
  count(Continent, Source2_clean) %>%
  complete(Continent, Source2_clean, fill = list(n = 0))

heatmap_plot <- heatmap_data %>%
  ggplot(aes(x = Continent, y = Source2_clean, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Count") +
  geom_text(aes(label = n), size = 3, family = "Times New Roman") +
  labs(title = "Source by Continent Heatmap",
       x = "Continent",
       y = "Source") +
  theme_minimal(base_family = "Times New Roman") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman"),
        axis.text.y = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 14, face = "bold", family = "Times New Roman", hjust = 0.5),
        axis.title = element_text(size = 12, family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        legend.title = element_text(family = "Times New Roman", face = "bold"))

ggsave("figs/source_continent_heatmap.png", heatmap_plot, width = 10, height = 8, dpi = 300)

# Network-like plot for Source, Country, Continent relationships
network_plot <- metadata %>%
  filter(!is.na(Source2_clean)) %>%
  count(Country, Source2_clean, Continent) %>%
  mutate(Source_Country = paste0(Source2_clean, " - ", Country)) %>%
  top_n(20, n) %>%
  ggplot(aes(x = reorder(Source_Country, n), y = n, fill = Continent)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Top 20 Source-Country Combinations",
       x = "Source - Country",
       y = "Count",
       fill = "Continent") +
  theme_minimal(base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 14, face = "bold", family = "Times New Roman", hjust = 0.5),
        axis.title = element_text(size = 12, family = "Times New Roman"),
        axis.text = element_text(size = 8, family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        legend.title = element_text(family = "Times New Roman", face = "bold"))

ggsave("figs/source_country_network.png", network_plot, width = 10, height = 8, dpi = 300)

cat("\nADDITIONAL PLOTS EXPORTED:\n")
cat("- source_continent_heatmap.png\n")
cat("- source_country_network.png\n")
cat("========================================\n")