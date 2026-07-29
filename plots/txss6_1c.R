## ============================================================
## Complete code:
## T6SS heatmap + PCoA
## With:
## 1. Times New Roman font
## 2. NPG color palette
## 3. Metadata annotations: collection_source and region
## ============================================================


## ============================================================
## 0. Install/load packages
## ============================================================

packages <- c(
  "pheatmap",
  "ggplot2",
  "ggsci",
  "showtext",
  "sysfonts",
  "RColorBrewer"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(pheatmap)
library(ggplot2)
library(ggsci)
library(showtext)
library(sysfonts)
library(RColorBrewer)


## ============================================================
## 1. Set Times New Roman font
## ============================================================

font_family <- "Times New Roman"

## ---------- Option A: Windows ----------
## If you are on Windows, this should work:
if (file.exists("C:/Windows/Fonts/times.ttf")) {
  font_add(
    family = "Times New Roman",
    regular = "C:/Windows/Fonts/times.ttf",
    bold = "C:/Windows/Fonts/timesbd.ttf",
    italic = "C:/Windows/Fonts/timesi.ttf",
    bolditalic = "C:/Windows/Fonts/timesbi.ttf"
  )
  showtext_auto()
}

## ---------- Option B: fallback if Times New Roman unavailable ----------
## This uses Google font "Tinos", which is metrically similar to Times New Roman.
## Uncomment if Times New Roman does not render correctly.
# font_add_google("Tinos", "Times New Roman")
# showtext_auto()

theme_set(theme_bw(base_family = font_family))


## ============================================================
## 2. Read full TXSS matrix
## ============================================================

txss <- read.table(
  "data/txss/tables/genome_component_matrix.tsv",
  header = TRUE,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

rownames(txss) <- trimws(txss$genome)
txss$genome <- NULL

txss[] <- lapply(txss, as.numeric)


## ============================================================
## 3. Extract T6SS columns
## ============================================================

t6ss_cols <- grep("^T6SS", colnames(txss), value = TRUE)

t6ss_mat <- as.matrix(txss[, t6ss_cols, drop = FALSE])

cat("T6SS columns extracted:\n")
print(t6ss_cols)

cat("T6SS matrix dimensions:\n")
print(dim(t6ss_mat))


## ============================================================
## 4. Read metadata
## ============================================================

metadata <- read.table(
  "data/phylogeny/metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
head(metadata)

## Your metadata should contain at least:
## genome, collection_source, region

metadata$genome <- trimws(metadata$tip)
metadata$collection_source <- trimws(metadata$Source)
metadata$region <- trimws(metadata$Contenient)

rownames(metadata) <- metadata$genome
metadata_matched <- metadata[rownames(t6ss_mat), , drop = FALSE]

metadata_matched$collection_source[is.na(metadata_matched$collection_source)] <- "Unknown"
metadata_matched$region[is.na(metadata_matched$region)] <- "Unknown"

metadata$collection_source[metadata$collection_source == ""] <- NA
metadata$region[metadata$region == ""] <- NA


## ============================================================
## 5. Match metadata to matrix genomes
## ============================================================

rownames(metadata) <- metadata$genome

missing_in_metadata <- setdiff(rownames(t6ss_mat), metadata$genome)

cat("Number of genomes missing metadata:\n")
print(length(missing_in_metadata))

if (length(missing_in_metadata) > 0) {
  cat("Genomes missing metadata:\n")
  print(missing_in_metadata)
}

metadata_matched <- metadata[rownames(t6ss_mat), , drop = FALSE]

metadata_matched$collection_source[is.na(metadata_matched$collection_source)] <- "Unknown"
metadata_matched$region[is.na(metadata_matched$region)] <- "Unknown"


## ============================================================
## 6. Cluster genomes and T6SS components
## ============================================================

row_dist <- dist(t6ss_mat, method = "euclidean")
row_clust <- hclust(row_dist, method = "ward.D2")

col_dist <- dist(t(t6ss_mat), method = "euclidean")
col_clust <- hclust(col_dist, method = "ward.D2")


## ============================================================
## 7. Define number of genome clusters
## ============================================================

## Change this if you want another number of clusters
k <- 4

genome_clusters <- cutree(row_clust, k = k)


## ============================================================
## 8. Create genome annotation table
## ============================================================

ann_row <- data.frame(
  collection_source = factor(metadata_matched$collection_source),
  region = factor(metadata_matched$region),
  T6SS_cluster = factor(genome_clusters[rownames(t6ss_mat)]),
  T6SS_total_count = rowSums(t6ss_mat),
  T6SS_components_present = rowSums(t6ss_mat > 0)
)

rownames(ann_row) <- rownames(t6ss_mat)

stopifnot(all(rownames(ann_row) == rownames(t6ss_mat)))


## ============================================================
## 9. NPG color palettes
## ============================================================

## Main NPG palette
npg_cols <- pal_npg("nrc")(10)

## Continuous NPG-based heatmap palette
## White = low/zero, NPG colors = higher counts
npg_heat_cols <- colorRampPalette(
  c("white", npg_cols)
)(100)

## Function to create NPG annotation colors
make_npg_colors <- function(x) {
  x <- factor(x)
  n <- length(levels(x))
  cols <- colorRampPalette(pal_npg("nrc")(10))(n)
  names(cols) <- levels(x)
  cols
}

ann_colors <- list(
  collection_source = make_npg_colors(ann_row$collection_source),
  region = make_npg_colors(ann_row$region),
  T6SS_cluster = make_npg_colors(ann_row$T6SS_cluster)
)


## ============================================================
## 10. Plot T6SS heatmap
## ============================================================

## PDF heatmap
pheatmap(
  t6ss_mat,
  cluster_rows = row_clust,
  cluster_cols = col_clust,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  color = npg_heat_cols,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 5,
  fontsize_col = 9,
  fontfamily = font_family,
  main = paste0("T6SS components with metadata; k = ", k),
  filename = "T6SS_heatmap_Times_New_Roman_NPG.pdf",
  width = 12,
  height = 18
)

## Optional: also save as high-resolution PNG
png(
  "T6SS_heatmap_Times_New_Roman_NPG.png",
  width = 3600,
  height = 5400,
  res = 300
)

pheatmap(
  t6ss_mat,
  cluster_rows = row_clust,
  cluster_cols = col_clust,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  color = npg_heat_cols,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 5,
  fontsize_col = 9,
  fontfamily = font_family,
  main = paste0("T6SS components with metadata; k = ", k)
)

dev.off()


## ============================================================
## 11. Save T6SS matrix and annotations
## ============================================================

write.csv(
  data.frame(genome = rownames(t6ss_mat), t6ss_mat, check.names = FALSE),
  "T6SS_only_matrix.csv",
  row.names = FALSE
)

write.csv(
  data.frame(genome = rownames(ann_row), ann_row, check.names = FALSE),
  "T6SS_annotations_with_metadata.csv",
  row.names = FALSE
)

row_order <- rownames(t6ss_mat)[row_clust$order]
col_order <- colnames(t6ss_mat)[col_clust$order]

write.csv(
  data.frame(
    genome = row_order,
    t6ss_mat[row_order, col_order, drop = FALSE],
    check.names = FALSE
  ),
  "T6SS_matrix_clustered_order.csv",
  row.names = FALSE
)


## ============================================================
## 12. Create T6SS presence/absence matrix for PCoA
## ============================================================

t6ss_pa <- ifelse(t6ss_mat > 0, 1, 0)

write.csv(
  data.frame(genome = rownames(t6ss_pa), t6ss_pa, check.names = FALSE),
  "T6SS_presence_absence_matrix.csv",
  row.names = FALSE
)


## ============================================================
## 13. PCoA using binary/Jaccard distance
## ============================================================

dist_pa <- dist(t6ss_pa, method = "binary")

pcoa <- cmdscale(
  dist_pa,
  eig = TRUE,
  k = 2
)

pcoa_df <- data.frame(
  genome = rownames(t6ss_pa),
  Axis1 = pcoa$points[, 1],
  Axis2 = pcoa$points[, 2],
  ann_row,
  check.names = FALSE
)

pcoa_var <- round(
  100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]),
  1
)

write.csv(
  pcoa_df,
  "T6SS_PCoA_coordinates_with_metadata.csv",
  row.names = FALSE
)


## ============================================================
## 14. PCoA plot: color by T6SS cluster, shape by region
## ============================================================

p_pcoa_cluster_region <- ggplot(
  pcoa_df,
  aes(
    x = Axis1,
    y = Axis2,
    color = T6SS_cluster,
    shape = region
  )
) +
  geom_point(
    size = 3.8,
    alpha = 0.9
  ) +
  scale_color_npg(name = "T6SS cluster") +
  theme_bw(base_family = font_family) +
  theme(
    text = element_text(family = font_family, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  ) +
  xlab(paste0("PCoA1 (", pcoa_var[1], "%)")) +
  ylab(paste0("PCoA2 (", pcoa_var[2], "%)")) +
  ggtitle("PCoA of T6SS presence/absence profiles")

print(p_pcoa_cluster_region)

ggsave(
  "T6SS_PCoA_cluster_region_Times_New_Roman_NPG.pdf",
  plot = p_pcoa_cluster_region,
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  "T6SS_PCoA_cluster_region_Times_New_Roman_NPG.png",
  plot = p_pcoa_cluster_region,
  width = 8,
  height = 6,
  dpi = 300
)


## ============================================================
## 15. Alternative PCoA plot: color by region
## ============================================================

p_pcoa_region <- ggplot(
  pcoa_df,
  aes(
    x = Axis1,
    y = Axis2,
    color = region
  )
) +
  geom_point(
    size = 3.8,
    alpha = 0.9
  ) +
  scale_color_npg(name = "Region") +
  theme_bw(base_family = font_family) +
  theme(
    text = element_text(family = font_family, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  ) +
  xlab(paste0("PCoA1 (", pcoa_var[1], "%)")) +
  ylab(paste0("PCoA2 (", pcoa_var[2], "%)")) +
  ggtitle("PCoA of T6SS profiles by region")

ggsave(
  "T6SS_PCoA_region_Times_New_Roman_NPG.pdf",
  plot = p_pcoa_region,
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  "T6SS_PCoA_region_Times_New_Roman_NPG.png",
  plot = p_pcoa_region,
  width = 8,
  height = 6,
  dpi = 300
)


## ============================================================
## 16. Alternative PCoA plot: color by collection source
## ============================================================

p_pcoa_source <- ggplot(
  pcoa_df,
  aes(
    x = Axis1,
    y = Axis2,
    color = collection_source
  )
) +
  geom_point(
    size = 3.8,
    alpha = 0.9
  ) +
  scale_color_npg(name = "Collection source") +
  theme_bw(base_family = font_family) +
  theme(
    text = element_text(family = font_family, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  ) +
  xlab(paste0("PCoA1 (", pcoa_var[1], "%)")) +
  ylab(paste0("PCoA2 (", pcoa_var[2], "%)")) +
  ggtitle("PCoA of T6SS profiles by collection source")

ggsave(
  "T6SS_PCoA_collection_source_Times_New_Roman_NPG.pdf",
  plot = p_pcoa_source,
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  "T6SS_PCoA_collection_source_Times_New_Roman_NPG.png",
  plot = p_pcoa_source,
  width = 8,
  height = 6,
  dpi = 300
)


## ============================================================
## 17. Optional: filled-point PCoA, publication style
## ============================================================

p_pcoa_filled <- ggplot(
  pcoa_df,
  aes(
    x = Axis1,
    y = Axis2,
    fill = T6SS_cluster,
    shape = region
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9,
    color = "black",
    stroke = 0.4
  ) +
  scale_fill_npg(name = "T6SS cluster") +
  scale_shape_manual(
    values = rep(c(21, 22, 23, 24, 25), length.out = length(levels(pcoa_df$region))),
    name = "Region"
  ) +
  theme_bw(base_family = font_family) +
  theme(
    text = element_text(family = font_family, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  ) +
  xlab(paste0("PCoA1 (", pcoa_var[1], "%)")) +
  ylab(paste0("PCoA2 (", pcoa_var[2], "%)")) +
  ggtitle("PCoA of T6SS profiles")

ggsave(
  "T6SS_PCoA_filled_points_Times_New_Roman_NPG.pdf",
  plot = p_pcoa_filled,
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  "T6SS_PCoA_filled_points_Times_New_Roman_NPG.png",
  plot = p_pcoa_filled,
  width = 8,
  height = 6,
  dpi = 300
)


## ============================================================
## Done
## ============================================================

cat("\nFinished.\n")
cat("Main files generated:\n")
cat("1. T6SS_heatmap_Times_New_Roman_NPG.pdf\n")
cat("2. T6SS_heatmap_Times_New_Roman_NPG.png\n")
cat("3. T6SS_PCoA_cluster_region_Times_New_Roman_NPG.pdf\n")
cat("4. T6SS_PCoA_cluster_region_Times_New_Roman_NPG.png\n")
cat("5. T6SS_PCoA_coordinates_with_metadata.csv\n")
cat("6. T6SS_annotations_with_metadata.csv\n")