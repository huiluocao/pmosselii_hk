install.packages("ggVennDiagram")
install.packages("ComplexUpset")
library(tidyverse)
library(ggVennDiagram)
library(ComplexUpset)

outdir <- "mummer_PseMo1_vs_3genomes"

query_genomes <- c(
  "GCA_016008855.1",
  "GCA_016008825.1",
  "SRR26921031"
)

# -------------------------------
# Function to read show-snps output
# -------------------------------

read_mummer_snps <- function(query) {
  
  file <- file.path(
    outdir,
    paste0("Pse_Mo1_vs_", query, ".snps.tsv")
  )
  
  lines <- readLines(file, warn = FALSE)
  
  # keep only data lines beginning with a number
  data_lines <- lines[grepl("^\\s*[0-9]+", lines)]
  
  if (length(data_lines) == 0) {
    return(tibble())
  }
  
  snp_tbl <- read.table(
    text = data_lines,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = ""
  )
  
  # show-snps -ClrT usually gives these columns
  colnames(snp_tbl)[1:min(ncol(snp_tbl), 12)] <- c(
    "ref_pos",
    "ref_base",
    "query_base",
    "query_pos",
    "buff",
    "dist",
    "ref_repeat",
    "query_repeat",
    "ref_len",
    "query_len",
    "frame",
    "tags"
  )[1:min(ncol(snp_tbl), 12)]
  
  snp_tbl <- snp_tbl %>%
    mutate(
      query_genome = query,
      ref_pos = as.integer(ref_pos),
      query_pos = as.integer(query_pos),
      variant_type = case_when(
        ref_base == "." | query_base == "." ~ "indel",
        TRUE ~ "SNP"
      ),
      variant_id_position = paste0("Pse_Mo1:", ref_pos),
      variant_id_allele = paste0(
        "Pse_Mo1:",
        ref_pos,
        ":",
        ref_base,
        ">",
        query_base
      )
    )
  
  return(as_tibble(snp_tbl))
}

# -------------------------------
# Read all SNP/indel files
# -------------------------------

all_variants <- map_dfr(query_genomes, read_mummer_snps)

cat("\nVariant table dimensions:\n")
print(dim(all_variants))

cat("\nVariant counts by genome:\n")
print(
  all_variants %>%
    count(query_genome, variant_type)
)

write_csv(
  all_variants,
  file.path(outdir, "Pse_Mo1_vs_3genomes_all_MUMmer_variants.csv")
)

pairwise_summary <- all_variants %>%
  group_by(query_genome) %>%
  summarise(
    total_variants = n(),
    SNPs = sum(variant_type == "SNP"),
    indels = sum(variant_type == "indel"),
    .groups = "drop"
  )

print(pairwise_summary)

write_csv(
  pairwise_summary,
  file.path(outdir, "Pse_Mo1_vs_3genomes_pairwise_variant_summary.csv")
)


shared_all_three <- all_variants %>%
  distinct(query_genome, variant_id_allele, .keep_all = TRUE) %>%
  group_by(variant_id_allele) %>%
  summarise(
    n_genomes = n_distinct(query_genome),
    genomes = paste(sort(unique(query_genome)), collapse = ";"),
    ref_pos = first(ref_pos),
    ref_base = first(ref_base),
    query_base = first(query_base),
    variant_type = first(variant_type),
    .groups = "drop"
  ) %>%
  filter(n_genomes == length(query_genomes)) %>%
  arrange(ref_pos)

print(shared_all_three)

write_csv(
  shared_all_three,
  file.path(outdir, "Variants_shared_by_all_3_genomes_vs_Pse_Mo1.csv")
)

variant_membership <- all_variants %>%
  distinct(query_genome, variant_id_allele) %>%
  group_by(variant_id_allele) %>%
  summarise(
    n_genomes = n_distinct(query_genome),
    genomes = paste(sort(unique(query_genome)), collapse = ";"),
    .groups = "drop"
  )

unique_to_one_query <- variant_membership %>%
  filter(n_genomes == 1)

write_csv(
  unique_to_one_query,
  file.path(outdir, "Variants_unique_to_one_comparison_genome.csv")
)

for (g in query_genomes) {
  
  unique_g <- unique_to_one_query %>%
    filter(genomes == g) %>%
    left_join(
      all_variants,
      by = "variant_id_allele"
    ) %>%
    filter(query_genome == g)
  
  write_csv(
    unique_g,
    file.path(outdir, paste0("Variants_unique_to_", g, "_vs_Pse_Mo1.csv"))
  )
}

variant_sets <- all_variants %>%
  distinct(query_genome, variant_id_allele) %>%
  group_by(query_genome) %>%
  summarise(
    variants = list(variant_id_allele),
    .groups = "drop"
  ) %>%
  deframe()

p_venn <- ggVennDiagram(variant_sets) +
  scale_fill_gradient(low = "white", high = "#E64B35FF") +
  theme_void() +
  ggtitle("Shared and unique variants relative to Pse_Mo1")

print(p_venn)

ggsave(
  file.path(outdir, "Pse_Mo1_vs_3genomes_variant_Venn.pdf"),
  p_venn,
  width = 6,
  height = 5
)

ggsave(
  file.path(outdir, "Pse_Mo1_vs_3genomes_variant_Venn.png"),
  p_venn,
  width = 6,
  height = 5,
  dpi = 300
)

variant_binary <- all_variants %>%
  distinct(variant_id_allele, query_genome) %>%
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from = query_genome,
    values_from = value,
    values_fill = FALSE
  )

p_upset <- ComplexUpset::upset(
  variant_binary,
  intersect = query_genomes,
  name = "Variant set",
  width_ratio = 0.2
) +
  ggtitle("Variant intersections relative to Pse_Mo1")

print(p_upset)

ggsave(
  file.path(outdir, "Pse_Mo1_vs_3genomes_variant_UpSet.pdf"),
  p_upset,
  width = 8,
  height = 5
)

ggsave(
  file.path(outdir, "Pse_Mo1_vs_3genomes_variant_UpSet.png"),
  p_upset,
  width = 8,
  height = 5,
  dpi = 300
)

p_counts <- all_variants %>%
  count(query_genome, variant_type) %>%
  ggplot(
    aes(
      x = query_genome,
      y = n,
      fill = variant_type
    )
  ) +
  geom_col(
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "SNP" = "#4DBBD5FF",
      "indel" = "#E64B35FF"
    )
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Comparison genome",
    y = "Number of variants relative to Pse_Mo1",
    fill = "Variant type",
    title = "Pairwise SNP/indel differences relative to Pse_Mo1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_counts)

ggsave(
  file.path(outdir, "Pse_Mo1_vs_3genomes_variant_count_barplot.pdf"),
  p_counts,
  width = 7,
  height = 5
)

ggsave(
  file.path(outdir, "Pse_Mo1_vs_3genomes_variant_count_barplot.png"),
  p_counts,
  width = 7,
  height = 5,
  dpi = 300
)


read_mummer_diff <- function(query) {
  
  file <- file.path(
    outdir,
    paste0("Pse_Mo1_vs_", query, ".ref_differences.tsv")
  )
  
  if (!file.exists(file)) {
    return(tibble())
  }
  
  lines <- readLines(file, warn = FALSE)
  lines <- lines[lines != ""]
  
  if (length(lines) == 0) {
    return(tibble())
  }
  
  diff_tbl <- read.table(
    text = lines,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = ""
  )
  
  diff_tbl <- as_tibble(diff_tbl)
  colnames(diff_tbl)[1] <- "diff_type"
  
  diff_tbl <- diff_tbl %>%
    mutate(query_genome = query)
  
  return(diff_tbl)
}

all_diffs <- map_dfr(query_genomes, read_mummer_diff)

write_csv(
  all_diffs,
  file.path(outdir, "Pse_Mo1_vs_3genomes_all_showdiff_regions.csv")
)

diff_summary <- all_diffs %>%
  count(query_genome, diff_type)

print(diff_summary)

write_csv(
  diff_summary,
  file.path(outdir, "Pse_Mo1_vs_3genomes_showdiff_summary.csv")
)
