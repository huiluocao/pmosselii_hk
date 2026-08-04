library(tidyverse)

roary_file <- "gene_presence_absence.csv"
metadata_file <- "../phylogeny/metadata.txt"

pa <- readr::read_csv(roary_file, show_col_types = FALSE)

metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(tip = as.character(tip))
hk_isolates <- c(
  "Pse_Mo1",
  "Pse_Mo2",
  "Pse_Mo3",
  "Pse_Mo4",
  "Pse_Mo6",
  "Pse_Mo7",
  "Pse_Mo8"
)
missing_hk <- setdiff(hk_isolates, colnames(pa))
print(missing_hk)

roary_info_cols <- c(
  "Gene",
  "Non-unique Gene name",
  "Annotation",
  "No. isolates",
  "No. sequences",
  "Avg sequences per isolate",
  "Genome Fragment",
  "Order within Fragment",
  "Accessory Fragment",
  "Accessory Order with Fragment",
  "QC",
  "Min group size nuc",
  "Max group size nuc",
  "Avg group size nuc"
)

genome_cols <- setdiff(colnames(pa), roary_info_cols)

cat("Number of genome columns in Roary:\n")
print(length(genome_cols))

cat("HK isolates found in Roary:\n")
print(intersect(hk_isolates, genome_cols))

hk_cols <- intersect(hk_isolates, genome_cols)
non_hk_cols <- setdiff(genome_cols, hk_cols)

pa_binary <- pa %>%
  mutate(
    across(
      all_of(genome_cols),
      ~ if_else(!is.na(.x) & .x != "", 1, 0)
    )
  )

pa_summary <- pa_binary %>%
  mutate(
    n_HK_present = rowSums(across(all_of(hk_cols))),
    n_nonHK_present = rowSums(across(all_of(non_hk_cols))),
    prop_HK_present = n_HK_present / length(hk_cols),
    prop_nonHK_present = n_nonHK_present / length(non_hk_cols)
  )


hk_specific_strict <- pa_summary %>%
  filter(
    n_HK_present == length(hk_cols),
    n_nonHK_present == 0
  ) %>%
  select(
    Gene,
    `Non-unique Gene name`,
    Annotation,
    `No. isolates`,
    `No. sequences`,
    n_HK_present,
    n_nonHK_present,
    prop_HK_present,
    prop_nonHK_present,
    all_of(hk_cols)
  )

cat("Strict HK-specific genes:\n")
print(nrow(hk_specific_strict))

write_csv(
  hk_specific_strict,
  "HK_specific_genes_present_all_HK_absent_all_nonHK.csv"
)

hk_enriched <- pa_summary %>%
  filter(
    n_HK_present >= 5,
    prop_nonHK_present <= 0.10
  ) %>%
  arrange(desc(prop_HK_present), prop_nonHK_present) %>%
  select(
    Gene,
    `Non-unique Gene name`,
    Annotation,
    `No. isolates`,
    `No. sequences`,
    n_HK_present,
    n_nonHK_present,
    prop_HK_present,
    prop_nonHK_present,
    all_of(hk_cols)
  )

cat("HK-enriched genes:\n")
print(nrow(hk_enriched))

write_csv(
  hk_enriched,
  "HK_enriched_genes_present_5of7_HK_rare_nonHK.csv"
)

hk_any_specific <- pa_summary %>%
  filter(
    n_HK_present >= 1,
    n_nonHK_present == 0
  ) %>%
  arrange(desc(n_HK_present)) %>%
  select(
    Gene,
    `Non-unique Gene name`,
    Annotation,
    n_HK_present,
    n_nonHK_present,
    prop_HK_present,
    prop_nonHK_present,
    all_of(hk_cols)
  )

cat("Genes present in at least one HK isolate and absent from all non-HK:\n")
print(nrow(hk_any_specific))

write_csv(
  hk_any_specific,
  "Genes_present_any_HK_absent_all_nonHK.csv"
)


for (hk in hk_cols) {
  
  other_cols <- setdiff(genome_cols, hk)
  
  unique_hk <- pa_binary %>%
    mutate(
      present_in_hk = .data[[hk]],
      present_in_others = rowSums(across(all_of(other_cols)))
    ) %>%
    filter(
      present_in_hk == 1,
      present_in_others == 0
    ) %>%
    select(
      Gene,
      `Non-unique Gene name`,
      Annotation,
      `No. isolates`,
      `No. sequences`,
      all_of(hk)
    )
  
  write_csv(
    unique_hk,
    paste0(hk, "_unique_genes_absent_all_others.csv")
  )
  
  cat(hk, "unique genes:", nrow(unique_hk), "\n")
}

env_isolates <- metadata %>%
  filter(Source == "Environment") %>%
  pull(tip)

env_cols <- intersect(env_isolates, genome_cols)

cat("Environmental genomes found in Roary:\n")
print(length(env_cols))

pa_env_summary <- pa_binary %>%
  mutate(
    n_HK_present = rowSums(across(all_of(hk_cols))),
    n_ENV_present = rowSums(across(all_of(env_cols))),
    prop_HK_present = n_HK_present / length(hk_cols),
    prop_ENV_present = n_ENV_present / length(env_cols)
  )

hk_vs_env_enriched <- pa_env_summary %>%
  filter(
    n_HK_present >= 5,
    prop_ENV_present <= 0.10
  ) %>%
  arrange(desc(prop_HK_present), prop_ENV_present) %>%
  select(
    Gene,
    `Non-unique Gene name`,
    Annotation,
    n_HK_present,
    n_ENV_present,
    prop_HK_present,
    prop_ENV_present,
    all_of(hk_cols),
    all_of(env_cols)
  )

write_csv(
  hk_vs_env_enriched,
  "HK_blood_enriched_vs_environment_genes.csv"
)

hk_vs_env_strict <- pa_env_summary %>%
  filter(
    n_HK_present == length(hk_cols),
    n_ENV_present == 0
  )

write_csv(
  hk_vs_env_strict,
  "HK_blood_present_all_absent_environment.csv"
)


gene_tests <- pa_summary %>%
  rowwise() %>%
  mutate(
    fisher_p = fisher.test(
      matrix(
        c(
          n_HK_present,
          length(hk_cols) - n_HK_present,
          n_nonHK_present,
          length(non_hk_cols) - n_nonHK_present
        ),
        nrow = 2,
        byrow = TRUE
      )
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    fisher_q = p.adjust(fisher_p, method = "BH")
  ) %>%
  arrange(fisher_p) %>%
  select(
    Gene,
    `Non-unique Gene name`,
    Annotation,
    n_HK_present,
    n_nonHK_present,
    prop_HK_present,
    prop_nonHK_present,
    fisher_p,
    fisher_q
  )

write_csv(
  gene_tests,
  "HK_gene_enrichment_Fisher_tests.csv"
)

head(gene_tests, 20)


top_genes <- hk_enriched %>%
  slice_head(n = 50) %>%
  pull(Gene)

heat_df <- pa_binary %>%
  filter(Gene %in% top_genes) %>%
  select(Gene, all_of(c(hk_cols, non_hk_cols))) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Order genomes: HK first, then non-HK
heat_df <- heat_df[, c(hk_cols, non_hk_cols), drop = FALSE]

pheatmap::pheatmap(
  heat_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = c("white", "#E64B35FF"),
  border_color = "grey80",
  main = "Top HK-enriched genes",
  filename = "Top_HK_enriched_genes_heatmap.pdf",
  width = 12,
  height = 10
)

library(ComplexUpset)

hk_binary <- pa_binary %>%
  select(Gene, all_of(hk_cols))

p_upset <- ComplexUpset::upset(
  hk_binary,
  intersect = hk_cols,
  name = "Gene set",
  width_ratio = 0.2
)

ggsave(
  "HK_isolates_gene_presence_UpSet.pdf",
  p_upset,
  width = 10,
  height = 6
)
