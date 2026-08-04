library(tidyverse)

snippy_dir <- "snippy_PseMo1_vs_3genomes"

query_genomes <- c(
  "GCA_016008855.1",
  "GCA_016008825.1",
  "SRR26921031"
)

read_snippy_tab <- function(q) {
  
  file <- file.path(
    snippy_dir,
    paste0(q, "_vs_Pse_Mo1"),
    "snps.tab"
  )
  
  if (!file.exists(file)) {
    warning("File not found: ", file)
    return(tibble())
  }
  
  readr::read_tsv(file, show_col_types = FALSE) %>%
    mutate(
      query_genome = q,
      variant_id = paste(CHROM, POS, REF, ALT, sep = ":")
    )
}

snippy_all <- map_dfr(query_genomes, read_snippy_tab)

write_csv(
  snippy_all,
  file.path(snippy_dir, "Pse_Mo1_vs_3genomes_annotated_snps_indels.csv")
)

# Summary by genome and mutation effect
effect_summary <- snippy_all %>%
  count(query_genome, TYPE, EFFECT, sort = TRUE)

write_csv(
  effect_summary,
  file.path(snippy_dir, "Pse_Mo1_vs_3genomes_variant_effect_summary.csv")
)

print(effect_summary)


shared_all_three <- snippy_all %>%
  distinct(query_genome, variant_id, .keep_all = TRUE) %>%
  group_by(variant_id) %>%
  summarise(
    n_genomes = n_distinct(query_genome),
    genomes = paste(sort(unique(query_genome)), collapse = ";"),
    CHROM = first(CHROM),
    POS = first(POS),
    TYPE = first(TYPE),
    REF = first(REF),
    ALT = first(ALT),
    EFFECT = first(EFFECT),
    LOCUS_TAG = first(LOCUS_TAG),
    GENE = first(GENE),
    PRODUCT = first(PRODUCT),
    .groups = "drop"
  ) %>%
  filter(n_genomes == length(query_genomes)) %>%
  arrange(CHROM, POS)

write_csv(
  shared_all_three,
  file.path(snippy_dir, "Variants_shared_by_all_3_genomes_vs_Pse_Mo1_annotated.csv")
)

print(shared_all_three)

variant_membership <- snippy_all %>%
  distinct(query_genome, variant_id) %>%
  group_by(variant_id) %>%
  summarise(
    n_genomes = n_distinct(query_genome),
    genomes = paste(sort(unique(query_genome)), collapse = ";"),
    .groups = "drop"
  )

unique_to_one_query <- variant_membership %>%
  filter(n_genomes == 1)

unique_annotated <- unique_to_one_query %>%
  left_join(snippy_all, by = "variant_id")

write_csv(
  unique_annotated,
  file.path(snippy_dir, "Variants_unique_to_one_comparison_genome_annotated.csv")
)
