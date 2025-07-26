library(tidyverse)
library(readr)
library(rlang)


bac_file <- "bac120_metadata_r226.tsv"
arc_file <- "ar53_metadata_r226.tsv"
output_file <- "ncbi_to_gtdb_mapped_raw16.csv"
output_file <- "ncbi_to_gtdb_mapped_lit.csv"
# =====================================
allowed_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Load and merge GTDB
read_and_parse_gtdb <- function(file) {
  read_tsv(file, show_col_types = FALSE) %>%
    separate(gtdb_taxonomy,
             into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
             sep = ";\\s*", fill = "right", remove = FALSE) %>%
    separate(ncbi_taxonomy,
             into = paste0("ncbi_", c("domain", "phylum", "class", "order", "family", "genus", "species")),
             sep = ";\\s*", fill = "right", remove = FALSE) %>%
    mutate(across(all_of(c("domain", "phylum", "class", "order", "family", "genus", "species")),
                  ~str_remove(., "^[a-z]__"))) %>%
    mutate(across(starts_with("ncbi_"),
                  ~str_remove(., "^[a-z]__")))
}

gtdb_bac <- read_and_parse_gtdb(bac_file)
gtdb_arc <- read_and_parse_gtdb(arc_file)
gtdb_df <- bind_rows(gtdb_bac, gtdb_arc)



process_one <- function(ncbi_name, rank) {
  rank_col <- rank
  ncbi_rank_col <- paste0("ncbi_", rank_col)
  
  if (!(rank_col %in% allowed_ranks)) {
    return(tibble(ncbi_name = ncbi_name, rank = rank, gtdb_counts = NA_character_))
  }
  if (!(ncbi_rank_col %in% colnames(gtdb_df))) {
    return(tibble(ncbi_name = ncbi_name, rank = rank, gtdb_counts = NA_character_))
  }
  
  # Filter using {rank} of the ncbi_ broken down in the GTDB file
  subset_df <- gtdb_df %>%
    filter(.data[[ncbi_rank_col]] == ncbi_name)
  
  if (nrow(subset_df) == 0) {
    return(tibble(ncbi_name = ncbi_name, rank = rank, gtdb_counts = NA_character_))
  }
  
  # Statistics are made for names and numbers under the same rank in GTDB
  rank_sym <- sym(rank_col)
  count_table <- subset_df %>%
    group_by(name = !!rank_sym) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n)) %>%
    mutate(out = paste0(name, ":", n)) %>%
    pull(out)
  
  out_string <- paste(count_table, collapse = ";")
  tibble(ncbi_name = ncbi_name, rank = rank, gtdb_counts = out_string)
}




input_df <- read_excel('16S/Result.xlsx')
input_df <- input_df %>%
  select(Microbiota, Classification) %>%
  rename(ncbi_name = Microbiota, rank = Classification)
results_16s <- pmap_dfr(input_df, process_one)
identical(results_16s$ncbi_name,raw16_result$Microbiota)


input_df <- read_excel('Literature/Result.xlsx')
input_df <- input_df %>%
  select(Microbiota, Classification) %>%
  rename(ncbi_name = Microbiota, rank = Classification)
results_lit <- pmap_dfr(input_df, process_one)


