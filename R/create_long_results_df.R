library(tidyverse)
library(DESeq2)

dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

extract_res_df = function(res_list, lfc_thres, padj_thres){

  df_list = list()

  for(res_name in names(res_list)){
    lfc_df = res_list[[res_name]] %>%
      as_tibble(rownames = 'gene') %>%
      # filter(abs(log2FoldChange) > lfc_thres,
      #        padj < padj_thres) %>%
      dplyr::select(gene, log2FoldChange) %>%
      dplyr::rename(!!res_name := log2FoldChange)

    lfc_df = tibble(gene = rownames(dds)) %>%
      left_join(lfc_df)

    df_list[[res_name]] = lfc_df
  }

  df_list
}

res_df = map(shrunken_res_lists, extract_res_df, 1, .05)

plus_lys_df = plyr::join_all(res_df$plus_lys, type = "left", by = "gene")
colnames(plus_lys_df)[2:ncol(plus_lys_df)] = paste0(colnames(plus_lys_df)[2:ncol(plus_lys_df)], "_plus_lys")

minus_lys_df = plyr::join_all(res_df$minus_lys, type = "left", by = "gene")
colnames(minus_lys_df)[2:ncol(minus_lys_df)] = paste0(colnames(minus_lys_df)[2:ncol(minus_lys_df)], "_minus_lys")

long_res_df =
  plus_lys_df %>% left_join(minus_lys_df, by = "gene") %>%
  pivot_longer(-gene, names_to = "genotype", values_to = "shrunk_lfc") %>%
  mutate(condition = ifelse(str_detect(genotype, "minus"), "minus_lys", "plus_lys")) %>%
  mutate(genotype = str_remove_all(genotype, "plus_lys|minus_lys")) %>%
  dplyr::rename(id = gene)



# interaction_df = plyr::join_all(res_df$compare_interaction, type = "left", by = "gene")
