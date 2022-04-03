library(tidyverse)
library(purrr)
library(DESeq2)
library(here)

# # CURRENTLY THE FILTER IS SET TO abs(lfc) >.5 and padj < .05
#
# dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))
#
#
#
# extract_res_df = function(res_list, lfc_thres, padj_thres){
#
#   df_list = list()
#
#   for(res_name in names(res_list)){
#     lfc_df = res_list[[res_name]] %>%
#       as_tibble(rownames = 'gene') %>%
#       filter(abs(log2FoldChange) > lfc_thres,
#              padj < padj_thres) %>%
#       dplyr::select(gene, log2FoldChange) %>%
#       dplyr::rename(!!res_name := log2FoldChange)
#
#     lfc_df = tibble(gene = rownames(dds)) %>%
#       left_join(lfc_df)
#
#     df_list[[res_name]] = lfc_df
#   }
#
#   df_list
# }
#
# res_df = map(shrunken_res_lists, extract_res_df, .5, .05)
#
# plus_lys_df = plyr::join_all(res_df$plus_lys, type = "left", by = "gene")
# colnames(plus_lys_df)[2:ncol(plus_lys_df)] = paste0(colnames(plus_lys_df)[2:ncol(plus_lys_df)], "_plus_lys")
#
# minus_lys_df = plyr::join_all(res_df$minus_lys, type = "left", by = "gene")
# colnames(minus_lys_df)[2:ncol(minus_lys_df)] = paste0(colnames(minus_lys_df)[2:ncol(minus_lys_df)], "_minus_lys")
#
# long_res_df =
#   plus_lys_df %>% left_join(minus_lys_df, by = "gene") %>%
#   pivot_longer(-gene, names_to = "genotype", values_to = "shrunk_lfc") %>%
#   mutate(condition = ifelse(str_detect(genotype, "minus"), "minus_lys", "plus_lys")) %>%
#   mutate(genotype = str_remove_all(genotype, "plus_lys|minus_lys")) %>%
#   dplyr::rename(id = gene) %>%
#   mutate(genotype = str_remove(genotype, "_$"))
#
# x = reduce(list(
#   as_tibble(shrunken_res_lists$plus_lys$EDS1, rownames = 'id') %>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_eds1_plus_lys")),
#   as_tibble(shrunken_res_lists$minus_lys$EDS1, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_eds1_minus_lys")),
#   as_tibble(shrunken_res_lists$plus_lys$RGT1, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste1( ., "_rgt1_plus_lys")),
#   as_tibble(shrunken_res_lists$minus_lys$RGT1, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_rgt1_minus_lys")),
#   as_tibble(shrunken_res_lists$plus_lys$LYS14, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_lys14_plus_lys")),
#   as_tibble(shrunken_res_lists$minus_lys$LYS14, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_lys14_minus_lys")),
#   as_tibble(shrunken_res_lists$plus_lys$EDS1_RGT1, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_eds1_rgt1_plus_lys")),
#   as_tibble(shrunken_res_lists$minus_lys$EDS1_RGT1, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_eds1_rgt1_minus_lys")),
#   as_tibble(shrunken_res_lists$plus_lys$EDS1_LYS14, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_eds1_lys14_plus_lys")),
#   as_tibble(shrunken_res_lists$minus_lys$EDS1_LYS14, rownames = 'id')%>%
#     rename_at(vars(!starts_with("id")), ~paste0( ., "_eds1_lys14_minus_lys"))
# ), left_join, by = "id") %>%
#   pivot_longer(-id) %>%
#   mutate(metric = str_extract(name, '(.*?)\\_')) %>%
#   mutate(name = str_remove(name, '(.*?)\\_'),
#          metric = str_remove(metric, "\\_")) %>%
#   mutate(ko = str_extract(name, '(.*?)\\_plus|(.*?)\\_minus')) %>%
#   mutate(ko = str_remove(ko, "_plus|_minus")) %>%
#   mutate(lysine = ifelse(str_detect(name, "plus"), 'plus', 'minus')) %>%
#   dplyr::select(-name) %>%
#   pivot_wider(names_from = metric) %>%
#   left_join(cc_df %>%
#               pivot_longer(-c(id,symbol), values_to = 'cc_pval') %>%
#               mutate(lysine =
#                        ifelse(str_detect(name, 'plus'),
#                               'plus',
#                               'minus'),
#                      ko = str_extract(name, '(.*?)\\_')) %>%
#               mutate(ko = str_remove(ko, "_$")) %>%
#               dplyr::select(-name)) %>%
#   filter(complete.cases(.))


#'
#' add a suffix to the metric column names of the res table
#' @description add a suffix to the metric column names of the res table
#'   to facilitate joining multiple results tables
#' @param res_table a DESeq results table object
#' @param suffix a string to add to each metric column name
#'
#' @return a tibble with rownames converted to a column named gene,
#'   and the suffix added to each metric column name
#'
#'
add_suffix_to_colnames = function(res_table, suffix){

  as_tibble(res_table, rownames = 'gene') %>%
    rename_at(vars(!starts_with("gene")), ~paste0( ., "_", suffix))

}

#'
#' convert a set of results table objects to a tidy long df
#' @param res_list a named list of results objects
#' @param condition_name the column name which will hold the names of the
#'   list, eg 'genotype'
#' @param set_name a condition for this set of tables, eg treatment1
#' @param set_condition the value to enter in the set_name column
#'
#' @return a tidy long dataframe of the tables in the results list
#'
reduce_res_list = function(res_list, condition_name, set_condition){
  purrr::reduce(map(names(res_list),
             ~add_suffix_to_colnames(res_list[[.]], .)),
         left_join, by = "gene") %>%
  pivot_longer(-gene) %>%
  mutate(metric = str_extract(name, '(.*?)\\_')) %>%
  mutate(name = str_remove(name, '(.*?)\\_'),
         metric = str_remove(metric, "\\_")) %>%
  dplyr::rename(!!condition_name := name) %>%
  pivot_wider(names_from = metric) %>%
  mutate(contrast  := set_condition)
}

lysine_results_long = function(shrunken_res_lists, cc_df){

  rbind(
    reduce_res_list(shrunken_res_lists$effect_of_minus_lys_on_del, 'ko', 'effect_of_minus_lysine_on_del') %>%
          left_join(filter(cc_df)),

  reduce_res_list(shrunken_res_lists$effect_of_del_in_minus_lys, 'ko', 'effect_of_del_in_minus_lysine') %>%
    left_join(filter(cc_df)),

  reduce_res_list(shrunken_res_lists$effect_of_del_in_plus_lys, 'ko', 'effect_of_del_in_plus_lysine') %>%
    left_join(filter(cc_df))

  )
}



# shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))
#
# long_res_df = lysine_results_long(shrunken_res_lists, cc_df)

# interaction_df = plyr::join_all(res_df$compare_interaction, type = "left", by = "gene")
