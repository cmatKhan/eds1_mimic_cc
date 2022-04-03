library(tidyverse)
library(here)

parse_cc_data = function(cc_path, cc_thres){

  cc_df_path = cc_path

  read_csv(cc_df_path) %>%
    pivot_longer(-c(id,symbol), names_to = 'tmp', values_to = 'cc_pval') %>%
    mutate(ko = str_extract(tmp, "^\\w+\\d+"),
           lysine = ifelse(str_detect(tmp, "plus"), 'plus', 'minus')) %>%
    dplyr::select(-tmp) %>%
    mutate(ko = toupper(ko)) %>%
    dplyr::rename(gene = id) %>%
    group_by(ko, gene) %>%
    summarize(min_pval = min(cc_pval), .groups = 'keep') %>%
    filter(min_pval < cc_thres)

}
