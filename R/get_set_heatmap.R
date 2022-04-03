library(tidyverse)


#'
#' @param long_res_df a results table in long format with columns id, shrunk_lfc,
#'   genotype, condition
#' @param pathway_genes_df a set of genes of interest with columns id and symbol.
#'   id corresponds to id column in long_res_df
#'
#'
getSetHeatMap = function(long_res_df, pathway_genes_df, select_genotypes,
                         select_conditions, title, text_size = 20){

  long_res_df %>%
    filter(id %in% pathway_genes_df$id) %>%
    mutate(shrunk_lfc = replace_na(shrunk_lfc, 0)) %>%
    filter(genotype %in% select_genotypes) %>%
    filter(condition %in% select_conditions) %>%
    left_join(pathway_genes_df) %>%
    ggplot() +
    geom_tile(aes(genotype, symbol, fill = shrunk_lfc  ), color = "black") +
    scale_fill_gradient2(low = "#075AFF",
                         mid = "#FFFFCC",
                         high = "#FF0000") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          text = element_text(size = text_size)) +
    ggtitle(title) +
    facet_grid(cols = vars(condition))
}
