library(DESeq2)
library(tidyverse)
library(here)

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))
cc_df = read_csv("../Combined_Eds1_Rgt1_Lys14_CCTargets.csv") %>%
  dplyr::rename(gene = id)

output_dir = here("data/dto_input")
dir.create(output_dir)

dtoDEtable = function(tf_name, res){

  as_tibble(res, rownames = 'gene') %>%
    dplyr::rename(!!quo_name(tf_name) := log2FoldChange) %>%
    dplyr::select(gene,tf_name)

}

dto_input = list(

    plus_lys_de = map(names(shrunken_res_lists$plus_lys),
                      ~dtoDEtable(., shrunken_res_lists$plus_lys[[.]])) %>%
      plyr::join_all() %>%
      setNames(tolower(colnames(.))) %>%
      select(gene, eds1, rgt1, lys14) %>%
      filter(gene %in% cc_df$gene),

    minus_lys_de = map(names(shrunken_res_lists$minus_lys),
                       ~dtoDEtable(., shrunken_res_lists$minus_lys[[.]])) %>%
      plyr::join_all() %>%
      setNames(tolower(colnames(.))) %>%
      select(gene, eds1, rgt1, lys14) %>%
      filter(gene %in% cc_df$gene)

)

dto_input$plus_lys_cc = select(cc_df, matches("gene|plus"))  %>%
  setNames(gsub("_plus_lys_pval", "", colnames(.))) %>%
  setNames(tolower(colnames(.))) %>%
  filter(gene %in% dto_input$plus_lys_de$gene)

dto_input$minus_lys_cc = select(cc_df, matches("gene|minus")) %>%
  setNames(gsub("_minus_lys_pval", "", colnames(.))) %>%
  setNames(tolower(colnames(.))) %>%
  # intentionally filter by plus lys -- ensures same genes in each (which there
  # would by if filterin with minus lys, but this is for consistency)
  filter(gene %in% dto_input$plus_lys_de$gene)

map(names(dto_input),
    ~write_csv(dto_input[[.]], file.path(output_dir, paste0(.,".csv"))))

write_csv(select(dto_input$plus_lys_de, gene), file.path(output_dir, "gene_universe.csv"), col_names = FALSE)

