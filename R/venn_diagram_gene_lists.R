library(tidyverse)
library(DESeq2)
library(here)

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

extractGenes = function(res, direction, lfc_thres, padj_thres=.05){

  switch (direction,
          'up' =   res %>%
            as_tibble(rownames = 'gene') %>%
            filter(log2FoldChange > lfc_thres, padj < padj_thres) %>%
            pull(gene),
          'down' =   res %>%
            as_tibble(rownames = 'gene') %>%
            filter(log2FoldChange < lfc_thres, padj < padj_thres) %>%
            pull(gene)
  )
}

extractVennResults = function(res_list){

  out_list = list(
    up = list(
      singles = map(res_list[c(1,4,5)],
                    extractGenes, 'up', 1),
      eds1_vs_lys14 = map(res_list[c(1,2,4)],
                          extractGenes, 'up', 1),
      eds1_vs_rgt1 = map(res_list[c(1,3,5)],
                         extractGenes, 'up', 1)
    ),
    down = list(
      singles = map(res_list[c(1,4,5)],
                    extractGenes, 'down', -1),

      eds1_vs_lys14 = map(res_list[c(1,2,4)],
                          extractGenes, 'down', -1),

      eds1_vs_rgt1 = map(res_list[c(1,3,5)],
                         extractGenes, 'down', -1)
    )
  )

  names(out_list$up$singles) = names(res_list[c(1,4,5)])
  names(out_list$up$eds1_vs_lys14) = names(res_list[c(1,2,4)])
  names(out_list$up$eds1_vs_rgt1) = names(res_list[c(1,3,5)])

  names(out_list$down$singles) = names(res_list[c(1,4,5)])
  names(out_list$down$eds1_vs_lys14) = names(res_list[c(1,2,4)])
  names(out_list$down$eds1_vs_rgt1) = names(res_list[c(1,3,5)])

  out_list

}

createVennPlotList = function(compare_names, gene_lists){
  venn_diagram_list = list()

  for(compare in compare_names){

    for(dir in names(gene_lists)){

      dir_list = gene_lists[[dir]]
      compare_col = if(dir == 'up') 'Reds' else 'Blues'
      p = ggVennDiagram(dir_list[[compare]], label_alpha = 0, edge_size = 0) +
        scale_fill_distiller(palette = compare_col, direction = 1) +
        ggtitle(paste0("comparison: ", compare, "; direction: ", dir),
                subtitle = "shrunken abs(log2FC) > 1, padj < .05")

      venn_diagram_list[[dir]][[compare]] = p
    }
  }

  venn_diagram_list
}



venn = list(
  genes = list(
    minus_lys = extractVennResults(shrunken_res_lists$minus_lys),
    plus_lys = extractVennResults(shrunken_res_lists$plus_lys))
  )

venn[['plots']] = list(
  # note the compare names of up and down are the same
  minus_lys = createVennPlotList(names(venn$genes$minus_lys$up), venn$genes$minus_lys),
  plus_lys = createVennPlotList(names(venn$genes$plus_lys$up), venn$genes$plus_lys)
)


