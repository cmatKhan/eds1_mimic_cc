library(tidyverse)
library(DESeq2)
library(ggVennDiagram)
library(here)

extractGenes = function(res, direction, lfc_thres, padj_thres=.05){

  switch (direction,
          'up' =   res %>%
            as_tibble(rownames = 'gene') %>%
            filter(log2FoldChange > lfc_thres, padj < padj_thres) %>%
            pull(gene),
          'down' =   res %>%
            as_tibble(rownames = 'gene') %>%
            filter(log2FoldChange < -lfc_thres, padj < padj_thres) %>%
            pull(gene)
  )
}

extractVennResults = function(res_list, abs_thres, padj_thres){

  abs_thres = abs(abs_thres)

  out_list = list(
    up = list(
      singles = map(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "LYS14"))],
                    extractGenes, 'up', abs_thres, padj_thres),
      eds1_vs_lys14 = map(res_list[which(names(res_list) %in% c("EDS1", "LYS14", "EDS1_LYS14"))],
                          extractGenes, 'up', abs_thres, padj_thres),
      eds1_vs_rgt1 = map(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "EDS1_RGT1"))],
                         extractGenes, 'up', abs_thres, padj_thres)
    ),
    down = list(
      singles = map(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "LYS14"))],
                    extractGenes, 'down', abs_thres, padj_thres),

      eds1_vs_lys14 = map(res_list[which(names(res_list) %in% c("EDS1", "LYS14", "EDS1_LYS14"))],
                          extractGenes, 'down', abs_thres, padj_thres),

      eds1_vs_rgt1 = map(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "EDS1_RGT1"))],
                         extractGenes, 'down', abs_thres, padj_thres)
    )
  )

  names(out_list$up$singles) = names(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "LYS14"))])
  names(out_list$up$eds1_vs_lys14) = names(res_list[which(names(res_list) %in% c("EDS1", "LYS14", "EDS1_LYS14"))])
  names(out_list$up$eds1_vs_rgt1) = names(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "EDS1_RGT1"))])

  names(out_list$down$singles) = names(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "LYS14"))])
  names(out_list$down$eds1_vs_lys14) = names(res_list[which(names(res_list) %in% c("EDS1", "LYS14", "EDS1_LYS14"))])
  names(out_list$down$eds1_vs_rgt1) = names(res_list[which(names(res_list) %in% c("EDS1", "RGT1", "EDS1_RGT1"))])

  out_list

}

createVennPlotList = function(compare_names, gene_lists, lys_cond, abs_thres, padj_thres){
  venn_diagram_list = list()

  for(compare in compare_names){

    for(dir in names(gene_lists)){

      dir_list = gene_lists[[dir]]
      compare_col = if(dir == 'up') 'Reds' else 'Blues'
      p = ggVennDiagram(dir_list[[compare]],
                        label = "count",
                        label_alpha = 0,
                        edge_size = 0,
                        label_geom = c("label", "text"),
                        label_size = 9) +
        scale_fill_distiller(palette = compare_col, direction = 1) +
        ggtitle(paste0("shrunken abs(log2FC) > ", abs_thres," padj < ", padj_thres)) +
        theme(legend.position="none",
              plot.title = element_text(size = 15))


      venn_diagram_list[[dir]][[compare]] = p
    }
  }

  venn_diagram_list
}

getVennResults = function(shrunken_res_lists, abs_lfc_thres = 1, padj_thres = .05){

  venn = list(
    genes = list(
      minus_lys = extractVennResults(shrunken_res_lists$minus_lys, abs_lfc_thres, padj_thres),
      plus_lys = extractVennResults(shrunken_res_lists$plus_lys, abs_lfc_thres, padj_thres))
  )

  venn[['plots']] = list(
    # note the compare names of up and down are the same
    minus_lys = createVennPlotList(names(venn$genes$minus_lys$up), venn$genes$minus_lys, "Minus Lysine", abs_lfc_thres, padj_thres),
    plus_lys = createVennPlotList(names(venn$genes$plus_lys$up), venn$genes$plus_lys, "Plus Lysine", abs_lfc_thres, padj_thres)
  )

  venn

}

