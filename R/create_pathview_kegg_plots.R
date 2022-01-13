library(pathview)
library(DESeq2)
library(tidyverse)
library(here)

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))


mapYeastPathway = function(name, res, pathway_id, lfc_thres){

  fltr_res = res %>%
    as.data.frame() %>%
    filter(abs(log2FoldChange) > lfc_thres,
           padj < .05) %>%
    select(log2FoldChange)


  pathview(
    gene.data = fltr_res,
    gene.idtype = 'kegg',
    pathway.id = pathway_id,
    species = "sce",
    out.suffix = name,
    map.symbol = FALSE
  )

}



pathways = list(
  lysine_pathway = "00300",
  biosyn_aa = "01230",
  glycolysis_gluconeogensis = "00010",
  tca_cycle = "00020",
  meiosis = '04113',
  carbon_meta = '01200',
  pentose_phosphate = '00030'
)


setwd("plots/plus_lys/lysine_pathway_kegg/")

map(names(shrunken_res_lists$plus_lys),
    ~mapYeastPathway(., shrunken_res_lists$plus_lys[[.]], pathway_id = pathways$lysine_pathway, lfc_thres = .05))
