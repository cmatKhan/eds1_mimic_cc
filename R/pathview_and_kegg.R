library(pathview)
library(DESeq2)
library(tidyverse)
library(here)

# from scer KEGG
pathways = list(
  carbon_metabolism = '01200',
  cell_cycle = '04111',
  biosyn_aa = '01230',
  glycolysis_gluconeogensis = "00010",
  tca_cycle = "00020",
  meiosis = '04113',
  carbon_meta = '01200',
  nitrogen_metabolism = '00910',
  lysine_pathway = "00300",
  lysine_deg = '00310',
  aa_and_nuc_sugar_metablism = '00520',
  pyruvate_metabolism = '00620',
  glyoxylate_and_dicarboxylate_metabolism = '00630',
  pentose_phosphate = '00030'
)

# use pathview::pathview to vis the lysine pathway + DE data
mapKEGGpathway = function(name, res, pathway_id, lfc_thres, padj_thres = .05, species = 'sce'){


  fltr_res = res %>%
    as.data.frame() %>%
    filter(abs(log2FoldChange) > lfc_thres &
             padj < padj_thres) %>%
    dplyr::select(log2FoldChange)


  pathview(
    gene.data = fltr_res,
    gene.idtype = 'kegg',
    # kegg.native = FALSE,
    # map.symbol = TRUE,
    expand.node = TRUE,
    pathway.id = pathway_id,
    species = species,
    out.suffix = paste0(name, "_", names(pathways[pathways == pathway_id]))
  )

}
