library(tidyverse)
library(DESeq2)
library(WebGestaltR)
library(here)

<<<<<<< HEAD
PADJ_THRES = 100
LFC_THRES = 0
MIN_NUM_GENES_IN_GO_CATEGORY = 5
MAX_NUM_GENES_IN_GO_CATEGORY = 2000
=======
PADJ_THRES = .05
LFC_THRES = 1
MIN_NUM_GENES_IN_GO_CATEGORY = 5
MAX_NUM_GENES_IN_GO_CATEGORY = 500
>>>>>>> 8fa902069a054d712991ad24f35a0a2007b11ed5

PROJECT_NAME = "mimic_cc"
ENRICH_DB = listGeneSet(organism = "scerevisiae")[c(2,4,6,7,8,9,10,11,12), 'name']

res_list = readRDS(here("data/shrunken_res_lists.rds"))

webgestalt_input = here("data/webgestalt_input")

dir.create(here(webgestalt_input))

for(res_list_name in names(res_list)){

  tp_res = res_list[[res_list_name]]

  for(geno in names(tp_res)){

    geno_at_tp_res = tp_res[[geno]] %>%
      as_tibble(rownames = 'gene') %>%
      filter(padj < PADJ_THRES, abs(log2FoldChange) > LFC_THRES) %>%
      select(gene, log2FoldChange)

    filename = paste0(geno, ".rnk")

    tp_dir = file.path(webgestalt_input, res_list_name)
    dir.create(tp_dir)

    outpath = file.path(tp_dir, filename)

     if(nrow(geno_at_tp_res) > 20){
       write_tsv(geno_at_tp_res, outpath, col_names = FALSE)
     }


  }
}

# note: ran into errors form webgestalt. posted issue. removed some results
# tables by hand -- lys14 from minus lys, rgt1 and lys14 from plus lys
# continues to cause errors with interaction terms

for(res_list_name in names(res_list)){

  message(res_list_name)

  interest_gene_folder =  file.path(webgestalt_input, res_list_name)

  output_dir = file.path(here("data/webgestalt_out"), res_list_name)
  dir.create(output_dir, recursive = TRUE)

  WebGestaltRBatch(
    interestGeneFolder = interest_gene_folder,
    enrichMethod = "GSEA",
    organism = "scerevisiae",
    interestGeneType="ensembl_gene_id",
<<<<<<< HEAD
    isParallel = FALSE,
    nThreads = 1,
=======
    isParallel = TRUE,
    nThreads = 10,
>>>>>>> 8fa902069a054d712991ad24f35a0a2007b11ed5
    enrichDatabase = ENRICH_DB,
    collapseMethod = "mean",
    minNum = MIN_NUM_GENES_IN_GO_CATEGORY,
    maxNum = MAX_NUM_GENES_IN_GO_CATEGORY,
    sigMethod = "fdr",
    fdrMethod = "BH",
    fdrThr = 0.05,
    topThr = 10,
    reportNum = 20,
    perNum = 1000,
    gseaP = 1,
    isOutput = TRUE,
    outputDirectory = output_dir,
    projectName = PROJECT_NAME,
    dagColor = "continuous",
    saveRawGseaResult = FALSE,
    gseaPlotFormat = 'png',
    setCoverNum = 10
  )
}

#
# ####### GSEA example #########
# rankFile <- system.file("extdata", "GeneRankList.rnk", package="WebGestaltR")
# outputDirectory <- getwd()
# enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
#                             enrichDatabase="pathway_KEGG", interestGeneFile=rankFile,
#                             interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
#                             outputDirectory=outputDirectory)
#
# ####### NTA example #########
# enrichResult <- WebGestaltR(enrichMethod="NTA", organism="hsapiens",
#                             enrichDatabase="network_PPI_BIOGRID", interestGeneFile=geneFile,
#                             interestGeneType="genesymbol", sigMethod="top", topThr=10,
#                             outputDirectory=getwd(), highlightSeedNum=10,
#                             networkConstructionMethod="Network_Retrieval_Prioritization")
