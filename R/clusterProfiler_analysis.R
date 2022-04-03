library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(org.Sc.sgd.db)

# https://hbctraining.github.io/DGE_workshop/lessons/09_functional_analysis.html

anno_df <- select(org.Sc.sgd.db,
               keys=rownames(dds_lrt),
               columns=c("SGD", "COMMON")) %>%
  as_tibble() %>%
  distinct(ORF, .keep_all = TRUE) %>%
  dplyr::rename(gene_name = COMMON,
                gene_id = ORF) %>%
  dplyr::select(-SGD)



clusterProfilerOutput = function(res_obj, anno_df, lfc_thres = 2){

  res_df = res_obj %>%
    as_tibble(rownames = "gene_id") %>%
    left_join(anno_df) %>%
    dplyr::filter(abs(log2FoldChange) > lfc_thres) %>%
    dplyr::filter(!is.na(gene_name))

  allOE_genes <- as.character(res_df$gene_name)
  sigOE_genes <- res_df %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::pull(gene_name)

  enrichGO(gene = sigOE_genes,
           universe = allOE_genes,
           keyType = "GENENAME",
           OrgDb = org.Sc.sgd.db,
           ont = "BP",
           pAdjustMethod = "BH",
           qvalueCutoff = 0.05,
           readable = TRUE)


}

x = clusterProfilerOutput(shrunken_res_lists$plus_lys$RGT1, anno_df)
