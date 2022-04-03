library(DESeq2)
library(tidyverse)
library(org.Sc.sgd.db)
library(here)
library(topGO)
library(pcaExplorer)
library(org.Sc.sgd.db)

yeast_orgdb = org.Sc.sgd.db

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

dds_lrt = readRDS(here("data/mimic_cc_dds_lrt.rds"))

dds_vst <- DESeq2::vst(dds_lrt, blind = FALSE)

anno <- select(yeast_orgdb,
               keys=rownames(dds_lrt),
               columns=c("ENSEMBL", "COMMON")) %>%
  as_tibble() %>%
  dplyr::rename(gene_id = ENSEMBL) %>%
  dplyr::rename(gene_name = "COMMON") %>%
  dplyr::select(-SGD, -ORF) %>%
  as.data.frame() %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  filter(!is.na(gene_id))


rownames(anno) = anno$gene_id

# pca2go_at <- pca2go(dds_vst,
#                     pca_ngenes = 1000,
#                     annotation = anno,
#                     annopkg = "org.Sc.sgd.db",
#                     ensToGeneSymbol = FALSE)


pcaExplorer(dds = dds_lrt, dst = dds_vst, annotation = anno)

BiocManager::install("org.At.tair.db")
library("org.At.tair.db")

anno_at <- get_annotation_orgdb(dds_at,"org.At.tair.db", idtype = "TAIR")
