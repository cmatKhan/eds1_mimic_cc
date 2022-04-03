library("GeneTonic")
library("DESeq2")
library("topGO")
library(tidyverse)
library(org.Sc.sgd.db)
library(here)

res_list = readRDS(here("data/shrunken_res_lists.rds"))
dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))


res_list$plus_lys$EDS1$id <- rowData(dds)$id

anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Sc.sgd.db, keys = rownames(dds), column = 'ORF', keytype = 'ORF'),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)


library("AnnotationDbi")
de_ids <- deseqresult2df(res_list$plus_lys$EDS1, FDR = 0.05)$id
de_ids = de_ids[de_ids %in% anno_df$gene_name[!is.na(anno_df$gene_name)]]
bg_ids <- rownames(dds)[rowSums(counts(dds)) > 0 & rownames(dds) %in% anno_df$gene_name[!is.na(anno_df$gene_name)]]

topgoDE <-
  pcaExplorer::topGOtable(de_ids,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Sc.sgd.db",
                          geneID = "ENSEMBL",
                          topTablerows = 500)

data("res_enrich_macrophage")
head(topgoDE_macrophage_IFNg_vs_naive, 2)

res_enrich_macrophage <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
colnames(res_enrich_macrophage)

gene_list = res_list$plus_lys$EDS1 %>%
  as_tibble(rownames = 'gene') %>%
  filter(abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  pull(gene)

ggo <- groupGO(gene     = gene_list,
               OrgDb    = org.Sc.sgd.db,
               ont      = "BP",
               keyType = "ORF",
               level    = 3,
               readable = TRUE)

head(ggo)

# new( "topGOdata", ontology=onts[i], allGenes = algSxT, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID = "SYMBOL" )

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = gene_list, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, mapping="org.Sc.sgd.db", ID = "ENSEMBL")
