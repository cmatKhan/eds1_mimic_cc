library(GeneTonic)
# library(tidyverse)
# library(pcaExplorer)
# library(DESeq2)
# library(here)
# library(org.Sc.sgd.db)
#
# res_list = readRDS(here("data/shrunken_res_lists.rds"))
#
# yeast_orgdb = org.Sc.sgd.db
#
# x = get(ID, org.Sc.sgdGO)
#
# select(yeast_orgdb, keys=ID,columns=c("ENTREZID", "GO"))
#
# get_go_sets = function(org_db, gene_set){
#   gene_to_go = select(yeast_orgdb,
#                       keys=gene_set,
#                       columns=c("ENTREZID", "GO"))
#
#   unmapped = gene_to_go %>%
#     filter(is.na(ONTOLOGY))
#
#   gene_to_go = gene_to_go %>%
#     filter(!is.na(ONTOLOGY)) %>%
#     group_by(ONTOLOGY) %>%
#     group_split()
#
#   reshapeGoToList = function(go_table, go_category){
#     x = list()
#     x[[go_category]] =
#       split(go_table[,1],go_table[,4]) %>%
#       map(pull, "ORF")
#     x
#   }
#
#   names(gene_to_go) = unique(bind_rows(gene_to_go)$ONTOLOGY)
#
#   out_list = unlist(map(names(gene_to_go),
#              ~reshapeGoToList(gene_to_go[[.]], .)),
#          recursive = FALSE)
#   out_list$unmapped = pull(unmapped, 'ORF')
#   out_list
# }
#
# yeast_gs = get_go_sets(yeast_orgdb, rownames(dds_lrt))

library("macrophage")

data("gse", package = "macrophage")

dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
# changing the ids to Ensembl instead of the Gencode used in the object
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)


data("res_de_macrophage")
head(res_macrophage_IFNg_vs_naive)

library("AnnotationDbi")
de_symbols_IFNg_vs_naive <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds_macrophage)$SYMBOL[rowSums(counts(dds_macrophage)) > 0]

library("topGO")
topgoDE_macrophage_IFNg_vs_naive <-
  pcaExplorer::topGOtable(de_symbols_IFNg_vs_naive,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "symbol",
                          topTablerows = 500)

res_enrich_macrophage <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
colnames(res_enrich_macrophage)

library("org.Hs.eg.db")
anno_df <- data.frame(
  gene_id = rownames(dds_macrophage),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds_macrophage), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds_macrophage)
)

res_enrich_macrophage <- get_aggrscores(res_enrich = res_enrich_macrophage,
                                        res_de = res_macrophage_IFNg_vs_naive,
                                        annotation_obj = anno_df,
                                        aggrfun = mean)

GeneTonic(dds = dds_macrophage,
          res_de = res_macrophage_IFNg_vs_naive,
          res_enrich = res_enrich_macrophage,
          annotation_obj = anno_df,
          project_id = "GT1")
