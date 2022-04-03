library(tidyverse)
library(DESeq2)
library(tidyverse)
library(topGO)
library(org.Sc.sgd.db)
library(here)
library(foreach)
library(doParallel)

source(here("R/venn_diagram_gene_lists.R"))

cc_df = read_csv("../shared_data/Combined_Eds1_Rgt1_Lys14_CCTargets.csv")
dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))
shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

venn_lists = list(
  sets = list(
    plus_lys = list(
    EDS1 = pull(filter(cc_df, eds1_plus_lys_pval < .05), id),
    LYS14 = pull(filter(cc_df, lys14_plus_lys_pval < .05), id),
    RGT1 = pull(filter(cc_df, rgt1_plus_lys_pval < .05), id)
  ),
  minus_lys = list(
    EDS1 = pull(filter(cc_df, eds1_minus_lys_pval < .05), id),
    LYS14 = pull(filter(cc_df, Lys14_minus_lys_pval < .05), id),
    RGT1 = pull(filter(cc_df, rgt1_minus_lys_pval < .05), id)
  )))

  venn_lists$plots = list(
  plots = list(
    plus_lys = ggVennDiagram(venn_lists$sets$plus_lys,
                  label = "count",
                  label_alpha = 0,
                  edge_size = 0,
                  label_geom = c("label", "text"),
                  label_size = 9) +
  scale_fill_distiller(palette = 'PuBuGn', direction = 1) +
  ggtitle("Plus Lysine",
          subtitle = paste0("padj < .05")) +
  theme(legend.position="none",
        plot.title = element_text(size = 25)),

  minus_lys = ggVennDiagram(venn_lists$sets$minus_lys,
                            label = "count",
                            label_alpha = 0,
                            edge_size = 0,
                            label_geom = c("label", "text"),
                            label_size = 9) +
  scale_fill_distiller(palette = 'YlOrBr', direction = 1) +
  ggtitle("Minus Lysine",
          subtitle = paste0("padj < .05")) +
  theme(legend.position="none",
        plot.title = element_text(size = 25)),

 eds1_plus_minus_lys =
   ggVennDiagram(list(plus_lys = venn_lists$sets$plus_lys$EDS1,
                      minus_lys = venn_lists$sets$minus_lys$EDS1),
              label = "count",
              label_alpha = 0,
              edge_size = 0,
              label_geom = c("label", "text"),
              label_size = 9) +
  scale_fill_distiller(palette = 'YlGn', direction = -1) +
  ggtitle("EDS1 in plus and minus lysine",
          subtitle = paste0("padj < .05")) +
  theme(legend.position="none",
        plot.title = element_text(size = 25)),

 rgt1_plus_minus_lys =
   ggVennDiagram(list(plus_lys = venn_lists$sets$plus_lys$RGT1,
                      minus_lys = venn_lists$sets$minus_lys$RGT1),
                 label = "count",
                 label_alpha = 0,
                 edge_size = 0,
                 label_geom = c("label", "text"),
                 label_size = 9) +
   scale_fill_distiller(palette = 'YlGn', direction = -1) +
   ggtitle("RGT1 in plus and minus lysine",
           subtitle = paste0("padj < .05")) +
   theme(legend.position="none",
         plot.title = element_text(size = 25)),

 lys14_plus_minus_lys =
   ggVennDiagram(list(plus_lys = venn_lists$sets$plus_lys$LYS14,
                      minus_lys = venn_lists$sets$minus_lys$LYS14),
                 label = "count",
                 label_alpha = 0,
                 edge_size = 0,
                 label_geom = c("label", "text"),
                 label_size = 9) +
   scale_fill_distiller(palette = 'YlGn', direction = -1) +
   ggtitle("LYS14 in plus and minus lysine",
           subtitle = paste0("padj < .05")) +
   theme(legend.position="none",
         plot.title = element_text(size = 25))
  ))

venn_for_topgo = getVennResults(shrunken_res_lists, abs_lfc_thres = 1, padj_thres = .05)

up_plus_lys_eds1_rgt1 = list(
  venn_for_topgo$genes$plus_lys$up[['eds1_vs_rgt1']]
)

up_plus_lys_eds1_rgt1 = up_plus_lys_eds1_rgt1[[1]]
names(up_plus_lys_eds1_rgt1) = paste0(names(up_plus_lys_eds1_rgt1), "_de")

up_plus_lys_eds1_rgt1$eds1_cc = pull(filter(cc_df, eds1_plus_lys_pval < .05), id)
up_plus_lys_eds1_rgt1$rgt1_cc = pull(filter(cc_df, rgt1_plus_lys_pval < .05), id)

df = tibble(
  gene = rownames(dds),
  eds1_de = ifelse(rownames(dds) %in% venn_for_topgo$genes$plus_lys$up$singles$EDS1, TRUE, FALSE),
  rgt1_de = ifelse(rownames(dds) %in% venn_for_topgo$genes$plus_lys$up$singles$RGT1, TRUE, FALSE),
  eds1_rgt1_de = ifelse(rownames(dds) %in% venn_for_topgo$genes$plus_lys$up$eds1_vs_rgt1$EDS1_RGT1, TRUE, FALSE),
  eds1_cc = ifelse(rownames(dds) %in% pull(filter(cc_df, eds1_plus_lys_pval < .05), id), TRUE, FALSE),
  rgt1_cc = ifelse(rownames(dds) %in% pull(filter(cc_df, rgt1_plus_lys_pval < .05), id), TRUE, FALSE)
) %>%
  pivot_longer(-gene, names_to = 'assay', values_to = 'membership')

# UpSetR::upset(UpSetR::fromList(up_plus_lys_eds1_rgt1),
#       order.by = "degree",
#       intersections = list(list('rgt1_cc', 'RGT1_de'),
#                            list('eds1_cc', 'EDS1_de'),
#                            list('rgt1_cc', 'eds1_cc', 'EDS1_RGT1_de')))
#
# ComplexUpset::upset(df, 'membership', name = 'assay')

gff = read_tsv("~/ref/S288C_R64/GCF_000146045.2_R64_genomic.gff", col_names = FALSE, skip = 8) %>%
  dplyr::select(X1, X9) %>%
  mutate(id =
           str_extract(X9, "(?=locus_tag=)()(.*?)(?=;)")) %>%
  mutate(id = str_remove(id, "locus_tag=")) %>%
  dplyr::rename(chr = X1) %>%
  dplyr::select(-X9)

x = cc_df %>%
  pivot_longer(-c(id, symbol),
               names_to = 'cond',
               values_to = "pval") %>%
  left_join(gff)




