library(brentlabRnaSeqTools)
library(DESeq2)
library(tidyverse)
library(here)

source("R/utils.R")

dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

lys_only_dds = dds[,str_detect(dds$aminoAcid, "Lys", negate = TRUE)]

design(lys_only_dds) = formula(~genotype)

colData(lys_only_dds) = colData(lys_only_dds) %>%
  as_tibble() %>%
  droplevels() %>%
  DataFrame()

lys_only_dds = DESeq(lys_only_dds)

res_lists =
  map(
    resultsNames(lys_only_dds),
    ~lfcShrink(lys_only_dds,
               coef = .,
               type = "ashr",
               parallel = TRUE,
               BPPARAM = BiocParallel::MulticoreParam(4)))

names(res_lists) = resultsNames(lys_only_dds)

outdir = "plots/minus_lys_sanityCheck"
dir.create(outdir)
setwd(outdir)

map(names(res_lists), ~mapKEGGpathway(., res_lists[[.]],pathways$lysine_pathway))

setwd(here())

# full_res_list_shrunken = readRDS(here("data/shrunken_res_lists.rds"))

# full_res_list_shrunken$plus_lys$EDS1 %>%
#   as_tibble(rownames = 'gene') %>%
#   filter(abs(log2FoldChange) >.5, padj < .05, gene %in%
#            c("YBR115C", 'YDL131W', 'YDL182W', 'YDR234W', 'YGL202W', 'YIL094C',
#              'YIR034C', 'YJL200C', 'YNR050C'))
