library(DESeq2)
library(tidyverse)
library(topGO)
library(org.Sc.sgd.db)
library(here)

source(here("R/create_long_results_df.R"))
long_res_df = long_res_df %>%
  dplyr::rename(gene = id) %>%
  mutate(genotype = toupper(str_remove(genotype, "_$")))

dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

dds_vst = varianceStabilizingTransformation(dds, blind=FALSE)

colnames(dds_vst) = paste0('sample', dds_vst$id)
dds_vst$sample = paste0('sample', dds_vst$id)

vst_long = assays(dds_vst)[[1]] %>%
  as_tibble(rownames = 'gene') %>%
  pivot_longer(-gene, names_to = 'sample',
               values_to = 'vst') %>%
  left_join(as_tibble(colData(dds_vst)))

# code from james.mcdonald https://support.bioconductor.org/p/9142025/#9142035
gn2go = mapIds(org.Sc.sgd.db, keys(org.Sc.sgd.db), "GOALL", "ORF", multiVals = "list")
# filter out NA mappings
gn2go = gn2go[!sapply(gn2go, function(x) all(is.na(x)))]

go2gn = inverseList(gn2go)

topgo_res = list(
  res_lists = readRDS(here("data/topGO_out/go_results_objects.rds")),
  venn_sets = readRDS(here("data/topGO_out/go_results_objects_by_venn_set.rds"))
)

go_term_plot = function(go_list, go2gn, long_res_df){

  tibble(
    gene = unlist(go2gn[go_list]),
    go_term = unlist(map(go_list, ~rep(., length(go2gn[[.]]))))) %>%
    left_join(long_res_df) %>%
    filter(genotype %in% c("EDS1", "RGT1", "EDS1_RGT1")) %>%
    group_by(genotype) %>%
    mutate(genotype = factor(genotype, levels = c('EDS1', 'RGT1', 'EDS1_RGT1'))) %>%
    ggplot(aes(genotype, shrunk_lfc, group = gene)) +
    geom_point() +
    geom_line() +
    facet_grid(rows = vars(condition), cols = vars(go_term))
}

x = topgo_res$venn_sets$bp$plus_lys$eds1_rgt1_intersect_rgt1Effect$results$summary$GO.ID
plt = go_term_plot(x, go2gn, long_res_df)
