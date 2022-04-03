library(DESeq2)
library(tidyverse)
library(patchwork)

# see https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


plotGeneCountCustom = function(dds, gene_id, gene_symbol, ylim_min = 0, ylim_max = 15){

  interaction_term_res = results(dds)

  gene_index = which(rownames(interaction_term_res) == gene_id) #

  gene_plt <- plotCounts(dds, gene_index,
                         intgroup = c("genotype","aminoAcid"),
                         returnData = TRUE,
                         normalized = TRUE) %>%
    mutate(count = log2(count))

  plt1 = gene_plt %>%
    filter(genotype %in% c("WT", "EDS1")) %>%
    ggplot(
      aes(x = genotype, y = count, color = aminoAcid, group = aminoAcid)) +
    geom_point() +
    stat_summary(fun=mean, geom="line") +
    ylim(ylim_min,ylim_max)+
    labs(y = "log2(norm_count)", x = "")

  plt2 = gene_plt %>%
    filter(genotype %in% c("WT", "RGT1")) %>%
    ggplot(
      aes(x = genotype, y = count, color = aminoAcid, group = aminoAcid)) +
    geom_point() +
    stat_summary(fun=mean, geom="line") +
    ylim(ylim_min,ylim_max)+
    labs(y = "", x = "genotype")

  plt3 = gene_plt %>%
    filter(genotype %in% c("WT", "LYS14")) %>%
    ggplot(
      aes(x = genotype, y = count, color = aminoAcid, group = aminoAcid)) +
    geom_point() +
    stat_summary(fun=mean, geom="line") +
    ylim(ylim_min,ylim_max)+
    labs(y = "", x = "genotype")

  plt4 = gene_plt %>%
    filter(genotype %in% c("WT", "EDS1_RGT1")) %>%
    ggplot(
      aes(x = genotype, y = count, color = aminoAcid, group = aminoAcid)) +
    geom_point() +
    stat_summary(fun=mean, geom="line") +
    ylim(ylim_min,ylim_max)+
    labs(y = "", x = "")

  plt5 = gene_plt %>%
    filter(genotype %in% c("WT", "EDS1_LYS14")) %>%
    ggplot(
      aes(x = genotype, y = count, color = aminoAcid, group = aminoAcid)) +
    geom_point() +
    stat_summary(fun=mean, geom="line") +
    ylim(ylim_min,ylim_max)+
    labs(y = "", x = "")

  ptchwrk = plt1 + plt2+ plt3 +plt4 + plt5+ plot_layout(guides = 'collect')

  ptchwrk = ptchwrk + plot_annotation(
    title = gene_symbol,
    subtitle = gene_id
  ) & theme_minimal()

  # ggsave(here("plots/condition_specific_THI4.png"), ptchwrk, device = "pdf")

  list(
    plts = list(
      eds1 = plt1,
      rgt1 = plt2,
      lys14 = plt3,
      eds1_rgt1 = plt4,
      eds1_lys14 = plt5
    ),
    combined = ptchwrk
  )
}
