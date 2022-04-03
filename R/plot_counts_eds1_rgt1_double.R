library(tidyverse)
library(patchwork)
library(DESeq2)

dds_lrt = readRDS(here("data/mimic_cc_dds_lrt.rds"))

interaction_term_res = results(dds)

plotEds1Rgt1DoubleCounts = function(gene_name, gene_id){

gene_index =which(rownames(interaction_term_res) == gene_id) #which.min(interaction_term_res$padj) #which(rownames(interaction_term_res) == "YAL038W") #

gene_plt <- plotCounts(dds, gene_index,
                       intgroup = c("genotype","aminoAcid"),
                       returnData = TRUE,
                       normalized = TRUE) %>%
  mutate(count = log2(count))

ylim_min = 0
ylim_max = 15

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

plt4 = gene_plt %>%
  filter(genotype %in% c("WT", "EDS1_RGT1")) %>%
  ggplot(
    aes(x = genotype, y = count, color = aminoAcid, group = aminoAcid)) +
  geom_point() +
  stat_summary(fun=mean, geom="line") +
  ylim(ylim_min,ylim_max)+
  labs(y = "", x = "")

ptchwrk = plt1 + plt2+ plt4 + plot_layout(guides = 'collect')

ptchwrk = ptchwrk + plot_annotation(
  title = gene_name,
  subtitle = gene_id
) & theme_minimal()

# ggsave(here("plots/condition_specific_THI4.png"), ptchwrk, device = "pdf")

ptchwrk

}
