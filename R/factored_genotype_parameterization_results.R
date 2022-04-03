library(DESeq2)
library(tidyverse)
library(here)

# note, this simply re-does what already exists in this DDS object for the sake
# of completeness. If the statistical test had not been run, this would work to
# create the same data set.
dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

dds$aminoAcid = factor(dds$aminoAcid,
                       levels = c("LysHisMetLeuUra", "HisMetLeuUra"),
                       labels = c("plusLys", "minusLys"))

colData(dds) = colData(dds) %>%
  as_tibble() %>%
  mutate(EDS1 = ifelse(genotype == "EDS1", 'del', NA)) %>%
  mutate(EDS1 = ifelse(is.na(EDS1) & genotype == "WT", "WT", EDS1)) %>%
  mutate(EDS1 = ifelse(is.na(EDS1) &
                         genotype == "EDS1_RGT1", "del", EDS1)) %>%
  mutate(EDS1 = ifelse(is.na(EDS1) &
                         genotype == "EDS1_LYS14", "del", EDS1)) %>%
  mutate(EDS1 = ifelse(is.na(EDS1) &
                         genotype == "LYS14", "WT", EDS1)) %>%
  mutate(EDS1 = ifelse(is.na(EDS1) &
                         genotype == "RGT1", "WT", EDS1)) %>%
  mutate(EDS1 = factor(EDS1, levels = c("WT", 'del'))) %>%
  # RGT1
  mutate(RGT1 = ifelse(genotype == "RGT1", "del", NA)) %>%
  mutate(RGT1 = ifelse(is.na(RGT1) &
                         genotype == "WT", "WT", RGT1)) %>%
  mutate(RGT1 = ifelse(is.na(RGT1) &
                         genotype == "EDS1_RGT1", "del", RGT1)) %>%
  mutate(RGT1 = ifelse(is.na(RGT1) &
                         genotype == "EDS1_LYS14", "WT", RGT1)) %>%
  mutate(RGT1 = ifelse(is.na(RGT1) &
                         genotype == "LYS14", "WT", RGT1)) %>%
  mutate(RGT1 = ifelse(is.na(RGT1) &
                         genotype == "EDS1", "WT", RGT1)) %>%
  mutate(RGT1 = factor(RGT1, levels = c("WT", 'del'))) %>%
  # LYS14
  mutate(LYS14 = ifelse(genotype == "LYS14", "del", NA)) %>%
  mutate(LYS14 = ifelse(is.na(LYS14) &
                         genotype == "WT", "WT", LYS14)) %>%
  mutate(LYS14 = ifelse(is.na(LYS14) &
                         genotype == "EDS1_RGT1", "WT", LYS14)) %>%
  mutate(LYS14 = ifelse(is.na(LYS14) &
                         genotype == "EDS1_LYS14", "del", LYS14)) %>%
  mutate(LYS14 = ifelse(is.na(LYS14) &
                         genotype == "RGT1", "WT", LYS14)) %>%
  mutate(LYS14 = ifelse(is.na(LYS14) &
                         genotype == "EDS1", "WT", LYS14)) %>%
  mutate(LYS14 = factor(LYS14, levels = c("WT", 'del'))) %>%
  # double
  mutate(double = ifelse(genotype == "EDS1_RGT1", "EDS1_RGT1", NA)) %>%
  mutate(double = ifelse(is.na(EDS1_RGT1) &
                              genotype == "WT", "WT", EDS1_RGT1)) %>%
  mutate(double = ifelse(is.na(EDS1_RGT1) &
                              genotype == "EDS1_LYS14", "EDS1_LYS14", EDS1_RGT1)) %>%
  mutate(double = ifelse(is.na(EDS1_RGT1) &
                              genotype == "LYS14", "WT", EDS1_RGT1)) %>%
  mutate(double = ifelse(is.na(EDS1_RGT1) &
                              genotype == "EDS1", "del", EDS1_RGT1)) %>%
  mutate(double = ifelse(is.na(EDS1_RGT1) &
                              genotype == "RGT1", "del", EDS1_RGT1)) %>%
  mutate(double = factor(EDS1_RGT1, levels = c("WT", 'del'))) %>%

  # EDS1_RGT1 Double
  mutate(EDS1_RGT1 = ifelse(genotype == "EDS1_RGT1", "del", NA)) %>%
  mutate(EDS1_RGT1 = ifelse(is.na(EDS1_RGT1) &
                              genotype == "WT", "WT", EDS1_RGT1)) %>%
  mutate(EDS1_RGT1 = ifelse(is.na(EDS1_RGT1) &
                              genotype == "EDS1_LYS14", "del", EDS1_RGT1)) %>%
  mutate(EDS1_RGT1 = ifelse(is.na(EDS1_RGT1) &
                              genotype == "LYS14", "WT", EDS1_RGT1)) %>%
  mutate(EDS1_RGT1 = ifelse(is.na(EDS1_RGT1) &
                              genotype == "EDS1", "del", EDS1_RGT1)) %>%
  mutate(EDS1_RGT1 = ifelse(is.na(EDS1_RGT1) &
                              genotype == "RGT1", "del", EDS1_RGT1)) %>%
  mutate(EDS1_RGT1 = factor(EDS1_RGT1, levels = c("WT", 'del'))) %>%
  # EDS1_LYS14 Double
  mutate(EDS1_LYS14 = ifelse(genotype == "EDS1_LYS14", "del", NA)) %>%
  mutate(EDS1_LYS14 = ifelse(is.na(EDS1_LYS14) &
                              genotype == "WT", "WT", EDS1_LYS14)) %>%
  mutate(EDS1_LYS14 = ifelse(is.na(EDS1_LYS14) &
                              genotype == "EDS1_RGT1", "WT", EDS1_LYS14)) %>%
  mutate(EDS1_LYS14 = ifelse(is.na(EDS1_LYS14) &
                              genotype == "LYS14", "del", EDS1_LYS14)) %>%
  mutate(EDS1_LYS14 = ifelse(is.na(EDS1_LYS14) &
                              genotype == "EDS1", "del", EDS1_LYS14)) %>%
  mutate(EDS1_LYS14 = ifelse(is.na(EDS1_LYS14) &
                              genotype == "RGT1", "WT", EDS1_LYS14)) %>%
  mutate(EDS1_LYS14 = factor(EDS1_LYS14, levels = c("WT", 'del'))) %>%
  DataFrame()


design(dds) = ~EDS1 + RGT1 + LYS14 + EDS1_RGT1 + EDS1:EDS1_RGT1
DESeq(dds)

# run deseq
design(dds) = ~ genotype + aminoAcid + genotype:aminoAcid

dds = DESeq(
  dds,
  test = "LRT",
  reduced = ~ genotype + aminoAcid
)

# write_rds(dds, here("data/mimic_cc_dds_lrt.rds"))

res_lists = list()

# effect of genotype in plus lys (base) condition ------------------------------

genotype_effect_in_plus_lys_list = resultsNames(dds)[1:6]

names(genotype_effect_in_plus_lys_list) =
  str_remove_all(genotype_effect_in_plus_lys_list,
                 "genotype_|_vs_WT")


res_lists$plus_lys =
  map(genotype_effect_in_plus_lys_list,
      ~results(dds,
               name = .,
               test = "Wald"))

# effect of minus lys AND genotype ---------------------------------------------
# minus lys + genome - minus lys

genotype_interaction_coef_list = resultsNames(dds)[8:12]
names(genotype_interaction_coef_list) =
  str_remove_all(genotype_interaction_coef_list,
                 "genotype|\\.aminoAcidHisMetLeuUra")

res_lists$minus_lys =
  map(genotype_interaction_coef_list,
      ~results(dds,
               contrast = list(., "aminoAcid_HisMetLeuUra_vs_LysHisMetLeuUra"),
               test = "Wald")
  )

res_lists$minus_lys$wt = results(dds,
                                 contrast = list("aminoAcid_HisMetLeuUra_vs_LysHisMetLeuUra", 'Intercept'),
                                 test = "Wald")

# interaction term comparison results ------------------------------------------

# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II

getInteractionCompareResults = function(row){
  results(dds,
          contrast = list(row[['Var1']], row[['Var2']]),
          test = "Wald")
}

genotype_permutations =
  expand.grid(genotype_interaction_coef_list, genotype_interaction_coef_list) %>%
  filter(Var1 != Var2)

compare_interaction_res_list =
  apply(genotype_permutations, 1, getInteractionCompareResults)

names(compare_interaction_res_list) =
  genotype_permutations %>%
  mutate(Var1 = str_remove_all(Var1, "genotype|\\.aminoAcidHisMetLeuUra"),
         Var2 = str_remove_all(Var2, "genotype|\\.aminoAcidHisMetLeuUra")) %>%
  mutate(tmp = paste(Var1, Var2 , sep = "_vs_")) %>%
  pull(tmp)

res_lists$compare_interaction = compare_interaction_res_list

res_lists$compare_interaction$eds1_rgt1_double_vs_single_effects =
  results(
    dds,
    contrast = c(0,-1,0,1,0,-1,0,-1,0,1,0,-1),
    test = "Wald"
  )

# shrink results ---------------------------------------------------------------

shrinkList = function(res_list){
  map(
    res_list,
    ~lfcShrink(dds,
               res = .,
               type = "ashr",
               parallel = TRUE,
               BPPARAM = BiocParallel::MulticoreParam(4)))
}

shrunken_res_lists = map(res_lists, shrinkList)

write_rds(shrunken_res_lists, here("data/shrunken_res_lists.rds"))
