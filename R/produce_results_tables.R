library(DESeq2)
library(tidyverse)
library(here)

# note, this simply re-does what already exists in this DDS object for the sake
# of completeness. If the statistical test had not been run, this would work to
# create the same data set.
dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

# set aa base level to +lys
dds$aminoAcid = relevel(dds$aminoAcid, "LysHisMetLeuUra")

# run deseq
design(dds) = ~ genotype + aminoAcid + genotype:aminoAcid

dds = DESeq(
  dds,
  test = "LRT",
  reduced = ~ genotype + aminoAcid
)

write_rds(dds, here("data/mimic_cc_dds_lrt.rds"))

res_lists = list()

# effect of genotype in plus lys (base) condition ------------------------------

genotype_effect_in_plus_lys_list = resultsNames(dds)[2:6]
names(genotype_effect_in_plus_lys_list) =
  str_remove_all(genotype_effect_in_plus_lys_list,
                 "genotype_|_vs_WT")


res_lists$plus_lys =
  map(genotype_effect_in_plus_lys_list,
      ~results(dds,
               name = .,
               test = "Wald"))

# effect of minus lys AND genotype ---------------------------------------------
# this is the main effect (minus aminoAcid) *plus* the interaction term

genotype_interaction_coef_list = resultsNames(dds)[8:12]
names(genotype_interaction_coef_list) =
  str_remove_all(genotype_interaction_coef_list,
                 "genotype|\\.aminoAcidHisMetLeuUra")

res_lists$minus_lys =
  map(genotype_interaction_coef_list,
      ~results(dds,
               contrast = list("aminoAcid_HisMetLeuUra_vs_LysHisMetLeuUra", .),
               test = "Wald")
  )

getInteractionCompareResults = function(row){
  results(dds,
          contrast = list(row[['Var1']], row[['Var2']]),
          test = "Wald")
}

# interaction term comparison results ------------------------------------------

# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II

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
