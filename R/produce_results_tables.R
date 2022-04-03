library(DESeq2)
library(tidyverse)
library(here)
library(BiocParallel)

# register(BiocParallel::MulticoreParam(10))

# note, this simply re-does what already exists in this DDS object for the sake
# of completeness. If the statistical test had not been run, this would work to
# create the same data set.
# dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

# move this over into the data processing scripts
# set aa base level to +lys
# dds$aminoAcid = factor(as.character(dds$aminoAcid),
#                        levels = c("LysHisMetLeuUra", "HisMetLeuUra"),
#                        labels = c('plus_lys', 'minus_lys'))

# dds$sample = paste(dds$genotype, dds$aminoAcid, sep = "_")
# colnames(dds) = dds$sample

# run deseq
<<<<<<< HEAD
# design(dds) = ~ genotype + aminoAcid + genotype:aminoAcid
#
# dds = DESeq(
#   dds,
#   test = "LRT",
#   reduced = ~ genotype + aminoAcid,
#   parallel = TRUE
# )
=======
design(dds) = ~ genotype + aminoAcid + genotype:aminoAcid

dds = DESeq(
  dds,
  test = "LRT",
  reduced = ~ genotype + aminoAcid
)
>>>>>>> 8fa902069a054d712991ad24f35a0a2007b11ed5

# write_rds(dds, here("data/mimic_cc_dds_lrt.rds"))

# shrink results ---------------------------------------------------------------

shrinkList = function(res_list){
  map(
    res_list,
    ~lfcShrink(dds,
               res = .,
               type = "ashr", parallel = TRUE,
               BPPARAM = BiocParallel::MulticoreParam(4)))
}

createShrunkenResList = function(dds){
  res_lists = list()

  # effect of genotype KO in +lys------------------------------

  plus_lys_contrast_list = list()

  for(i in seq(2,6)){
    contrast_vector = rep(0, length(resultsNames(dds)))
    contrast_vector[i] = 1
    plus_lys_contrast_list[[i-1]] = contrast_vector
  }

  names(plus_lys_contrast_list) = str_remove_all(resultsNames(dds)[2:6],  "genotype_|_vs_WT")

  res_lists$plus_lys =
    map(plus_lys_contrast_list,
        ~results(dds,
                 contrast = .,
                 test = "Wald"))

  # effect of minus lys AND genotype ---------------------------------------------
  # minus lys + genome - minus lys

  contrast_list = list()

  for(i in seq(2,6)){
    contrast_vector = rep(0, length(resultsNames(dds)))
    contrast_vector[i] = 1
    contrast_vector[i+6] = 1
    contrast_list[[i-1]] = contrast_vector
  }

  names(contrast_list) = str_remove_all(resultsNames(dds)[2:6],  "genotype_|_vs_WT")

  # effect of deletion in minus lys is (beta_1 + beta_3) - beta_2
  res_lists$minus_lys =
    map(contrast_list,
        ~results(dds,
                 contrast = .,
                 test = "Wald")
    )

  res_lists$minus_lys$wt = results(dds,
                                   contrast = list('Intercept', "aminoAcid_minus_lys_vs_plus_lys"),
                                   test = "Wald")

  # interaction term comparison results ------------------------------------------

  # the interaction term for condition effect in genotype III vs genotype II.
  # this tests if the condition effect is different in III compared to II

  genotype_interaction_coef_list = resultsNames(dds)[8:12]

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
    mutate(Var1 = str_remove_all(Var1, "genotype|\\.aminoAcidminus_lys"),
           Var2 = str_remove_all(Var2, "genotype|\\\\.aminoAcidminus_lys")) %>%
    mutate(tmp = paste(Var1, Var2 , sep = "_vs_")) %>%
    pull(tmp)

  res_lists$compare_interaction = compare_interaction_res_list

  res_lists$compare_interaction$eds1_rgt1_double_vs_single_effects =
    results(
      dds,
      contrast = c(0,-1,0,1,0,-1,0,-1,0,1,0,-1),
      test = "Wald"
    )

  shrunken_res_lists = map(res_lists, shrinkList)

  shrunken_res_lists
}



# write_rds(shrunken_res_lists, here("data/shrunken_res_lists.rds"))

# Genetic Interaction Tables ---------------------------------------------------

setInsigFcToZero = function(res, set_name){

  res_df = as_tibble(res, rownames = 'genes') %>%
    filter(!is.na(log2FoldChange)) %>%
    mutate(log2FoldChange = ifelse(padj > .05, 0, log2FoldChange)) %>%
    mutate(lfc_dir = ifelse(log2FoldChange > 0, "+", "-")) %>%
    mutate(lfc_dir = ifelse(log2FoldChange == 0, "0", lfc_dir)) %>%
    dplyr::rename(!!set_name := lfc_dir)
}

extractGeneticInteractionTables = function(dds){

  ## rgt1 and lys14 vs eds1
  contrast_list = list(
    plus_lys = list(
      eds1_vs_wt         = c('genotype', 'EDS1', 'WT'),
      rgt1_vs_eds1       = list('genotype_RGT1_vs_WT', 'genotype_EDS1_vs_WT'),
      lys14_vs_eds1      = list('genotype_LYS14_vs_WT', 'genotype_EDS1_vs_WT'),
      eds1rgt1_vs_rgt1   = list("genotype_EDS1_RGT1_vs_WT" , 'genotype_RGT1_vs_WT'),
      eds1lys14_vs_lys14 = list("genotype_EDS1_LYS14_vs_WT" , 'genotype_LYS14_vs_WT')
    )
  )


  res_lists = list()
  res_lists$plus_lys =
    map(contrast_list$plus_lys,
        ~results(dds,
                 contrast = .,
                 test = "Wald")
    )

  res_lists = map(res_lists, shrinkList)

  shrunken_res_df_lists = map(names(res_lists$plus_lys),
                                       ~setInsigFcToZero(res_lists$plus_lys[[.]],
                                                         .))
  names(shrunken_res_df_lists) = names(contrast_list$plus_lys)

  pattern_df = reduce(map(names(shrunken_res_df_lists),
                          ~dplyr::select(shrunken_res_df_lists[[.]], genes, .)),
                      left_join)

  pattern_df = pattern_df[complete.cases(pattern_df),]
  pattern_df = pattern_df  %>%
  unite('pattern_rgt1', c("eds1_vs_wt",
                          "rgt1_vs_eds1",
                          "eds1rgt1_vs_rgt1"),
        sep="_", remove = FALSE) %>%
  unite('pattern_lys14', c("eds1_vs_wt",
                             "lys14_vs_eds1",
                             "eds1lys14_vs_lys14"),
          sep="_", remove = FALSE)

  list(
    shrunk_res_list = res_lists,
    shrunken_res_df = shrunken_res_df_lists,
    interaction_pattern = pattern_df
  )


}
