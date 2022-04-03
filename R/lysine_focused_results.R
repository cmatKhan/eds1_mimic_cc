library(tidyverse)
library(DESeq2)
library(here)
library(org.Sc.sgd.db)

dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

# extract results --------------------------------------------------------------

getLysineFocusedResults = function(dds){

  geno_names = str_remove_all(resultsNames(dds)[8:12],
                              "genotype|\\.aminoAcidminus_lys")

  res_list = list(
    effect_of_minus_lys_on_wt = results(
      dds,
      contrast = c("aminoAcid", "minus_lys", "plus_lys"),
      test = "Wald"
    ),
    # if minus lys is more than plus lys in the KO, positive fold change
    effect_of_minus_lys_on_del = map(resultsNames(dds)[8:12],
                      ~results(dds,
                               contrast = list(., "aminoAcid_minus_lys_vs_plus_lys"),
                               test = "Wald") ),

    effect_of_del_in_minus_lys = map2(resultsNames(dds)[2:6],
                                      resultsNames(dds)[8:12],
                                      ~results(dds, contrast = list(c(.x,.y)),
                                               test = "Wald")),

    effect_of_del_in_plus_lys = map(resultsNames(dds)[2:6],
                                    ~results(dds, name = ., test = "Wald"))
  )
  names(res_list$effect_of_minus_lys_on_del) = geno_names

  names(res_list$effect_of_del_in_minus_lys) = geno_names

  names(res_list$effect_of_del_in_plus_lys) = str_remove_all(resultsNames(dds)[2:6], "genotype_|_vs_WT")


  res_list
}

shrunken_res_lists = getLysineFocusedResults(dds)

shrunken_res_lists$effect_of_minus_lys_on_wt = lfcShrink(dds,
              res = shrunken_res_lists$effect_of_minus_lys_on_wt,
              type = "ashr")

shrunken_res_lists$effect_of_minus_lys_on_del =
  map(shrunken_res_lists$effect_of_minus_lys_on_del,
      ~lfcShrink(dds, res = ., type = "ashr"))

shrunken_res_lists$effect_of_del_in_minus_lys =
  map(shrunken_res_lists$effect_of_del_in_minus_lys,
      ~lfcShrink(dds, res = ., type = "ashr"))

shrunken_res_lists$effect_of_del_in_plus_lys =
  map(shrunken_res_lists$effect_of_del_in_plus_lys,
      ~lfcShrink(dds, res = ., type = "ashr"))







