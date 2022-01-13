library(STRINGdb)
library(DESeq2)
library(tidyverse)
library(here)

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

# see https://stringdb-static.org/download/species.v11.5.txt
scer_tax_id = 4932

string_db = STRINGdb$new(version = "11",
                         species = scer_tax_id,
                         score_threshold = 200, # default is 400. 200 was in vignette
                         input_directory = here("data/STRINGdb")) # do not specify this and data will be tmp/cached

string_db

#'
#' @param res_table_name name of the results table
#' @param res_direction results are parsed into 'up' and 'down' corresponding to
#'   direction of lfc
#' @param string_db_mapped the string_db_mapped object
#'
printNetworks = function(res_table_name, res_direction, string_db_mapped, dir_outpath){

  message(paste0("working on direction: ", res_direction))

  plot_dirpath = file.path(dir_outpath, res_table_name, "stringdb_plots", res_direction )

  dir.create(plot_dirpath, recursive = TRUE)

  payload_id = string_db$post_payload(string_db_mapped$STRING_id,
                                      colors = string_db_mapped$color )

  clustersList = string_db$get_clusters(string_db_mapped$STRING_id)

  for(i in seq(length(clustersList))){
    if(length(clustersList[[i]]) > 1){

      outpath = file.path(plot_dirpath, paste0(res_table_name, "_", i, ".pdf"))

      pdf(outpath)
      string_db$plot_network(clustersList[[i]], payload_id = payload_id)
      dev.off()
    }
  }

}

createNetworkMaps = function(res_table_name,
                             res,
                             dir_outpath = "",
                             padj_thres = .05,
                             log2fc_thres = 2){

  message(paste0("working on table: ", res_table_name))

  res = res %>%
    as_tibble(rownames = 'gene') %>%
    filter(padj < padj_thres, abs(log2FoldChange) > log2fc_thres) %>%
    select(gene, log2FoldChange, padj) %>%
    as.data.frame()

  if(nrow(res) > 2){
    res_mapped = string_db$map(
      res,
      "gene",
      removeUnmappedRows = TRUE )

    # full = string_db$map(res, "gene", removeUnmappedRows = TRUE )

    mapped = list(

      up = string_db$add_diff_exp_color(
        subset(res_mapped,
               log2FoldChange > 0),
        logFcColStr="log2FoldChange"),

      down = string_db$add_diff_exp_color(
        subset(res_mapped,
               log2FoldChange < 0),
        logFcColStr="log2FoldChange")
    )

    map(names(mapped), ~printNetworks(res_table_name,
                                      .,
                                      mapped[[.]],
                                      dir_outpath))
  } else{
    message("less than 2 genes at these thresholds: ", res_table_name)
  }

}

mapResultsTables = function(res_list_name, res_list){

  message(paste0("working on res_list: ", res_list_name))

  # if(res_list_name == "compare_interaction"){
  #   message("debug here")
  # }

  dir_outpath = here("data", res_list_name)
  dir.create(dir_outpath)

  map(names(res_list),
      ~createNetworkMaps(.,
                         res_list[[.]],
                         dir_outpath))
}

map(names(shrunken_res_lists), ~mapResultsTables(., shrunken_res_lists[[.]]))
