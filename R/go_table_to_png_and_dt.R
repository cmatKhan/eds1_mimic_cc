library(gridExtra)
library(grid)
library(here)
library(tidyverse)

goTableToPngAndDT = function(dir_prefix, go_cat, cond, tablename, df, height=300, width=1500){

  # create output path
  out_dirpath = file.path(dir_prefix, go_cat, cond)
  dir.create(out_dirpath, recursive = TRUE)
  outpath = file.path(out_dirpath, paste0(tablename, ".png"))

  # print data table to file
  png(outpath, height = height, width = width, bg = 'white')
    df %>%
    tableGrob(rows = NULL) -> g
    grid.newpage()
    ttheme_default(base_size = 30)
    grid.draw(g)
  dev.off()
}

byGeno = function(out_dir_prefix, go_cat, cond, geno_list){

  for(i in names(geno_list)){
    goTableToPngAndDT(out_dir_prefix,
                      go_cat,
                      cond,
                      i,
                      geno_list[[i]]$results$summary)
  }
  # map(names(geno_list)
  #     ~goTableToPngAndDT(out_dir_prefix,
  #                        go_cat,
  #                        cond,
  #                        .,
  #                        geno_list[[.]]$results$summary))
}

byCond = function(out_dir_prefix, go_cat, cond_list){

  x = 0
  for(i in names(cond_list)){
    byGeno(out_dir_prefix,
           go_cat,
           i,
           cond_list[[i]])
  }
  # map(names(cond_list)
  #     ~byGeno(out_dir_prefix,
  #             go_cat,
  #             .,
  #             cond_list[[.]]))

}

topgo_res = list(
  res_lists = readRDS(here("data/topGO_out/go_results_objects.rds")),
  venn_sets = readRDS(here("data/topGO_out/go_results_objects_by_venn_set.rds"))
)

venn_set_table_dir = here("plots/venn_set/")

map(names(topgo_res$venn_sets), ~byCond(venn_set_table_dir, ., topgo_res$venn_sets[[.]]))
