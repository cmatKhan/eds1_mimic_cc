library(pathview)
library(DESeq2)
library(tidyverse)
library(here)

source(here("R/utils.R"))

# read in results
shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

# note: not changed from testing -- the results in the actual plots results are
# still from the presentation 2021 01 17
outdir = "plots/minus_lys/tca_pathway"
dir.create(outdir)
setwd(outdir)

# THERE IS AN ERROR IN THE MAPPING FROM ID TO SYMBOL IN AT LEAST THE TCA CYCLE

# note: this was changed since presentation to the utils function. thershold was also
# accidently set to .05. changed to .5. doesn't affect results * much * -- removed
# one gene which had lfc of ~-.4
x = map(names(shrunken_res_lists$minus_lys),
  ~mapKEGGpathway(., shrunken_res_lists$minus_lys[[.]], pathway_id = pathways$tca_cycle, lfc_thres = 1))

setwd(here())
