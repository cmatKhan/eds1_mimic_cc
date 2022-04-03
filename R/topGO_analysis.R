library(DESeq2)
library(tidyverse)
library(topGO)
library(org.Sc.sgd.db)
library(here)
library(foreach)
library(doParallel)

source(here("R/venn_diagram_gene_lists.R"))

# QUESTIONS
## should allGene be set to the same background set for all sets, or should it
## be the background set for a given set (meaning the geneSelectionFunction
## would select a subset of the given set to call significant as opposed to
## allGenes being all genes in the yeast genome and the selectionFunction
## subsetting a subset of that since the results are pre-filtered by padj)

## After a lot of exploration, setting the function resultsFilterGetGenes()
## such that it does not filter on log2FC or padj and just returns a named
## vector where the names are padj and values are padj, and then setting the
## selectGene function to filter based on a padj thres, returns results that
## are most interpretable (meaning terms that you would expect exist in the
## tables). This is a worthwhile observation to share.

# define functions -------------------------------------------------------------

# input results obj, output named list of genes
# filter by a threshold on padj
# output named list, values are padj, names are gene ids
resultsFilterGetGenes = function(res_obj, padj_thres = .05){

  # fltr_df = res_obj %>%
  # as_tibble(rownames = 'gene') %>%
  #   filter(padj < padj_thres)
  #
  # gene_names = fltr_df$gene
  #
  # fltr_df %>%
  #   pull(log2FoldChange) %>%
  #   stats::setNames(gene_names)

  fltr_df = res_obj %>%
    as_tibble(rownames = 'gene')

  gene_names = fltr_df$gene

  fltr_df %>%
    pull(padj) %>%
    replace_na(100) %>%
    stats::setNames(gene_names)

}

#' helper function for topGO to filter gene_lists
#' see the function below or the topGO docs
#' @param lfc this is intended to function on each element in an array, so
#'   padj is a adjusted p_value from, for example, gene_list (see topGO function
#'   below)
#' @param lfc_thres threshold by which to filter -- set default so this
#'   function works in topGO constructor
topDiffGenes = function(metric, metric_thres = 1){
  # return(abs(lfc) > lfc_thres)
  return(metric < .05)
}

#' @param gene_list a named vector where the values are padj. filter the list
#'   for abs(log2FC) above some threshold first
#' @param gn2go a gene to go mapping -- see the code above for how to create
#'   this
#' @param ontology_category one of c("BP", "CC", "MF")
#' @param res_name name of the gene_list, eg if from EDS1 ko, maybe 'eds1'
#' @param plot_outdir path to the out directory for plots
#' @param go_topNodes an integer value, default to 10. See the topGO docs on
#'   how this is used
#'
#' @return a list including the topGO obj, results tables including a summary,
#'   some diagnostics and the GO network plot
createTopGOobject = function(gene_list,
                             gn2go,
                             ontology_category,
                             res_name,
                             plot_outdir,
                             go_topNodes = 10){

  out_list = list()

  out_list$topGO_obj = new("topGOdata",
                            ontology = ontology_category,
                            allGenes = gene_list[names(gene_list) %in% names(gn2go)],
                            geneSel = topDiffGenes,
                            gene2GO = gn2go,
                            annot = annFUN.gene2GO,
                            )

  getResults = function(out_list, go_topNodes){

    out_list$results = list(
      fisher = runTest(out_list$topGO_obj,
                        algorithm = "classic",
                        statistic = "fisher"),
      ks = runTest(out_list$topGO_obj,
                    algorithm = "classic",
                    statistic = "ks"),
      ks_elim = runTest(out_list$topGO_obj,
                      algorithm = "elim",
                      statistic = "ks")
    )

    out_list$results$summary = GenTable(out_list$topGO_obj,
                                         classicFisher = out_list$results$fisher,
                                         classicKS = out_list$results$ks,
                                         elimKS = out_list$results$ks_elim,
                                         orderBy = "elimKS",
                                         ranksOf = "classicFisher",
                                         topNodes = go_topNodes,
                                        numChar=1000)

    #Defined colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
    # found here -- thank you stranger! https://rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot
    colMap = function(x) {
      .col = rep(rev(heat.colors(length(unique(x)))), time = table(x))
      return(.col[match(1:length(x), order(x))])
    }

    pValue.classic = score(out_list$results$ks)
    pValue.elim = score(out_list$results$ks_elim)[names(pValue.classic)]
    gstat = termStat(out_list$topGO_obj, names(pValue.classic))
    gSize = gstat$Annotated / max(gstat$Annotated) * 4
    gCol = colMap(gstat$Significant)

    diagnostic_plot_outpath =
      file.path(plot_outdir,paste0(res_name,"_classic_vs_elim.pdf"))
    pdf(diagnostic_plot_outpath)
      plot(pValue.classic,
           pValue.elim,
           xlab = "p-value classic",
           ylab = "p-value elim",
           pch = 19,
           cex = gSize,
           col = gCol)
    dev.off()

    sel.go = names(pValue.classic)[pValue.elim < pValue.classic]

    out_list$diagnostics$elim_lessThan_classic_terms =
      cbind(termStat(out_list$topGO_obj, sel.go),
            elim = pValue.elim[sel.go],
            classic = pValue.classic[sel.go])

    # out_list$go_network = showSigOfNodes(out_list$topGO_obj,
    #                                      score(out_list$results$ks_elim),
    #                                      firstSigNodes = go_topNodes,
    #                                      useInfo = 'all')

    ontology_network_plot_prefix =
      file.path(plot_outdir,paste0(res_name,"_ontology_network"))
    printGraph(out_list$topGO_obj,
               out_list$results$ks_elim,
               firstSigNodes = go_topNodes,
               fn.prefix = ontology_network_plot_prefix,
               useInfo = "all",
               pdfSW = TRUE)

    out_list
  }

  tryCatch({
    expr = getResults(out_list, go_topNodes)},
    error = function(e){
      print(paste0("error: ", e, ". Check the topGO obj. It may be that there ",
                   "are no significant genes, in which case this is correctly ",
                   "skipped."))
    },
    finally = {
      out_list
    }
  )
}

# load data --------------------------------------------------------------------

# yeast_orgdb = org.Sc.sgd.db

shrunken_res_lists = readRDS(here("data/shrunken_res_lists.rds"))

# dds_lrt = readRDS(here("data/mimic_cc_dds_lrt.rds"))

# code from james.mcdonald https://support.bioconductor.org/p/9142025/#9142035
gn2go = mapIds(org.Sc.sgd.db, keys(org.Sc.sgd.db), "GOALL", "ORF", multiVals = "list")
# filter out NA mappings
gn2go = gn2go[!sapply(gn2go, function(x) all(is.na(x)))]

# filter results and run GO enrichment -----------------------------------------

fltr_gene_list = map(shrunken_res_lists, ~map(., resultsFilterGetGenes,
                                              padj_thres = .05))

getAllTopGo = function(condition_list, condition_name, ontology_cat){

  outdir = file.path(here("data/topGO_out"), condition_name, ontology_cat)

  dir.create(outdir, recursive = TRUE)

  foreach::foreach(
    i = names(condition_list)
  ) %do% {
    createTopGOobject(condition_list[[i]],
                      gn2go,
                      ontology_cat,
                      i,
                      outdir)
  }

  # map(names(condition_list),
  #     ~createTopGOobject(condition_list[[.]],
  #                        gn2go,
  #                        ontology_cat,
  #                        .,
  #                        outdir))
}

go_results = list(
  bp = map(names(fltr_gene_list), ~getAllTopGo(fltr_gene_list[[.]], ., "BP")),
  cc = map(names(fltr_gene_list), ~getAllTopGo(fltr_gene_list[[.]], ., "CC")),
  mf = map(names(fltr_gene_list), ~getAllTopGo(fltr_gene_list[[.]], ., "MF"))
)

names(go_results$bp) = names(fltr_gene_list)
names(go_results$bp$plus_lys) = names(fltr_gene_list$plus_lys)
names(go_results$bp$minus_lys) = names(fltr_gene_list$minus_lys)
names(go_results$bp$compare_interaction) = names(fltr_gene_list$compare_interaction)

names(go_results$cc) = names(fltr_gene_list)
names(go_results$cc$plus_lys) = names(fltr_gene_list$plus_lys)
names(go_results$cc$minus_lys) = names(fltr_gene_list$minus_lys)
names(go_results$cc$compare_interaction) = names(fltr_gene_list$compare_interaction)

names(go_results$mf) = names(fltr_gene_list)
names(go_results$mf$plus_lys) = names(fltr_gene_list$plus_lys)
names(go_results$mf$minus_lys) = names(fltr_gene_list$minus_lys)
names(go_results$mf$compare_interaction) = names(fltr_gene_list$compare_interaction)

write_rds(go_results, here("data/topGO_out/go_results_objects.rds"))

# get venn set enrichment ------------------------------------------------------

venn_for_topgo = getVennResults(shrunken_res_lists, abs_lfc_thres = 0, padj_thres = 100)

venn_set_res = list(
  plus_lys = list(

    eds1_only = fltr_gene_list$plus_lys$EDS1[names(fltr_gene_list$plus_lys$EDS1) %in% c(setdiff(venn_for_topgo$genes$plus_lys$up$singles$EDS1,
                                                           venn_for_topgo$genes$plus_lys$up$singles$RGT1),
                                                   setdiff(venn_for_topgo$genes$plus_lys$down$singles$EDS1,
                                                           venn_for_topgo$genes$plus_lys$down$singles$RGT1))],

    rgt1_only = fltr_gene_list$plus_lys$RGT1[names(fltr_gene_list$plus_lys$RGT1) %in% c(setdiff(venn_for_topgo$genes$plus_lys$up$singles$RGT1,
                                                           venn_for_topgo$genes$plus_lys$up$singles$EDS1),
                                                   setdiff(venn_for_topgo$genes$plus_lys$down$singles$RGT1,
                                                           venn_for_topgo$genes$plus_lys$down$singles$EDS1))],


    eds1_rgt1_intersect_eds1Effect = fltr_gene_list$plus_lys$EDS1[names(fltr_gene_list$plus_lys$EDS1) %in% c(intersect(venn_for_topgo$genes$plus_lys$up$singles$RGT1,
                                                           venn_for_topgo$genes$plus_lys$up$singles$EDS1),
                                                             intersect(venn_for_topgo$genes$plus_lys$down$singles$RGT1,
                                                           venn_for_topgo$genes$plus_lys$down$singles$EDS1))],

    eds1_rgt1_intersect_rgt1Effect = fltr_gene_list$plus_lys$RGT1[names(fltr_gene_list$plus_lys$RGT1) %in% c(intersect(venn_for_topgo$genes$plus_lys$up$singles$RGT1,
                                                                              venn_for_topgo$genes$plus_lys$up$singles$EDS1),
                                                                    intersect(venn_for_topgo$genes$plus_lys$down$singles$RGT1,
                                                                              venn_for_topgo$genes$plus_lys$down$singles$EDS1))]
  ),

  minus_lys = list(

    eds1_only = fltr_gene_list$minus_lys$EDS1[names(fltr_gene_list$minus_lys$EDS1) %in% c(setdiff(venn_for_topgo$genes$minus_lys$up$singles$EDS1,
                                                           venn_for_topgo$genes$minus_lys$up$singles$RGT1),
                                                   setdiff(venn_for_topgo$genes$minus_lys$down$singles$EDS1,
                                                           venn_for_topgo$genes$minus_lys$down$singles$RGT1))],

    rgt1_only = fltr_gene_list$minus_lys$RGT1[names(fltr_gene_list$minus_lys$RGT1) %in% c(setdiff(venn_for_topgo$genes$minus_lys$up$singles$RGT1,
                                                           venn_for_topgo$genes$minus_lys$up$singles$EDS1),
                                                   setdiff(venn_for_topgo$genes$minus_lys$down$singles$RGT1,
                                                           venn_for_topgo$genes$minus_lys$down$singles$EDS1))],

    eds1_rgt1_intersect_eds1Effect = fltr_gene_list$minus_lys$EDS1[names(fltr_gene_list$minus_lys$EDS1) %in% c(intersect(venn_for_topgo$genes$minus_lys$up$singles$RGT1,
                                                                       venn_for_topgo$genes$minus_lys$up$singles$EDS1),
                                                             intersect(venn_for_topgo$genes$minus_lys$down$singles$RGT1,
                                                                       venn_for_topgo$genes$minus_lys$down$singles$EDS1))],

    eds1_rgt1_intersect_rgt1Effect = fltr_gene_list$minus_lys$RGT1[names(fltr_gene_list$minus_lys$RGT1) %in% c(intersect(venn_for_topgo$genes$minus_lys$up$singles$RGT1,
                                                                                                     venn_for_topgo$genes$minus_lys$up$singles$EDS1),
                                                                                           intersect(venn_for_topgo$genes$minus_lys$down$singles$RGT1,
                                                                                                     venn_for_topgo$genes$minus_lys$down$singles$EDS1))]
  )

)

venn_go_results = list(
  bp = map(names(venn_set_res), ~getAllTopGo(venn_set_res[[.]], ., "BP")),
  cc = map(names(venn_set_res), ~getAllTopGo(venn_set_res[[.]], ., "CC")),
  mf = map(names(venn_set_res), ~getAllTopGo(venn_set_res[[.]], ., "MF"))
)

names(venn_go_results$bp) = names(venn_set_res)
names(venn_go_results$bp$plus_lys) = names(venn_set_res$plus_lys)
names(venn_go_results$bp$minus_lys) = names(venn_set_res$minus_lys)

names(venn_go_results$cc) = names(venn_set_res)
names(venn_go_results$cc$plus_lys) = names(venn_set_res$plus_lys)
names(venn_go_results$cc$minus_lys) = names(venn_set_res$minus_lys)

names(venn_go_results$mf) = names(venn_set_res)
names(venn_go_results$mf$plus_lys) = names(venn_set_res$plus_lys)
names(venn_go_results$mf$minus_lys) = names(venn_set_res$minus_lys)

write_rds(venn_go_results, here("data/topGO_out/go_results_objects_by_venn_set.rds"))
