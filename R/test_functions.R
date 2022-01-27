library(brentlabRnaSeqTools)
library(Rsamtools)
library(tidyverse)
library(rtracklayer)
library(here)

# include some functions here for use in this notebook

# note: copy of the one in brentlabRnaSeqTools, but needs to handle more genotypes
eds_createNovoalignHtseqQcSampleSheet =
  function(meta_df, bam_prefix,
           bam_suffix = "_sorted_aligned_reads_with_annote.bam"){

    meta_df %>%
      mutate(genotype1 = str_replace(genotype1, "_", "."),
             genotype2 = str_replace(genotype2, "_", "."),
             genotype3 = str_replace(genotype3, "_", "."),
             genotype4 = str_replace(genotype4, "_", ".")) %>%
      select(fastqFileNumber, sample,
             genotype1, perturbation1, marker1,
             genotype2, perturbation2, marker2,
             genotype3, perturbation3, marker3,
             genotype4, perturbation4, marker4) %>%
      mutate(bam_path =
               file.path(bam_prefix,
                         paste0(sample, bam_suffix))) %>%
      unite(geno1, genotype1, perturbation1, marker1, sep = "_") %>%
      unite(geno2, genotype2, perturbation2, marker2, sep = "_") %>%
      unite(geno3, genotype3, perturbation3, marker3, sep = "_") %>%
      unite(geno4, genotype4, perturbation4, marker4, sep = "_") %>%
      pivot_longer(c(geno1, geno2, geno3, geno4),
                   names_to = 'perturbation_id', values_to=  "tmp2") %>%
      mutate(tmp2 = str_remove(tmp2, "__$")) %>%
      filter(tmp2 != "") %>%
      separate(tmp2, into = c('locus', 'perturbation', 'marker'), sep = "_") %>%
      mutate(locus = str_replace(locus, "\\.", "_")) %>%
      filter(locus != "NA")

  }

test_locusCoverage = function(bam_path, granges, library_strandedness,
                              index_suffix = ".bai",
                              coverage_threshold = 0, ...){

  # construct index path
  bam_index_path = paste0(bam_path, index_suffix)
  # check arguments
  for(arg in c(bam_path, bam_index_path)){
    if(!file.exists(arg)){
      stop(paste0("path not valid: ", arg, ". If this ends in .bai, then
                  it means you must first index the bam file. Make sure the
                  index is saved in the same directory as the bam"))
    }
  }

  sbp = strandedScanBamParam(granges, library_strandedness, ...)
  sbp

  # set parameter options
  p_param <- PileupParam(min_mapq = sbp@mapqFilter,
                         min_nucleotide_depth = coverage_threshold,
                         distinguish_strands=TRUE,
                         distinguish_nucleotides=FALSE)

  message("reading bam file...")
  bamfile = BamFile(bam_path, bam_index_path)

  message("calculating coverage...")
  coverage_df = pileup(bamfile,
                       scanBamParam = sbp,
                       pileupParam = p_param)

  message("done")
  length(unique(coverage_df$pos))/sum(width(granges))
}

test_geneGRanges = function(annote_obj_path, gene_id, id_col = "ID",
                            feature_col = "type", feature = "cds"){

  # TODO generalize this to 'featureGranges'

  annot_obj = readRDS(annote_obj_path)

  if(!class(annot_obj) == "GRanges"){
    stop(paste0("annote_obj_path must lead to a GRanges object. ",
                "It should be created with a command similar to the ",
                "following: rtracklayer::import('path/to/gff')"))
  }

  regions = tryCatch(
    expr = {
      annot_obj[grepl(gene_id, annot_obj@elementMetadata[,id_col]) &
                  grepl(feature, annot_obj@elementMetadata[,feature_col],
                        ignore.case = TRUE),]
    },
    error = function(e){
      message(
        paste0('geneGRanges() Error: cannot create GRanges for gene_id: ',
               gene_id))
      print(e)
    },
    warning = function(w){
      message("geneGRanges() warning: ")
      print(w)
    },
    finally = {
      # none
    }
  )

  return(regions)
}

granges = featureGRanges(Sys.getenv("yeast_gff"), "YBR033W", 'locus_tag')
bam_path = "/mnt/htcf_scratch/chasem/eds1/results/star_salmon/YKL166C_YPL203W_Minimal_HisMetLeuUra_Glucose_2_90_1NMPP7_4.markdup.sorted.bam"
library_strandedness = 'reverse'

locusCoverage(bam_path, granges, library_strandedness)
