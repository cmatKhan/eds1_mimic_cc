# Purpose

This is a Rstudio project directory containing the analysis of the RNAseq data 
for the scerevisiae EDS1 'Mimic Calling Cards' condition experiment. You may 
pull this directory. If you have R >= 4.0.0, a correspondingly current version 
of Rstudio, and the packages listed at the top of any of the scripts you wish 
to run, then you can use the scripts in the `R` directory, and the notebooks in 
the `notebooks` directory, to generate the same analyses, including the plots 
in the `plots` directory. Below is a description of the directory structure:

# plots

The `plots` directory contains the following subdirectories

```
/plots
├── [4.0K]  compare_interaction
│   ├── [4.0K]  EDS1_LYS14_vs_EDS1
│   ├── [4.0K]  EDS1_LYS14_vs_EDS1_RGT1
│   ├── [4.0K]  EDS1_LYS14_vs_LYS14
│   ├── [4.0K]  EDS1_LYS14_vs_RGT1
│   ├── [4.0K]  EDS1_RGT1_vs_EDS1
│   ├── [4.0K]  EDS1_RGT1_vs_EDS1_LYS14
│   ├── [4.0K]  EDS1_RGT1_vs_LYS14
│   ├── [4.0K]  EDS1_RGT1_vs_RGT1
│   ├── [4.0K]  EDS1_vs_EDS1_LYS14
│   ├── [4.0K]  EDS1_vs_EDS1_RGT1
│   ├── [4.0K]  EDS1_vs_LYS14
│   ├── [4.0K]  EDS1_vs_RGT1
│   ├── [4.0K]  LYS14_vs_EDS1
│   ├── [4.0K]  LYS14_vs_EDS1_LYS14
│   ├── [4.0K]  LYS14_vs_EDS1_RGT1
│   ├── [4.0K]  LYS14_vs_RGT1
│   ├── [4.0K]  RGT1_vs_EDS1
│   ├── [4.0K]  RGT1_vs_EDS1_LYS14
│   ├── [4.0K]  RGT1_vs_EDS1_RGT1
│   └── [4.0K]  RGT1_vs_LYS14
├── [4.0K]  minus_lys
│   ├── [4.0K]  EDS1
│   ├── [4.0K]  EDS1_LYS14
│   ├── [4.0K]  EDS1_RGT1
│   ├── [4.0K]  LYS14
│   ├── [4.0K]  lysine_pathway_kegg
│   └── [4.0K]  RGT1
└── [4.0K]  plus_lys
    ├── [4.0K]  EDS1
    ├── [4.0K]  EDS1_LYS14
    ├── [4.0K]  EDS1_RGT1
    ├── [4.0K]  LYS14
    ├── [4.0K]  lysine_pathway_kegg
    └── [4.0K]  RGT1
```

The top level directories represent the condition. The subdirectories of each 
condition, eg EDS1, are results generated for a specific genotype. 

### compare_interaction

This contains KEGG pathway and STRINGdb visualizations for results tables which 
compare the interaction term for -lys condition. Paraphrased from the DESeq2::results
documentation: the interaction term for the -lys condition specific effect in 
del EDS1 vs del LYS14. This tests if the whether the -lys effect is different 
in del EDS1 compared to del LYS14.

### minus/plus lys

This contains KEGG pathway and STRINGdb visualizations for results tables of the 
condition specific effect of a given genotype, eg del EDS1. Paraphrased from the 
DESeq2::results documentation: this is the main effect (-lys) *plus* the 
interaction term. __NOTE__: my intent in creating these was to use them to compare 
the images in minus_lys to the same image in plus_lys. Another way of doing 
this would be to use the interaction term alone, but then that would not show 
the genes already affected by -lys.

__NOTE__: `minus_lys` and `plus_lys` each have a subdirectory, `lysine_pathway_kegg`, 
which contains the KEGG pathway PDF, with each node (gene) in the pathway highlighted 
by log2FC. Please note that the identifiers in the nodes are KEGG gene identifiers -- 
I am working on replacing those with gene IDs, though that may have to be done by 
hand if there are any pathway diagrams which are deemed 'evidence worthy' to support 
a hypothesis.

# R

R scripts used to generate items in the `data` and `plots` directories, and any 
supplementary functions to notebooks.

# docs

This is the directory served by github pages -- if and when I produce a notebook 
to display and summarize results, it will be rendered in this directory and 
served.

# data

This contains the data used in the `R` and `notebooks` directories -- all data 
necessary to generating plots, etc is here (including the DESeqDataObject and 
the results tables generated in `R/generate_results_tables.R`)

# notebooks

R Notebook(s) intended to display/summarize the findings. One of the notebooks 
will be used to generate the served documents in `docs`.
