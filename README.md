# Plots

The plots directory contains the following subdirectories

```
/home/chase/projects/eds1/mimic_cc_rnaseq/plots
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
