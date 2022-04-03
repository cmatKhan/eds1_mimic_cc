# https://huoww07.github.io/Bioinformatics-for-RNA-Seq/lessons/06_Pathway_Enrichment.html

# load necessary library ggplot2
library(ggplot2)

# add another column in the results table to label the significant genes using threshold of padj<0.05 and absolute value of log2foldchange >=1
res_table <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_table <- res_table %>%
  mutate(threshold_OE =  padj < 0.05 & abs(log2FoldChange) >= 1)
# you can view the modified table
view(res_table)
# make volcano plot, the significant genes will be labeled in red
ggplot(res_table) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  scale_color_manual(values=c("black", "red")) +  # black v.s. red dots
  ggtitle("WT v.s. SNF2") +                       # this line defines the title of the plot
  xlab("log2 fold change") +                      # this line defines the name of the x-axis
  ylab("-log10 adjusted p-value") +               # name of y-axis
  scale_x_continuous(limits = c(-7.5,7.5)) +      # the axis range is set to be from -7.5 to 7.5
  theme(legend.position = "none", #c(0.9, 0.9),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

Plot top 50 significant genes in a heatmap

Sort the rows from smallest to largest padj and take the top 50 genes:

  significant_results_sorted <- significant_results[order(significant_results$padj), ]
significant_genes_50 <- rownames(significant_results_sorted[1:50, ])

We now have a list of 50 genes with most significant padj value. But we need to find the counts corresponding to these genes. To extract the counts from the rlog transformed object:

  rld_counts <- assay(rld)

Select by row name using the list of genes:

  rld_counts_sig <- rld_counts[significant_genes_50, ]

Plot multiple genes in a heatmap:

  pheatmap(rld_counts_sig,
           cluster_rows = T,
           show_rownames = T,
           annotation = meta,
           border_color = NA,
           fontsize = 10,
           scale = "row",
           fontsize_row = 8,
           height = 20)
