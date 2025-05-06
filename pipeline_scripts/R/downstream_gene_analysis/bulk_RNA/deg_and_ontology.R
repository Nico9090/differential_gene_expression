#LIBRARIES______________________________________________________________________
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
suppressWarnings(suppressMessages(library(tximport)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
suppressWarnings(suppressMessages(library(AnnotationDbi)))
suppressWarnings(suppressMessages(library(GOplot)))
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(pheatmap)))
suppressWarnings(suppressMessages(library(extrafont)))
suppressWarnings(suppressMessages(library(sva)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(fgsea)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(EnhancedVolcano)))
library(biomaRt)
library(goseq)
library(dplyr)
#_______________________________________________________________________________
#LOAD DATA
counts<-read_csv("") #add the file with the counts
#_______________________________________________________________________________
#PCA set-up_____________________________________________________________________
make_PCA<-function(counts_data,design,group,...){
  counts_data<-counts_data %>% select(where(is.numeric))
  meta_data_variables<-list(...) #add variables such as diseased states and sex
  meta_data<- as.data.frame(meta_data_variables)
  rownames(meta_data)<-colnames(counts_data)
  dds <- DESeqDataSetFromMatrix(countData =counts_data,
                              colData = meta_data,
                              design= ~ design )
  dds <- DESeq(dds)
  resultsNames(dds)
  vsd <- vst(dds, blind=TRUE)
  pl <- plotPCA(vsd, 
              intgroup=c(group,"sizeFactor"), 
              returnData = TRUE)
  percentVar <- round(100 * attr(pl, "percentVar"))
  return(list(
    dds=dds,
    pl=pl,
    percentVar=percentVar))
} #do this first before plotting
#________________________________________________________________________________
#PCA PLOT________________________________________________________________________
plot_PCA<-function(pl,percentVar,title){
  ggplot(pl,aes(PC1,PC2,
                shape=sex,
                pointsize = 5)) + #plot
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(title) +
  theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 25,  face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )
}
#_____________________________________________________________________
#VOLCANO PLOT_________________________________________________________
plot_Volcano<-function(dds,title){
  ### Deseq2 differential expression calculation
  res <- results(dds)#, name="condition_WT_vs_L1")
  res_df <- as.data.frame(res)
  res_df$gene_name <- counts$external_gene_name
  res_df$Gene_ID<-counts$Gene_ID
  res_df$color <- ifelse(res_df$log2FoldChange > 0, "blue", 
                       ifelse(res_df$log2FoldChange < 0, "red", "gray"))
  top_genes <- head(res_df[order(res_df$padj), ], 10) #change to any number of wanted genes
  ggplot(res_df,
       aes(x = log2FoldChange, 
           y = -log10(pvalue))
       ) +
  geom_point(aes(color = color), 
             alpha = 0.6) +
  scale_color_manual(values = c("red", "blue", "gray")) +
  theme_minimal() +
  labs(title = title, 
       x = "log2 Fold Change", 
       y = "-log10 p-value"
       ) +
  theme(axis.title.x = element_text(size = 18, face = "bold", family = "Courier"),
        axis.title.y = element_text(size = 18, face = "bold", family = "Courier"),
        axis.text.x = element_text(size = 14, face = "plain", family = "Courier"),
        axis.text.y = element_text(size = 14, face = "plain", family = "Courier"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none")+ 
  geom_text_repel(data = top_genes, aes(label = top_genes$gene_name, 
                                        size = 3,
                                        family = "Courier"))
  return(res_df)
}
#______________________________________________________________________________________
#Gene Ontology Set-up__________________________________________________________________
DEG_table<-write_csv(res_df,"") #create csv file of DEG table
plot_GO<-function(ensembl_genome,table_DEG){
  ensembl <- useMart("ensembl", 
                     dataset = ensembl_genome)
  lengths <- getBM(
    attributes = c("ensembl_gene_id", 
                   "start_position", 
                   "end_position"
                  ),
  filters = "ensembl_gene_id",
  values =table_DEG$Gene_ID,
  mart = ensembl
  )
  lengths$dist=lengths$end_position-lengths$start_position
  table_DEG$length<-lengths$dist
  sigData <- as.integer(!is.na(table_DEG$padj) & table_DEG$padj < 0.01)
  names(sigData) <- table_DEG$Gene_ID
  pwf <- nullp(sigData, "mm10", "ensGene", bias.data = table_DEG$length)#probability weight function
  goResults <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))
  #PLOT
  goResults %>% 
    top_n(10, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, 
               y=term, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
        geom_point() +
        expand_limits(x=0) +
        labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
}
