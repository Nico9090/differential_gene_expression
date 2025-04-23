#LIBRARIES______________________________________________________________________
#_______________________________________________________________________________
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
#FASTQs are aligned to gtf 
counts_csv<-"" #should have gene names from ens_names.R
counts<-read_csv(counts_csv) 
#PCA PLOT_______________________________________________________________________
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
}
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
