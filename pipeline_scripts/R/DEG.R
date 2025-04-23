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
#ADD GENE NAMES_________________________________________________________________
counts_csv<-"" #should have gene names from ens_names.R
counts<-read_csv(counts_csv) 
#PCA PLOT_______________________________________________________________________
