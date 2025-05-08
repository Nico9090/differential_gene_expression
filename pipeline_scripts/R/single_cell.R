#Functions
#_______________________________________________________________________________
library(dplyr)
library(Seurat)
library(patchwork)
library(biomaRt)
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
library(goseq)
library(GO.db)
library(tidyr)
#_______________________________________________________________________________

load_data<-function(data_dir,project_name){
  read.data <- Read10X(data.dir=data_dir)
  read_data <- CreateSeuratObject(counts = mut1.data, 
                                project = project_name, 
                                min.cells =3, 
                                min.features = 200)
  
}

