library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)
#<---------------------------------------------------------->
#10X raw data : positions, aligned fiducials, scalefactors, 
#tissue low res, tissue high res, detected tissue
#<---------------------------------------------------------->
tisPos <- read_csv("rawData/kidneyInjStates/spatial/tissue_positions_list.csv",
                   col_names = FALSE)
colnames(tisPos) <- c("barcode", "in_tissue", "array_row", "array_col",
                      "pxl_row_in_fullres", "pxl_col_in_fullres")

h5Path <- "rawData/kidneyInjStates/filtered_feature_bc_matrix.h5"
counts <- Read10X_h5(h5Path)
seuratObj <- CreateSeuratObject(counts = counts)

tisPosFilt <- tisPos %>%
  dplyr::filter(barcode %in% colnames(seuratObj))

seuratObj@meta.data <- seuratObj@meta.data %>%
  tibble::rownames_to_column(var = "barcode") %>%
  left_join(tisPosFilt, by = "barcode") %>%
  tibble::column_to_rownames(var = "barcode")
head(seuratObj@meta.data)


imagePath <- "rawData/kidneyInjStates/spatial/" 
image <- Read10X_Image(imagePath)
seuratObj[["slice1"]] <- image
seuratObj <- SCTransform(seuratObj, assay = "RNA", verbose = FALSE)


seuratObj <- RunPCA(seuratObj)
seuratObj <- RunUMAP(seuratObj, dims = 1:30)


seuratObj <- FindNeighbors(seuratObj, dims = 1:30)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)

SpatialDimPlot(seuratObj, group.by = "seurat_clusters")
SpatialFeaturePlot(seuratObj, features = c("PKHD1L1"))
