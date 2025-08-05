#pre-process
library(rhdf5)
library(Matrix)
library(Seurat)

setwd("~/omics_DS/organoid/") #directory with h5 files
sc_data_paths<-dplyr::tibble(ko="GSM7821874_E17-filtered_feature_bc_matrix.h5",
                             wt="GSM7821876_C15-filtered_feature_bc_matrix.h5")
#KO------------>
h5ls(sc_data_paths[["ko"]]) #see the .h5 paths
sc_ko_barcodes<-h5read(sc_data_paths[["ko"]],"/matrix/barcodes")#bar codes
sc_ko_data<-h5read(sc_data_paths[["ko"]],"/matrix/data")#data
sc_ko_features<-h5read(sc_data_paths[["ko"]],"/matrix/features") #features
sc_ko_feature_id<-h5read(sc_data_paths[["ko"]],"/matrix/features/id") #gene ids
sc_ko_feature_gene_name<-h5read(sc_data_paths[["ko"]],"/matrix/features/name") #gene names

sc_ko_indices<-h5read(sc_data_paths[["ko"]],"/matrix/indices") #gene names
sc_ko_shape<-h5read(sc_data_paths[["ko"]],"/matrix/shape")
sc_ko_indptr <- h5read(sc_data_paths[["ko"]], "/matrix/indptr")

dcg_sc_matrix<-new("dgCMatrix",
                   Dim = as.integer(sc_ko_shape),
                   x = as.double(sc_ko_data),
                   i = as.integer(sc_ko_indices),
                   p = as.integer(sc_ko_indptr)) #input for seurat object


sc_ko_barcodes <- as.character(sc_ko_barcodes)
sc_ko_feature_gene_name <- as.character(sc_ko_feature_gene_name)


rownames(dcg_sc_matrix) <- make.unique(sc_ko_feature_gene_name)
colnames(dcg_sc_matrix) <- sc_ko_barcodes

#make seurat object
seurat_obj <- CreateSeuratObject(counts = dcg_sc_matrix)

#WT---------------------->
h5ls(sc_data_paths[["wt"]]) 
sc_wt_shape<-h5read(sc_data_paths[["wt"]],"/matrix/shape")
sc_wt_data<-h5read(sc_data_paths[["wt"]],"/matrix/data")#data
sc_wt_indices<-h5read(sc_data_paths[["wt"]],"/matrix/indices") #gene names
sc_wt_indptr <- h5read(sc_data_paths[["wt"]], "/matrix/indptr")
sc_wt_barcodes<-h5read(sc_data_paths[["wt"]],"/matrix/barcodes")
sc_wt_gene_name<-h5read(sc_data_paths[["wt"]],"/matrix/features/name")

dcg_sc_wt_matrix<-new("dgCMatrix",
                      Dim = as.integer(sc_wt_shape),
                      x = as.double(sc_wt_data),
                      i = as.integer(sc_wt_indices),
                      p = as.integer(sc_wt_indptr)) #input for seurat object

sc_ko_barcodes <- as.character(sc_wt_barcodes)
sc_ko_feature_gene_name <- as.character(sc_wt_gene_name)


rownames(dcg_sc_wt_matrix) <- make.unique(sc_wt_gene_name)
colnames(dcg_sc_wt_matrix) <- sc_wt_barcodes

#this works 
colnames(dcg_sc_wt_matrix) <- make.unique(sc_wt_barcodes)
rownames(dcg_sc_wt_matrix) <- make.unique(sc_wt_gene_name)
#make seurat object
wt_seurat_obj <- CreateSeuratObject(counts = dcg_sc_wt_matrix)
