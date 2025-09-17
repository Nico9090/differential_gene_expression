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

library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)
#<---------------------------------------------------------->
#10X raw data : positions, aligned fiducials, scalefactors, 
#tissue low res, tissue high res, detected tissue
#<---------------------------------------------------------->
tisPos #file with the tissue positions

#plot image
ggplot(tisPos, aes(x = gridX, y = gridY,
                   color = as.factor(inTissue))) +
  geom_point(size = 2) +
  xlab("x") +
  ylab("y") +
  scale_y_reverse() + 
  coord_fixed() +
  theme_minimal() +
  scale_color_manual(
    values = c("0" = "gray70", "1" = "red3"),
    labels = c("0" = "Background", "1" = "Tissue"),
    name = "Region"
  )+
  labs(title = "Kidney tissue section")

#<------------------------------------------------------------>

counts <- Read10X_h5(h5Path) #single cell data
seuratObj <- CreateSeuratObject(counts = counts)

tisPosFilt <- tisPos %>% #keep only the barcodes from the single cell
  dplyr::filter(barcode %in% colnames(seuratObj)) #id by barcode

#removes background
#------------------------------------------------------------->


seuratObj@meta.data <- seuratObj@meta.data %>%
  tibble::rownames_to_column(var = "barcode") %>%
  left_join(tisPosFilt, by = "barcode") %>%
  tibble::column_to_rownames(var = "barcode")
head(seuratObj@meta.data)


image <- Read10X_Image(imagePath)
seuratObj[["slice1"]] <- image
seuratObj <- SCTransform(seuratObj,
                         assay = "RNA",
                         verbose = FALSE)


seuratObj <- RunPCA(seuratObj)
seuratObj <- RunUMAP(seuratObj, dims = 1:30)


seuratObj <- FindNeighbors(seuratObj, dims = 1:30)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)
saveRDS(seuratObj,"healthyF1.rds")
SpatialDimPlot(seuratObj, group.by = "seurat_clusters")
SpatialFeaturePlot(seuratObj, features = c("PKHD1L1"))

healthyMarkers <- FindAllMarkers(seuratObj,min.pct = 0,
                                 logfc.threshold = 0)

clust4Cells <- healthyMarkers %>%
  dplyr::filter(cluster == 4)

makeVolcanoPlot(degTable = clust4Cells)
geneInformation <- bitr(geneID = clust4Cells$gene,
                        fromType = "SYMBOL",toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)


clust4Cells <- clust4Cells %>%
  dplyr::filter(gene %in% geneInformation$SYMBOL)

clust4Cells$geneId <- geneInformation$ENTREZID

clust4Go <- goSeqFunction(entrezIds = clust4Cells$geneId,
                          degTable = clust4Cells)

clust4GoPlot <- clust4Go %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=numInCat, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="Category size", size="Count")








#map go terms back to genes
# First get GO-to-ENTREZ mappings
go_to_entrez <- as.list(org.Hs.egGO2ALLEGS)
genes_in_your_data <- clust4Cells$geneId  # Entrez IDs

# Keep only genes from your data
go_to_entrez_filtered <- lapply(go_to_entrez, function(x) {
  intersect(x, genes_in_your_data)
})
sig_terms <- clust4Go$category[clust4Go$over_represented_pvalue < 0.05]

# Extract genes for each significant term
genes_per_term <- go_to_entrez_filtered[sig_terms]
genes_per_term[["GO:0005902"]]

entrez_to_symbol <- mapIds(org.Hs.eg.db,
                           keys = unique(unlist(genes_per_term)),
                           keytype = "ENTREZID",
                           column = "SYMBOL",
                           multiVals = "first")

# Convert each term's gene list
genes_per_term_symbols <- lapply(genes_per_term, function(x) na.omit(entrez_to_symbol[x]))

genes_per_term_symbols[["GO:0005902"]]


L1GoTerms <- bitr(geneID = "PKHD1L1",fromType = "SYMBOL",toType = "GO",
                  OrgDb = org.Hs.eg.db)

immuneResponseGenes <- genes_per_term_symbols[["GO:0006955"]]

immuneGenes <-unname(immuneResponseGenes)


immuneData <- clust4Cells %>%
  dplyr::filter(gene %in% immuneGenes)
library(ggrepel)

immuneDataPlot <- makeVolcanoPlot(immuneData)+
  geom_text_repel(data = top,aes(x=log2FoldChange,y=-log10(padj),
                           label = gene))+
  xlab(label = "Average log2 Fold Change")+
  labs(title = "Immune Response(GO:0006955) Genes")


top <- immuneData %>%
  dplyr::arrange(padj)%>%
  slice_head(n=10)





























topHealthy <- healthyMarkers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)
print(topHealthy)


library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
healthyGoIds <- bitr(geneID = topHealthy$gene,fromType = "SYMBOL",
                     toType = "GO",OrgDb = org.Hs.eg.db)
goTerms <- AnnotationDbi::select(GO.db, keys = unique(healthyGoIds$GO),
                                  columns = c("TERM", "ONTOLOGY"),
                                  keytype = "GOID") %>% distinct()
goTable <- left_join(healthyGoIds,
                     goTerms,
                     by = c("GO" = "GOID"))















metaData <- seuratObj@meta.data %>%
  dplyr::select(seurat_clusters, pixelX, pixelY)

ggplot(metaData, aes(x = pixelX, y = pixelY, color = seurat_clusters)) +
  geom_point(size = 1, shape = 0,
             stroke = 4) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  scale_color_manual(
    values = c("0" = "red4",  
               "1" = "blue4",  
               "2" = "green4",  
               "3" = "purple4",
               "4" = "steelblue",
               "5" = "orange",
               "6" = "black"), 
    name = "cluster"
  )+
  labs(title = "Cell Populations in Healthy Kidney Tissue")







humanPKD
significant <- humanPKD %>%
  dplyr::filter(padj<=0.05)













gene <- "WNT9B"

# Extract expression values (use scaled or raw counts — here we use raw log-normalized data)
expr <- FetchData(seuratObj, vars = gene)

metaData <- cbind(seuratObj@meta.data,
                  expr = expr[[gene]])

ggplot(metaData,
       aes(x = pixelX, y = pixelY, color = expr)) +
  geom_point(size = 1, shape = 0,
             stroke = 4) +
  scale_color_viridis_c(option = "B") +  
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  labs(color = "WNT9B")


kif26b
msc
wdr86
wnt4





klk7 + tyro3+ il1rl1












l1Expr <- ggplot(metaData,
       aes(x = pixelX, y = pixelY, color = expr)) +
  geom_point(size = 1, shape = 0,
             stroke = 4) +
  scale_color_viridis_c(option = "B") +  
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  labs(color = "PKHD1L1")








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
