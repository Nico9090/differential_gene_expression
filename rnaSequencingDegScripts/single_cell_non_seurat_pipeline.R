library(tidyverse)
library(Seurat)
sc_thyroid <- read_delim("GSE183963_Single_cell_mouse_ndata.csv",
                         delim = ";",
                         locale = locale(decimal_mark = ","))

fpc_genes <- sc_thyroid %>%
  dplyr::filter(`...1` %in% c("Pkhd1l1","Pkhd1","Actb","Gapdh"))

rownames(sc_thyroid)<- sc_thyroid$`...1`
matrix <- sc_thyroid[,-1]
rownames(matrix) <- rownames(sc_thyroid)
mat <- as.matrix(matrix)
head(rownames(mat))

mat_t <- t(mat)               # now cells × genes

pca <- prcomp(mat_t, scale. = TRUE)
saveRDS(pca,"mouse_thyroid_sc_pca.rds")
# View results
summary(pca)
plot(pca$x[,1], pca$x[,2], pch=16)

cell_ids <- colnames(mat)
samples <- sub("_(.*)", "", cell_ids)   # take everything before underscore
plot(pca$x[,1], pca$x[,2], 
     col = as.factor(samples), 
     pch = 16,
     xlab="PC1", ylab="PC2")
legend("topright", legend=unique(samples), col=1:length(unique(samples)), pch=16)



pca_coords <- pca$x
pc_use <- pca_coords[, 1:20]  # first 20 PCs
set.seed(123)
k <- 5  # choose number of clusters
clusters <- kmeans(pc_use, centers = k)$cluster

# clusters is a vector assigning each cell to a cluster
table(clusters)
plot(pc_use[,1], pc_use[,2], col=clusters, pch=16,
     xlab="PC1", ylab="PC2")
legend("topright", legend=1:k, col=1:k, pch=16)


cluster5_cells <- names(clusters[clusters == 5])
  



other_cells    <- names(clusters[clusters != 5])

avg_expr_cluster <- rowMeans(mat[, cluster5_cells])
avg_expr_other    <- rowMeans(mat[, other_cells])
logfc <- log2((avg_expr_cluster + 1e-6) / (avg_expr_other + 1e-6))
pvals <- apply(mat, 1, function(gene_expr) {
  wilcox.test(gene_expr[cluster5_cells],
              gene_expr[other_cells])$p.value
})

p_adj <- p.adjust(pvals, method = "BH")  # Benjamini-Hochberg
markers <- data.frame(
  gene = rownames(mat),
  avg_log2FC = logfc,
  p_val = pvals,
  p_val_adj = p_adj
)

markers <- markers[order(-markers$avg_log2FC), ] #for cluster 1
head(markers)
write.csv(markers, "mouse_thyroid_cluster5_markers.csv", row.names = FALSE)





cell_cluster3 <- read_csv("mouse_thyroid_cluster3_markers.csv")
cell_cluster1 <- read_csv("mouse_thyroid_cluster1_markers.csv")

cell_cluster2 <-read_csv("mouse_thyroid_cluster2_markers.csv") #upregulated pkhd1l1
cell_cluster5 <- read_csv("mouse_thyroid_cluster5_markers.csv") #upregulated pkhd1l1
cell_cluster4 <- read_csv("mouse_thyroid_cluster4_markers.csv") #upregulated pkhd1l1

