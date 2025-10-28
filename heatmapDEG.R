rawCounts <- read_csv("../cnric/counts/Yale_gene_counts_with_gene_names.csv")
names(del34LiverCounts)
del34LiverCounts <- rawCounts[,c(1,11,12,13,14,17,18,20,21,22,27)]
onlyCounts <- del34LiverCounts[,-c(1,11)]
meta_data <- data.frame(row.names = colnames(onlyCounts),
                        genotype = c(rep("WT",4),"Pkhd1 del 3-4/del 3-4","WT",
                                     rep("Pkhd1 del 3-4/del 3-4",3)))
meta_data$genotype <- factor(meta_data$genotype,
                             levels = c("WT","Pkhd1 del 3-4/del 3-4"))
dds <- DESeqDataSetFromMatrix(countData =onlyCounts,
                              colData = meta_data,
                              design= ~genotype )
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- as.data.frame(norm_counts)

norm_counts$gene_name <- del34LiverCounts$external_gene_name

toLabel <- go_dataMouse %>% dplyr::filter(name_1006 == "DNA binding") %>% dplyr::pull(external_gene_name)
sigCounts <- norm_counts %>%
  dplyr::filter(gene_name %in% sigDel34Liver$gene_name) %>%
  dplyr::filter(!is.na(gene_name))
#START
top200_counts <- sigCounts %>%
  dplyr::filter(gene_name %in% toLabel) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::distinct(gene_name,.keep_all = TRUE)


nrow(top200_counts)
rownames(top200_counts) <- top200_counts$gene_name
annot_col <- tibble(condition = factor(c(rep("Pkhd1 del 3-4/del 3-4",1),rep("WT",3),
                                         rep("Pkhd1 del 3-4/del 3-4",2),
                                         rep("WT",3)))
)
annot_col <- as.data.frame(annot_col)
mat <- as.matrix(top200_counts[ , -10])   # remove gene name column
mat <- apply(mat, 2, as.numeric)          # ensure all columns are numeric
rownames(annot_col) <- colnames(mat)





library(pheatmap)


rownames(mat) <- top200_counts$gene_name 
colnames(mat)
colnames(mat) <- c("Pkhd1 del 3-4/del 3-4 rep 1",
                   "WT rep 1",
                   "WT rep 2",
                   "WT rep 3",
                   "Pkhd1 del 3-4/del 3-4 rep 2",
                   "Pkhd1 del 3-4/del 3-4 rep 3",
                   "WT rep 4",
                   "WT rep 5",
                   "WT rep 6")
rownames(annot_col) <- colnames(mat)
mat_log <- log2(mat + 1) 
my_colors <- list(
  condition = c(
    "Pkhd1 del 3-4/del 3-4" = "steelblue",
    "WT" = "red4"
  )
)

my_palette <- colorRampPalette(c("white","orange","red"))(10)
pheatmap(mat_log[101:125,],
         annotation_col = annot_col,
         annotation_colors = my_colors,
         color = my_palette,
         labels_row = rownames(mat),
         main = "Pkhd1 del 3-4/del 3-4 Liver Top DEGs\nGenes in cell component: nucleus",
         fontsize_row = 6,
         fontsize_col = 8,
         angle_col = 45,
         cluster_cols = TRUE)


 # add space after rows 20 and 40


ensembl2 <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # change dataset for other species

# 2. Specify your gene (symbol or Ensembl ID)
genes <- top200_counts$gene_name  # replace with your gene of interest

# 3. Get GO terms
go_dataMouse <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "namespace_1003"),
  filters = "external_gene_name",
  values = genes,
  mart = ensembl2
)

View(go_dataMouse)
