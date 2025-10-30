library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
rawCounts <- read_delim("GSE242476_Raw_Gene_Counts.txt")
head(rawCounts,1)
onlyCounts <- rawCounts[,-1]
metaData <- tibble(samples = colnames(onlyCounts),
                   condition = c(rep("normal",4),
                                 rep("ARPKD",4)))
metaData$condition <- factor(metaData$condition,
                   levels = c("normal","ARPKD"))
dds <- DESeqDataSetFromMatrix(countData =onlyCounts,
                              colData = metaData,
                              design= ~condition )
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- as.data.frame(norm_counts)
norm_counts$geneId <- rawCounts$gene_id
all(rownames(dds) == rownames(onlyCounts))

geneInformation <- bitr(geneID = norm_counts$geneId,
                        fromType = "ENSEMBL",toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db)
nrow(geneInformation)
nrow(norm_counts)
norm_counts <- norm_counts %>%
  dplyr::filter(geneId %in% geneInformation$ENSEMBL)
geneInformation <- geneInformation %>%
  dplyr::distinct(ENSEMBL,.keep_all = TRUE)
all(rownames(geneInformation) == rownames(norm_counts))
norm_counts$gene_name <- geneInformation$SYMBOL


degTable <- read_csv("human_arpkd_deg.csv")
sigGenes <- degTable %>%
  dplyr::filter(padj <= 0.05,abs(log2FoldChange)>0.58)

downGenes <- sigGenes %>%
  dplyr::filter(log2FoldChange <0)

upGenes <- sigGenes %>%
  dplyr::filter(log2FoldChange >0)
#-
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")+
  labs(color = "condition")+
  ggtitle("Human ARPKD Kidney PCA")+
  theme(axis.text = element_text(face = "bold",
                                 size =10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",
                                  size =15))

upCounts <- norm_counts %>%
  dplyr::filter(gene_name %in% upGenes$gene_name) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::distinct(gene_name,.keep_all = TRUE) 


enrichGOHelperFunction <- function(geneSymbols,orgDb,ont){
  enrichGO(gene = geneSymbols,
           OrgDb = orgDb,
           keyType = "SYMBOL",
           ont = ont) %>%
    as.data.frame() %>%     
    as_tibble() %>%
    dplyr::select(Description, p.adjust, geneID) %>%
    dplyr::filter(p.adjust <= 0.05)
}
upGo <- enrichGOHelperFunction(geneSymbols = upCounts$gene_name,
                               orgDb = org.Hs.eg.db,ont = "BP")


pathway <- upGo %>%
  dplyr::filter(Description == "axonogenesis") %>%
  dplyr::select(geneID) %>%
  str_split("/") %>%
  unlist()

pathwayDeg <- upGenes %>%
  dplyr::filter(gene_name %in% pathway)

pathwayCounts <- upCounts %>%
  dplyr::filter(gene_name %in% pathway)
downCounts <- norm_counts %>%
  dplyr::filter(gene_name %in% downGenes$gene_name) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::distinct(gene_name,.keep_all = TRUE) 

tfCounts <- norm_counts %>%
  dplyr::filter(gene_name %in% matchedTF$gene_name)

rownames(upCounts) <- upCounts$gene_name
rownames(downCounts) <- downCounts$gene_name
rownames(pathwayCounts) <- pathwayCounts$gene_name
rownames(tfCounts) <- tfCounts$gene_name
annot_col <- tibble(condition = factor(c(rep("normal",4),rep("ARPKD",4)))
)
annot_col <- as.data.frame(annot_col)
mat <- as.matrix(tfCounts[ , -c(9,10)])   # remove gene name column
mat <- apply(mat, 2, as.numeric)          # ensure all columns are numeric
rownames(annot_col) <- colnames(mat)
library(pheatmap)

rownames(mat) <- tfCounts$gene_name 
colnames(mat) <- c("Normal 5 monthss",
                   "Normal 18 months",
                   "Normal 9 years",
                   "Normal 4 years",
                   "ARPKD 9 days",
                   "ARPKD 5 years rep 1",
                   "ARPKD 5 years rep 2",
                   "ARPKD 10 years")
rownames(annot_col) <- colnames(mat)
mat_log <- log2(mat + 1) 
my_colors <- list(
  condition = c(
    "ARPKD" = "steelblue",
    "normal" = "red4"
  )
)
my_palette <- colorRampPalette(c("white","orange","red"))(10)

pheatmap(mat_log,
         annotation_col = annot_col,
         annotation_colors = my_colors,
         color = my_palette,
         labels_row = rownames(mat_log),
         main = "ARPKD Human Kidney DEGs\n FGA transcription factors",
         fontsize_row = 6,
         fontsize_col = 8,
         angle_col = 45,
         cluster_cols = TRUE)
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", "description"),
  filters = "hgnc_symbol",
  values = sigGenes$gene_name,
  mart = mart
)

# Look for terms like "transcription factor" in description
tf_annot <- annotations[grep("transcription factor", annotations$description, ignore.case = TRUE), ]

annotations$is_receptor <- grepl("receptor", annotations$description, ignore.case = TRUE)
annotations$is_ligand <- grepl("ligand", annotations$description, ignore.case = TRUE)
receptors <- annotations %>%
  dplyr::filter(is_receptor == TRUE)


ligands <- annotations %>%
  dplyr::filter(is_ligand == TRUE)


# 2. Specify your gene (symbol or Ensembl ID)
genes <- receptors$hgnc_symbol  

# 3. Get GO terms
go_data <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "namespace_1003"),
  filters = "external_gene_name",
  values = genes,
  mart = mart
)

View(go_data)
biological_processes <- go_data %>%
  dplyr::filter(namespace_1003 %in% c("biological_process","molecular_function"))

















library(dplyr)
library(ggplot2)

top_go <- biological_processes %>%
  group_by(name_1006, namespace_1003) %>%
  summarise(n_genes = n_distinct(external_gene_name)) %>%
  arrange(desc(n_genes)) %>%
  slice_head(n = 15)

top_go <- top_go %>%
  dplyr::filter(n_genes > 10)
ggplot(top_go, aes(x = reorder(name_1006, n_genes), y = n_genes,fill = namespace_1003)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(x = "GO Term", y = "Number of Genes", title = "Top GO Terms by Gene Count") +
  theme(axis.text.y = element_text(size = 10))
upGo <- upGo %>%
  dplyr::mutate(count = sapply(str_split(geneID, "/"), length))
upGo %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::slice_head(n = 25) %>%
  ggplot()+
  geom_point(mapping = aes(x = count,y = Description,
                           colour = p.adjust,size = count))+
  labs(color = "p-adjust",size = "number of genes")+
  xlab(label = "number of genes")+
  scale_color_gradient(low = "orange",high = "red")+
  ggtitle("Human ARPKD Kidney",
          subtitle = "Go terms of upregulated genes")+
  theme(axis.text = element_text(face = "bold",
                                 size =10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",
                                  size =15))



library(enrichR)
dbs <- c("ChEA_2022", "ENCODE_TF_ChIP-seq_2023")
results <- enrichr("FGA", dbs)

# Example: get ChEA_2022 results
df <- results[["ChEA_2022"]]
df <- df %>%
  dplyr::mutate(gene = sapply(str_split(Term, " "), `[`, 1))
head(df)



results[["ChEA_2022"]] %>%
  mutate(gene = sapply(str_split(Term, " "), `[`, 1)) %>%
  ggplot(aes(x = reorder(gene, Adjusted.P.value), y = 1,
             color = Adjusted.P.value)) +
  scale_color_gradient(low = "orange",high = "purple")+
  coord_flip() +
  labs(
    x = "Transcription Factor",
    y = "p-adjusted",
    title = "TFs predicted to regulate ADAMTSL3"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())


matchedTF <- sigGenes %>%
  dplyr::filter(gene_name %in% df$gene)

library(ggrepel)
sigGenes %>%
  ggplot()+
  geom_point(mapping = aes(x = log2FoldChange,
                           y = -log10(padj),
                           color = log2FoldChange))+
  geom_text_repel(data = matchedTF,
                  mapping = aes(x = log2FoldChange,
                                y = -log10(padj),
                                label = gene_name),
                  size = 6)+
  scale_color_gradient(low = "steelblue",high = "red")+
  ylim(0,10)+
  xlim(-10,10)+
  xlab(label = "Log2-FoldChange")+
  ylab(label = "-Log10-P-adj")+
  ggtitle("Human ARPKD DEGs",
          subtitle = "Transcription Factors of FGA")+
  theme(axis.text = element_text(face = "bold",
                                 size =10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",
                                  size =15))






