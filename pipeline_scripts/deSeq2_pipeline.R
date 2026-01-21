raw <- read_csv("")
rs <- rowSums(raw[ , -c(1,27)], na.rm = TRUE)
raw_filtered <- raw[rs > 10, ]
only_counts <- raw_filtered[,-c(1,27)]
meta_data <- tibble(samples = names(only_counts),
                    genotype = c(rep("L1",5),"D1",rep("WT",7),
                                 rep("D1",3),"WT","L1",
                                 rep("D1",3),"L1",rep("WT",3),
                                 "L1"),
                    sex = c("M","F",rep("M",4),"F","M",
                            rep("F",3),"M","F",rep("M",4),
                            "F",rep("M",3),"F","M","F",
                            rep("M",2)
                    ),
                    tissue = c(rep("kidney",2),"liver",
                               rep("kidney",6),rep("liver",4),
                               rep("kidney",2),rep("liver",6),
                               rep("kidney",4),"liver"))
meta_data$genotype <- factor(
  meta_data$genotype,
  levels = c("WT", "D1", "L1")
)
meta_data$sex <- factor(meta_data$sex,
                        levels = c("F","M"))
meta_data$tissue <- factor(meta_data$tissue,
                           levels = c("kidney","liver"))

dds <- DESeqDataSetFromMatrix(countData =only_counts,
                              colData = meta_data,
                              design= ~ genotype + sex + tissue)
#correct for sex and tissue?
covariates <- model.matrix(~ 0 + sex + tissue, colData(dds))
vsd <- vst(dds, blind=FALSE) # accounts for design
vsd_corrected <- vsd
assay(vsd_corrected) <- removeBatchEffect(
  assay(vsd),
  covariates = covariates
)
pca <- prcomp(
  t(assay(vsd_corrected)),
  scale. = FALSE
)

loadings <- pca$rotation
loadings <- as.data.frame(loadings)
loadings$Gene_ID <- raw_filtered$Gene_ID
loadings$external_gene_name <- raw_filtered$external_gene_name

scores   <- pca$x

pc1 <- loadings[, "PC1"]

cutoff <- quantile(abs(pc1), 0.95)  # top 5%
loadings_filt <- loadings %>%
  dplyr::filter(abs(PC1) >= cutoff)

pos_genes <- loadings_filt %>%
  dplyr::filter(PC1 > 0)
neg_genes <- loadings_filt %>%
  dplyr::filter(PC1<0)
go_pos <- enrichGOHelperFunction(geneSymbols = pos_genes$external_gene_name,
                                 orgDb = org.Mm.eg.db,
                                 ont="BP")

go_neg <- go_neg %>%
  dplyr::mutate(count = sapply(str_split(geneID, "/"), length))

go_pos %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::slice_head(n = 25) %>%
  ggplot()+
  geom_point(mapping = aes(x = GeneRatio,y = Description,
                           colour = p.adjust,size = count))+
  labs(color = "p-adjust",size = "number of genes")+
  xlab(label = "pathway ratio")+
  scale_color_gradient(low = "orange",high = "red")+
  ggtitle("PC1 gene pathways of top 5%",
          subtitle = "")+
  theme(axis.text = element_text(face = "bold",
                                 size =10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",
                                  size =15))


go_neg <- enrichGOHelperFunction(geneSymbols = neg_genes$external_gene_name,
                                 orgDb = org.Mm.eg.db,
                                 ont="BP")

#correlation matrix
expr_matrix <- assay(vsd_corrected)
cor_mat <- cor(expr_matrix, method = "pearson")
heatmap(cor_mat, symm = TRUE)

pl <- plotPCA(vsd_corrected, 
              intgroup=c("genotype","sex","tissue"), 
              returnData = TRUE,
              pcsToUse = c(1,2))
percentVar <- round(100 * attr(pl, "percentVar"))
ggplot(pl, aes(PC1, PC2,
               color = genotype,
               shape = sex)) +
  geom_point(size = 4,stroke = 2) +
  geom_text_repel(aes(label = name),
                  size = 2)+
  stat_ellipse(aes(group = tissue),
               linetype = "dashed",
               color = "grey40")+
  xlab(paste0("PC1: (",percentVar[1],"%)")) +
  ylab(paste0("PC2: (",percentVar[2],"%)")) +
  scale_color_manual(values = c("D1" = "red",
                                "L1" = "steelblue",
                                "WT" = "green3"))+
  scale_shape_manual(
    values = c("F" = 21, "M" = 22)
  ) +
  #facet_wrap(~ tissue,scales = "fixed") +
  theme(
    plot.title = element_text(size = 25,face = "bold"),
    axis.title.x = element_text(size = 15,face = "bold"),
    axis.title.y = element_text(size = 15,face = "bold"),
    axis.text.x = element_text(size = 10,face = "bold"),
    axis.text.y = element_text(size = 10,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    panel.background = element_rect(fill = "white")
  )


dds <- DESeq(dds)
resultsNames(dds)
res_D1 <- results(dds, contrast = c("genotype", "D1", "WT"))
res_L1 <- results(dds, contrast = c("genotype", "L1", "WT"))

d1_vs_WT_res <- as.data.frame(res_D1) 
d1_vs_WT_res$gene_name <- raw_filtered$external_gene_name
d1_vs_WT_res$gene_id<- raw_filtered$Gene_ID
d1_vs_WT_res <- d1_vs_WT_res %>%
  dplyr::left_join(y = raw_filtered,
                   by = c("gene_id"="Gene_ID"))




d1_vs_WT_res <- read_csv("corForSexandTissue_D1WT_deg.csv")
sig <- d1_vs_WT_res %>%
  dplyr::filter(padj <= 0.05,abs(log2FoldChange)>=1)



L1_vs_WT_res <- as.data.frame(res_L1) 
L1_vs_WT_res$gene_name <- raw_filtered$external_gene_name
L1_vs_WT_res$gene_id<- raw_filtered$Gene_ID
L1_vs_WT_res <- L1_vs_WT_res %>%
  dplyr::left_join(y = raw_filtered,
                   by = c("gene_id"="Gene_ID"))

L1_vs_WT_res <- read_csv("corForSexandTissue_L1WT_deg.csv")
sig <- L1_vs_WT_res %>%
  dplyr::filter(padj <= 0.05,abs(log2FoldChange)>=1)




top_d1 <- d1_vs_WT_res %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10)
ggplot()+
  geom_point(data = d1_vs_WT_res,
             mapping=aes(x=log2FoldChange,
                         y=-log10(padj)))+
  geom_text_repel(data = top_d1,
                  mapping=aes(x=log2FoldChange,
                              y=-log10(padj),
                              label = gene_name))+
  ggtitle("")+
  labs(subtitle = "")+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 15,face = "bold"),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),  
  )+
  xlab(label = "Log2-FoldChange")+
  ylab(label = "p-adj")


top_L1 <- L1_vs_WT_res %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10)
ggplot()+
  geom_point(data = L1_vs_WT_res,
             mapping=aes(x=log2FoldChange,
                         y=-log10(padj)))+
  geom_text_repel(data = top_L1,
                  mapping=aes(x=log2FoldChange,
                              y=-log10(padj),
                              label = gene_name))+
  ggtitle("")+
  labs(subtitle = "")+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 15,face = "bold"),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),  
  )+
  xlab(label = "Log2-FoldChange")+
  ylab(label = "p-adj")











#gene weights to PCA
pca <- prcomp(t(expr_matrix))
loadings <- as.data.frame(pca$rotation)
loadings$Gene_ID <- raw_filtered$Gene_ID
loadings$external_gene_name <- raw_filtered$external_gene_name

write_csv(loadings,"all_d1L1WT_kidneyLiver_PCs.csv")
pcs <- read_csv("allSamples_d1L1WT_kidneyLiver_PCs.csv")
base::sum(pcs$PC1)
weights_genes <- loadings[order(abs(loadings$PC5),
                              decreasing = TRUE)[1:20], ]

pc5 <- loadings[, "PC5"]

cutoff <- quantile(abs(pc5), 0.99)  # top 1%

weights_genes <- loadings %>%
  dplyr::filter(abs(PC5)>=cutoff)
library(org.Mm.eg.db)
weighted_genes_pathways <- enrichGOHelperFunction(geneSymbols = sig$external_gene_name,
                                                  orgDb = org.Mm.eg.db,ont = "BP")

weighted_genes_pathways <- weighted_genes_pathways %>%
  dplyr::mutate(count = sapply(str_split(geneID, "/"), length))

enrichGOHelperFunction <- function(geneSymbols, orgDb, ont) {
  enrichGO(
    gene = geneSymbols,
    OrgDb = orgDb,
    keyType = "SYMBOL",
    ont = ont
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(Description, p.adjust, geneID, BgRatio,GeneRatio) %>%
    dplyr::mutate(
      BgRatio = {
        parts <- stringr::str_split(BgRatio, "/", simplify = TRUE)
        as.numeric(parts[, 1]) / as.numeric(parts[, 2])
      }
    ) %>%
    dplyr::mutate(
      GeneRatio = {
        parts <- stringr::str_split(GeneRatio, "/", simplify = TRUE)
        as.numeric(parts[, 1]) / as.numeric(parts[, 2])
      }
    ) %>%
    dplyr::filter(p.adjust <= 0.05)
}


weighted_genes_pathways %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::slice_head(n = 25) %>%
  ggplot()+
  geom_point(mapping = aes(x = GeneRatio,y = Description,
                           colour = p.adjust,size = count))+
  labs(color = "p-adjust",size = "number of genes")+
  xlab(label = "pathway ratio")+
  scale_color_gradient(low = "orange",high = "red")+
  ggtitle("PC1 gene pathways of top 5%",
          subtitle = "")+
  theme(axis.text = element_text(face = "bold",
                                 size =10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",
                                  size =15))

















#DEG table
res <- results(dds)
res_df <- as.data.frame(res) 
res_df$gene_name <- raw_filtered$external_gene_name
res_df$gene_id<- raw_filtered$Gene_ID
top <- res_df %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10)
library(ggrepel)
ggplot()+
  geom_point(data = res_df,
             mapping=aes(x=log2FoldChange,
                         y=-log10(padj)))+
  geom_text_repel(data = top,
                  mapping=aes(x=log2FoldChange,
                              y=-log10(padj),
                              label = gene_name))+
  ggtitle("")+
  labs(subtitle = "")+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 15,face = "bold"),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),  
  )+
  xlab(label = "Log2-FoldChange")+
  ylab(label = "p-adj")


res_counts_df <- res_df %>%
  dplyr::left_join(y = raw_filtered,
                   by = c("gene_id"="Gene_ID"))
write_csv(res_counts_df,"d1_L1_VSWT_maleVFemale_kidneyVLiverDeg.csv")

sig <- res_counts_df %>%
  dplyr::filter(padj <= 0.05,
                abs(log2FoldChange)>0.58)

sig_activated <- sig %>%
  dplyr::filter(log2FoldChange > 0)
sig_repressed <- sig %>%
  dplyr::filter(log2FoldChange < 0)

go_activated <- enrichGOHelperFunction(geneSymbols = sig_activated$gene_name,
                                       orgDb = org.Mm.eg.db,ont = "BP")
go_repressed <- enrichGOHelperFunction(geneSymbols = sig_repressed$gene_name,
                                       orgDb = org.Mm.eg.db,ont = "BP")


pathway_of_interest <- go_repressed %>%
  dplyr::filter(Description == "cilium organization") %>%
  dplyr::select(geneID) %>%
  stringr::str_split(pattern = "/")%>%
  unlist()


cilium_org_degs <- sig_repressed %>%
  dplyr::filter(gene_name %in% pathway_of_interest)
p_degs <- cilium_org_degs %>%
  dplyr::filter(padj <= 0.05,abs(log2FoldChange)>0.58)
