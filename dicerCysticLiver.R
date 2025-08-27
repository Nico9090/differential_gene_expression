library(rtracklayer)
library(dplyr)
#For mouse
#Mut1
mut1Gtf <- import("GSM7286688_LiverMut1.gtf") 
mut1ExprTable <- as.data.frame(mcols(mut1Gtf))[,c("gene_id","transcript_id","FPKM","TPM")]
geneExprMut1 <- mut1ExprTable %>%
  group_by(gene_id) %>%
  summarize(meanTPM = mean(as.numeric(TPM), na.rm=TRUE),
            meanFPKM = mean(as.numeric(FPKM), na.rm=TRUE))
#WT1
wt1Gtf <- import("GSM7286690_LiverWT1.gtf")
wt1ExprTable <- as.data.frame(mcols(wt1Gtf))[,c("gene_id","transcript_id","FPKM","TPM")]
geneExprWT1 <- wt1ExprTable %>%
  dplyr::group_by(gene_id) %>%
  summarize(meanTPM = mean(as.numeric(TPM), na.rm=TRUE),
            meanFPKM = mean(as.numeric(FPKM), na.rm=TRUE))
#Mut2
mut2Gtf <- import("GSM7286689_LiverMut2.gtf")
mut2ExprTable <- as.data.frame(mcols(mut2Gtf))[,c("gene_id","transcript_id","FPKM","TPM")]
geneExprMut2 <- mut2ExprTable %>%
  dplyr::group_by(gene_id) %>%
  summarize(meanTPM = mean(as.numeric(TPM), na.rm=TRUE),
            meanFPKM = mean(as.numeric(FPKM), na.rm=TRUE))
#WT2
wt2Gtf <- import("GSM7286691_LiverWT2.gtf")
wt2ExprTable <- as.data.frame(mcols(wt2Gtf))[,c("gene_id","transcript_id","FPKM","TPM")]
geneExprWT2 <- wt2ExprTable %>%
  dplyr::group_by(gene_id) %>%
  summarize(meanTPM = mean(as.numeric(TPM), na.rm=TRUE),
            meanFPKM = mean(as.numeric(FPKM), na.rm=TRUE))
#Join based on TPM
wt1 <- geneExprWT1 %>% select(gene_id, meanTPM) %>% rename(WT1 = meanTPM)
wt2 <- geneExprWT2 %>% select(gene_id, meanTPM) %>% rename(WT2 = meanTPM)
mut1 <- geneExprMut1 %>% select(gene_id, meanTPM) %>% rename(Mut1 = meanTPM)
mut2 <- geneExprMut2 %>% select(gene_id, meanTPM) %>% rename(Mut2 = meanTPM)

exprTable <- wt1 %>%
  inner_join(wt2, by="gene_id") %>%
  inner_join(mut1, by="gene_id") %>%
  inner_join(mut2, by="gene_id")
#<---------------------------------------------------------------->

limmaFit<-function(counts,
                   design,
                   meta_data,
                   sample_identifier_by){
  #<---sample_identifier_by= condition or column 2 of meta data
  colnames(design)<-levels(meta_data[[sample_identifier_by]])
  fit<-limma::lmFit(as.matrix(log2(counts+1)),design)
  #<--change PLDF - WTF when needed --->
  contrast_matrix <- limma::makeContrasts(Mut - WT,
                                          levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  deg_table <- limma::topTable(fit2, adjust = "fdr", number = Inf)
  return(deg_table)
}
onlyTPM <- exprTable[,-1]
metaData <- dplyr::tibble(samples=names(onlyTPM),
                          condition=c(rep("WT",2),
                                      rep("Mut",2)))
metaData$condition <- factor(x = metaData$condition,
                             levels = c("WT","Mut"))
design <- model.matrix(~0 + condition,
                       data = metaData)
names(design) <- c("WT","Mut")

degTable <- limmaFit(counts = onlyTPM,
                     design = design,
                     meta_data = metaData,
                     sample_identifier_by = "condition")

degTable$geneId <- exprTable$gene_id

#obtain gene names
library(clusterProfiler)
library(org.Mm.eg.db)
geneInformation <- clusterProfiler::bitr(geneID = degTable$geneId,
                                         fromType = "ENTREZID",
                                         toType = "SYMBOL",
                                         OrgDb = org.Mm.eg.db)
geneInformation <- geneInformation %>% 
  dplyr::filter(!is.na(SYMBOL))

degTable <- degTable[degTable$geneId %in% geneInformation$ENTREZID,
                     ]
degTable$geneSymbol <- geneInformation$SYMBOL
library(tidyverse)
write_csv(degTable,"dicerCysticLiverMouseDeg.csv")
degTable <- readr::read_csv("dicerCysticLiverMouseDeg.csv")
degTable<- degTable %>%
  mutate(color = case_when(
    
    adj.P.Val <= 0.05 & logFC >= 1.5  ~ "blue4",
    adj.P.Val <= 0.05 & logFC <= -1.5 ~ "red4",
    AveExpr >=5 ~"green4",
    TRUE                                  ~ "black" 
  ))
degTable <- degTable %>%
  mutate(group = case_when(
    
    adj.P.Val <= 0.05 & logFC >= 1.5  ~ "Up-regulated",
    adj.P.Val <= 0.05 & logFC <= -1.5 ~ "Down-regulated",
    AveExpr >=5 ~"High expression",
    
    TRUE                              ~ "Not significantly expressed"
  ))


ggplot2::ggplot(data = degTable,
                mapping = aes(x = logFC,
                              y = -log10(adj.P.Val),
                              color = group)) +
  xlim(-15,15)+
  geom_point(shape=3,size=2.5,stroke = 1.5)+ 
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = "dashed",  # or "dashed"
             color = "black",
             linewidth=0.8)+      # color of the line
  geom_hline(yintercept = 1.3,      # e.g., -log10(0.05) ≈ 1.3
             linetype = "dashed",
             color = "black",
             linewidth=0.8)+
  scale_color_manual(values = c("Up-regulated" = "blue4",
                                "Down-regulated" = "red4",
                                "High expression" = "green4",
                                "Not significantly expressed"="black"))+
  labs(color="Expression",title="Cystic Liver Mouse from Dicer Syndrome")


rankedGenes <- degTable$logFC
names(rankedGenes) <- degTable$geneSymbol
rankedGenes <- sort(rankedGenes,
                    decreasing = TRUE)
goSeqFunction <- function(entrezIds, degTable) {
  # Load libraries
  library(biomaRt)
  library(goseq)
  library(org.Mm.eg.db)
  
  # Connect to Ensembl BioMart
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Fetch gene lengths for Entrez IDs
  lengths <- getBM(
    attributes = c("entrezgene_id", "start_position", "end_position"),
    filters = "entrezgene_id",
    values = entrezIds,
    mart = ensembl
  )
  lengths$gene_length <- lengths$end_position - lengths$start_position + 1
  
  # Create named vector of gene lengths
  gene_lengths <- lengths$gene_length
  names(gene_lengths) <- as.character(lengths$entrezgene_id)
  
  # Create binary DE vector: 1 = significant, 0 = not
  sigData <- as.integer(!is.na(degTable$adj.P.Val) & degTable$adj.P.Val <= 0.05)
  names(sigData) <- as.character(degTable$geneId)  # must match Entrez IDs
  
  # Keep only genes present in both DE vector and lengths
  keep <- intersect(names(sigData), names(gene_lengths))
  sigData <- sigData[keep]
  gene_lengths <- gene_lengths[keep]
  
  # Remove any NAs
  good <- !is.na(sigData) & !is.na(gene_lengths)
  sigData <- sigData[good]
  gene_lengths <- gene_lengths[good]
  
  # Probability weighting function (length bias correction)
  pwf <- nullp(sigData, bias.data = gene_lengths, genome = "mm10", id = "entrez")
  
  # Fetch GO annotations manually for mouse
  GOcats <- as.list(org.Mm.egGO2ALLEGS)
  
  # Run GO enrichment (Biological Process)
  goResults <- goseq(pwf, gene2cat = GOcats, test.cats = "GO:BP")
  
  return(goResults)
}




goRes <- goSeqFunction(entrezIds = degTable$geneId,
                       degTable = degTable)
goRes %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")


GOcats = as.list(org.Mm.egGO2ALLEGS)

sigData <- as.integer(!is.na(degTable$adj.P.Val) & degTable$adj.P.Val <= 0.05)
names(sigData) <- as.character(degTable$geneId)


goRes$genes <- lapply(goRes$category, function(go_term) {
  intersect(names(sigData)[sigData==1], GOcats[[go_term]])
})


goRes <- goRes %>%
  rowwise() %>%
  mutate(avgLog2FC = mean(degTable$logFC[degTable$geneId %in% genes], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(direction = ifelse(avgLog2FC > 0, "Activated", "Suppressed"))

goRes %>%
  top_n(30, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc = numDEInCat*100/numInCat) %>%
  ggplot(aes(x = hitsPerc, 
             y = term, 
             colour = direction,
             size = numDEInCat)) +
  geom_point() +
  scale_colour_manual(values = c("Activated" = "blue", "Suppressed" = "red")) +
  expand_limits(x=0) +
  labs(x = "Hits (%)", y = "GO term", colour = "Direction", size = "Count",
       title = "Dicer Syndrome Cystic Liver")

















del34Liv<-readr::read_csv("~/omics_DS/rna_seq/cnric/deg/D1_WT_malelivDEG.csv")
del34GeneInformation <- bitr(geneID = del34Liv$gene_name,
                             fromType = "SYMBOL",
                             toType = "ENTREZID",
                             OrgDb = org.Mm.eg.db)
del34GeneInformation <- del34GeneInformation %>%
  dplyr::distinct(SYMBOL,.keep_all = TRUE)
del34Liv <- del34Liv %>% 
  dplyr:: distinct(gene_name,.keep_all = TRUE)
del34Liv <- del34Liv[del34Liv$gene_name %in%
                       del34GeneInformation$SYMBOL,]
del34Liv$entrezID <- del34GeneInformation$ENTREZID
del34Liv <- dplyr::rename(del34Liv,adj.P.Val = padj )
del34Liv <- dplyr::rename(del34Liv,geneId = entrezID )




del34GoRes <- goSeqFunction(entrezIds = del34Liv$geneId,
                            degTable = del34Liv)
del34GoRes %>% 
  top_n(30, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

GOcats = as.list(org.Mm.egGO2ALLEGS)

sigData <- as.integer(!is.na(del34Liv$adj.P.Val) & del34Liv$adj.P.Val <= 0.05)
names(sigData) <- as.character(del34Liv$geneId)


del34GoRes$genes <- lapply(del34GoRes$category, function(go_term) {
  intersect(names(sigData)[sigData==1], GOcats[[go_term]])
})


del34GoRes <- del34GoRes %>%
  rowwise() %>%
  mutate(avgLog2FC = mean(del34Liv$log2FoldChange[del34Liv$geneId %in% genes], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(direction = ifelse(avgLog2FC > 0, "Activated", "Suppressed"))

del34GoRes %>%
  top_n(30, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc = numDEInCat*100/numInCat) %>%
  ggplot(aes(x = hitsPerc, 
             y = term, 
             colour = direction,
             size = numDEInCat)) +
  geom_point() +
  scale_colour_manual(values = c("Activated" = "blue", "Suppressed" = "red")) +
  expand_limits(x=0) +
  labs(x = "Hits (%)", y = "GO term", colour = "Direction", size = "Count",
       title = "Pkhd1 del 3-4/del 3-4 liver")

mapkCascadeLiver <- del34GoRes %>%
  filter(category == "GO:0000165")
mapkGenes <- mapkCascadeLiver %>%
  pull(genes)
mapkGenes<- unlist(mapkGenes)

mapkGeneInformation <- bitr(geneID = mapkGenes,
                            fromType = "ENTREZID",
                            toType = "SYMBOL",
                            OrgDb = org.Mm.eg.db)
del34LivMapKinase <- del34Liv[del34Liv$gene_name %in%
                                mapkGeneInformation$SYMBOL,]

del34LivMapKinase<- del34LivMapKinase %>%
  mutate(color = case_when(
    baseMean >=1000 ~"green4",
    adj.P.Val <= 0.05 & log2FoldChange >= 1.5  ~ "blue4",
    adj.P.Val <= 0.05 & log2FoldChange <= -1.5 ~ "red4",
    
    TRUE                                  ~ "black" 
  ))
del34LivMapKinase <- del34LivMapKinase %>%
  mutate(group = case_when(
    baseMean >=1000 ~"High counts",
    adj.P.Val <= 0.05 & log2FoldChange >= 1.5  ~ "Up-regulated",
    adj.P.Val <= 0.05 & log2FoldChange <= -1.5 ~ "Down-regulated",
    TRUE                              ~ "Not significantly expressed"
  ))
top30 <- del34LivMapKinase %>% 
  dplyr::arrange(adj.P.Val) %>%
  slice_head(n = 30)

library(ggrepel)
ggplot2::ggplot(data = del34LivMapKinase,
                mapping = aes(x = log2FoldChange,
                              y = -log10(adj.P.Val),
                              color = group)) +
  xlim(-15,15)+
  geom_point(shape=3,size=2.5,stroke = 1.5)+ 
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = "dashed",  # or "dashed"
             color = "black",
             linewidth=0.8)+      # color of the line
  geom_hline(yintercept = 1.3,      # e.g., -log10(0.05) ≈ 1.3
             linetype = "dashed",
             color = "black",
             linewidth=0.8)+
  geom_text_repel(data = top7,aes(label=gene_name))+
  scale_color_manual(values = c("Up-regulated" = "blue4",
                                "Down-regulated" = "red4",
                                "High expression" = "green4",
                                "Not significantly expressed"="black"))+
  labs(color="Expression",title="Pkhd1 del 3-4 Liver MAPK Cascade")

ggplot(top30, aes(x = reorder(gene_name,
                                          log2FoldChange),
                              y = log2FoldChange, fill = log2FoldChange)) +
  geom_bar(stat = "identity") +               # Use identity for precomputed values
  coord_flip() +                              # Flip axes if labels are long
  labs(
    x = "Gene",
    y = "LogFC",
    title = "Pkhd1 del 3-4 Liver MapK Cascade"
  ) +
  theme_minimal() +                           # Clean theme
  theme(
    axis.text.y = element_text(size = 10),    # Adjust label size
    legend.position = "bottom"
  )
