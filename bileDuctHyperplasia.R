library(clusterProfiler)
library(org.Mm.eg.db)
library(GEOquery)
library(limma)
library(umap)
library(readxl)
library(DESeq2)
library(ggrepel)
library(STRINGdb)
library(biomaRt)
source("~/omicsDataSets/rnaSeq/cnric/scripts/r_scripts/functions.R")
#source : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280724
bileDuctObstruction <- read_excel("rawData/GSE280724_count_data_24hrs.xlsx",
                                  sheet = 2)
metaData <- tibble(samples = colnames(bileDuctObstruction[,-1]),
                   condition = c(rep("injured",5),
                                 rep("control",4)))
metaData$condition <- factor(metaData$condition,
                             levels= c("control","injured"))

deSeqRes <- DESeq2.fn(raw.counts = bileDuctObstruction[,-1],
                      meta.data = metaData,
                      design = "condition")

degTable <- deSeqRes$dds_res %>% as.data.frame
degTable$geneId <- bileDuctObstruction$ensembl_gene_id

#obtain gene names
geneInformation <- bitr(geneID = degTable$geneId,
                        fromType = "ENSEMBL",toType = "SYMBOL",
                        OrgDb = org.Mm.eg.db)

mappedDegTable <- degTable %>%
  dplyr::filter(geneId %in% geneInformation$ENSEMBL)

geneInformation <- geneInformation %>%
  dplyr::distinct(ENSEMBL,.keep_all = TRUE)

mappedDegTable$geneSymbol <- geneInformation$SYMBOL

#in bile duct obstruction cholangiocyte proliferation increases
#what are genes involved 
topGenes <- mappedDegTable %>%
  dplyr::arrange(padj) %>%
  slice_head(n = 10)

ggplot()+
  geom_point(data = mappedDegTable,mapping = aes(x = log2FoldChange,
                                                 y = -log10(padj) ))+
  geom_text_repel(data = topGenes,mapping = aes(x = log2FoldChange,
                                                y = -log10(padj),
                                                label = geneSymbol))

significant <- mappedDegTable %>%
  dplyr::filter(padj <= 0.05)

bileDuctObstruction_Go <- enrichGOHelperFunction(geneSymbols = significant$geneSymbol,
                                                 orgDb = org.Mm.eg.db,
                                                 ont = "BP")

wntSignalingGenes <- bileDuctObstruction_Go %>%
  dplyr::filter(Description == "Wnt signaling pathway") %>%
  dplyr::select(geneID) %>%
  str_split("/")%>%
  unlist

wntSignalingDeg <- mappedDegTable %>%
  dplyr::filter(geneSymbol %in% wntSignalingGenes)

wntSignalingGo_CC <- enrichGOHelperFunction(geneSymbols = wntSignalingDeg$geneSymbol,
                                            orgDb = org.Mm.eg.db,
                                            ont = "CC")


group_annotations <- tribble(
  ~cell_component_group, ~genes,
  "beta-catenin-TCF complex", "Ctnnb1/Tle4/Ldb1/Lef1/Tle3/Bcl9/Pygo2/Tcf7l1/Bcl9l",
  "beta-catenin destruction complex", "Axin2/Ctnnb1/Csnk1a1/Siah1a/Siah1b/Dact1/Gsk3a",
  "ATPase complex", "Ruvbl2/Vps4b/Hdac2/Atp6v1c2/Mbd2/Vcp/Hdac1/Ruvbl1/Brd",
  "euchromatin","Klf4/Ruvbl2/Ctr9/Ctnnb1/Myc/Smarca4/Skic8",
  "collagen-containing extracellular matrix","Col6a1/Nid1/Fgf10/Fgf9/Cela1/Myoc/Fgfr2/Gpc4/Sfrp1/Ndp",
  "ciliary basal body","Ift20/Csnk1a1/Csnk1d/Jade1/Pin1/Pkd2/Nphp4"
)


gene_group_map <- group_annotations %>%
  separate_rows(genes, sep = "/") %>%
  rename(gene = genes)

wntSignalingDeg_annotated <- wntSignalingDeg %>%
  dplyr::left_join(gene_group_map, by = c("geneSymbol" = "gene")) %>%
  mutate(cell_component_group = coalesce(cell_component_group,
                                         "other"))


toPlot <- wntSignalingDeg_annotated %>%
  dplyr::filter(cell_component_group != "other")









ggplot() +
  geom_bar(data = toPlot, 
           mapping = aes(x = reorder(geneSymbol, log2FoldChange), 
                         y = log2FoldChange, 
                         fill = log2FoldChange), 
           stat = "identity") +
  facet_wrap(~ cell_component_group, scales = "free_y") +
  theme_minimal(base_size = 10) +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red4", midpoint = 0) +
  labs(title = "Wnt signaling genes in liver hyperplasia",
       x = "Gene Symbol",
       y = expression(log[2]("Fold Change")),
       fill = "Log2FC") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 10),
        panel.spacing = unit(1, "lines"))



forString <- bileDuctObstruction_Go %>%
  dplyr::filter(Description == "Wnt signaling pathway")%>%
  dplyr::select(geneID) %>%
  str_split("/") %>%
  unlist

wntString <- stringDBHelperFunction(geneSymbols = forString)
plotVisNet(mapping = wntString$mapping,
           interactions = wntString$interactions)

target_gene <- "Siah1b"

target_ids <- wntString$mapping %>%
  filter(external_gene_name == target_gene) %>%
  pull(ensembl_peptide_id)

filtered_interactions <- wntString$interactions %>%
  mutate(from_clean = str_split(from, "\\.", simplify = TRUE)[,2]) %>%
  filter(from_clean %in% target_ids | to %in% target_ids)
proteins_in_network <- unique(c(filtered_interactions$from, filtered_interactions$to))

filtered_mapping <- wntString$mapping %>%
  filter(ensembl_peptide_id %in% str_split(proteins_in_network,
                                           "\\.",simplify = TRUE)[,2])
plotVisNet(mapping = filtered_mapping,
           interactions = filtered_interactions)
