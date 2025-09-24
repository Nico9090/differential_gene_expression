library(tidyverse)
library(limma)
library(GEOquery)
library(enrichplot)
library(org.Hs.eg.db)
gpl <- getGEO("GPL5639", AnnotGPL = TRUE) %>%
  Table()
head(colnames(gpl))

loadData <- read_delim("GSM3583748_QH58429.txt",skip = 7) %>%
  dplyr::select(G_Name,`635nm...7`,`532nm...8`,`635nm...9`,`532nm...10`)

annot <- gpl %>%
  dplyr::select(ID, Gene_Symbol, Ensembl_Gene_ID) %>%
  dplyr::rename(G_Name = ID, GeneSymbol = Gene_Symbol)

workData <- loadData %>%
  left_join(annot, by = "G_Name")


#unsure <------
RG <- list(
  R = workData$`635nm...7`,   # Cy5 (red)
  G = workData$`532nm...8`,   # Cy3 (green)
  genes = workData[, c("G_Name", "GeneSymbol")]
) #for replicate 1

class(RG) <- "RGList" 
RG <- backgroundCorrect(RG, method = "normexp")
MA <- normalizeWithinArrays(RG, method = "loess")
MA <- normalizeBetweenArrays(MA, method = "Aquantile")
#--->
names(workData)
onlyIntensities <- workData[,-c(1,6,7)]
metaData <- tibble(
  samples = names(onlyIntensities),
  Condition = c("Het","KO","Het","KO")
)

metaData$Condition <- factor(metaData$Condition, levels = c("Het","KO"))


limma_fit_for_deg<-function(counts,
                            design,
                            meta_data,
                            sample_identifier_by){
  #<---sample_identifier_by= condition or column 2 of meta data
  colnames(design)<-levels(meta_data[[sample_identifier_by]])
  fit<-limma::lmFit(as.matrix(log2(counts+1)),design)
  #<--change PLDF - WTF when needed --->
  contrast_matrix <- limma::makeContrasts(KO - Het,
                                          levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  deg_table <- limma::topTable(fit2, adjust = "fdr", number = Inf)
  return(deg_table)
}

design <- model.matrix(~0 + Condition,
                       data = metaData) 
deGenes <- limma_fit_for_deg(counts = onlyIntensities,
                             design = design,
                             meta_data = metaData,sample_identifier_by = "Condition")
deGenes$geneSymbol <- workData$GeneSymbol

ggplot()+
  geom_point(data = deGenes,mapping = aes(x=logFC,
                                          y=-log10(adj.P.Val) ))

sigDeGenes <- deGenes %>%
  dplyr::filter(adj.P.Val <= 0.05)

sigGo <- enrichGO(
  gene          = sigDeGenes$geneSymbol,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
) %>%
  as.data.frame() %>%     
  as_tibble() %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::filter(p.adjust <= 0.05)
