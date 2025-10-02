title: "Assay Dev"
output: html_document
date: "2025-09-30"
---
```{r libraries}
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(STRINGdb)
library(visNetwork)
library(biomaRt)
```

```{r DEGs}
del34LiverDeg <- read_csv("csvFiles/D1_WT_malelivDEG.csv")
sig <- del34LiverDeg %>%
  dplyr::filter(padj <= 0.05)

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
del34LiverGo <- enrichGOHelperFunction(geneSymbols = sig$gene_name,
                                       orgDb = org.Mm.eg.db,ont = "BP")
```

```{r pathway focus}
pathway <- del34LiverGo %>%
  dplyr::filter(Description == "regulation of epithelial cell proliferation") %>%
  dplyr::select(geneID) %>%
  str_split("/")%>%
  unlist

pathwayGenes <- del34LiverDeg %>%
  dplyr::filter(gene_name %in% pathway)

pathwayGenes_CC <- enrichGOHelperFunction(geneSymbols = pathwayGenes$gene_name,
                                            orgDb = org.Mm.eg.db,
                                            ont = "CC")

pathwayGenes_MF <- enrichGOHelperFunction(geneSymbols = pathwayGenes$gene_name,
                                            orgDb = org.Mm.eg.db,
                                            ont = "MF")
group_annotations <- tribble(
  ~cell_component_group, ~genes,
  "beta-catenin binding", pathwayGenes_MF %>%
    dplyr::filter(Description == "beta-catenin binding") %>%
    dplyr::select(geneID)%>%
    unlist,
  "lipid transmembrane transporter activity", pathwayGenes_MF %>%
    dplyr::filter(Description == "lipid transmembrane transporter activity") %>%
    dplyr::select(geneID)%>%
    unlist,
  "G protein-coupled receptor binding", pathwayGenes_MF %>%
    dplyr::filter(Description == "G protein-coupled receptor binding") %>%
    dplyr::select(geneID)%>%
    unlist,
  "Notch binding",pathwayGenes_MF %>%
    dplyr::filter(Description == "Notch binding") %>%
    dplyr::select(geneID)%>%
    unlist,
  "mRNA binding",pathwayGenes_MF %>%
    dplyr::filter(Description == "mRNA binding") %>%
    dplyr::select(geneID)%>%
    unlist
)


gene_group_map <- group_annotations %>%
  tidyr::separate_rows(genes, sep = "/") %>%
  dplyr::rename(gene = genes)

pathwayGenes_annotated <- pathwayGenes %>%
  dplyr::left_join(gene_group_map, by = c("gene_name" = "gene")) %>%
  mutate(cell_component_group = coalesce(cell_component_group,
                                         "other"))


toPlot <- pathwayGenes_annotated %>%
  dplyr::filter(cell_component_group != "other")






ggplot(toPlot, aes(x = reorder(gene_name, log2FoldChange),
                   y = log2FoldChange,
                   fill = log2FoldChange)) +
  geom_col() +
  facet_wrap(~ cell_component_group, scales = "free_y") +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red4", midpoint = 0) +
  labs(title = "regulation of epithelial cell proliferation DEGs",
       x = "Gene Symbol",
       y = expression(log[2]("Fold Change")),
       fill = "Log2FC") +
  theme_minimal(base_size = 10) +
  theme(
  legend.position = "none",
  strip.text = element_text(face = "bold"),
  axis.text.y = element_text(size = 10),   # keep this one
  panel.spacing = unit(1, "lines"),
  axis.title.x = element_text(margin = margin(t = 10)),
  axis.text.x = element_text(size = 9),
  strip.background = element_rect(fill = "grey90", color = NA)
)
```
