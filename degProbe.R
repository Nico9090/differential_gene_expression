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
```{r classify genes}
#pull out non-coding RNAs
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
annot <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol"),
               filters = "ensembl_gene_id",
               values = sig$Gene_ID,
               mart = mart)

sig_deg_annot <- merge(sig, annot, by.x="Gene_ID",
                       by.y="ensembl_gene_id")







ncRNAs <- subset(sig_deg_annot, gene_biotype != "protein_coding")

longncRNAs <- subset(ncRNAs,gene_biotype == "lncRNA")

lncRNAGo <- enrichGOHelperFunction(geneSymbols = longncRNAs$gene_name,
                                   orgDb = org.Mm.eg.db,ont = "BP")
head(ncRNAs)
```


```{python}
import pandas as pd
rows = []
for name, ids in proteinIds.items():
    for uid in ids:
      try:
        info = getProteinDescriptions(uid)
        rows.append({
            "proteinName": name,
            "uniprotId": uid,
            "description": info
        })
      except (KeyError,IndexError):
        rows.append({
            "proteinName": name,
            "uniprotId": uid,
            "description": "N/A"
        })

proteinTable = pd.DataFrame(rows)
```

```{python}
print(proteinTable.head(20))
proteinTable.to_csv("cellsubstrAdhesionGenesInDel34Liver.csv")
```
```{r}
geneDescriptions <- read_csv("epithDevGenesInDel34Liver.csv")
```
```{r}
stringDBHelperFunction <- function(geneSymbols){
  #orgId for mouse == 100900
  stringDB <- STRINGdb$new(version="11.5", species=10090, score_threshold=400)
  mapped <- stringDB$map(data.frame(gene=geneSymbols),
                         "gene",
                         removeUnmappedRows=TRUE)
  
  interactions <- stringDB$get_interactions(mapped$STRING_id)
  proteinIds <- unique(c(interactions$from,
                         interactions$to))
  proteinIdsClean <- gsub("10090\\.", "",
                          proteinIds) #change if needed
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mapping <- getBM(
    attributes = c("ensembl_peptide_id", "external_gene_name"),
    filters = "ensembl_peptide_id",
    values = proteinIdsClean,
    mart = ensembl
  )
  return(list(mapping = mapping,
              interactions = interactions))
  
}

plotVisNet <- function(mapping, interactions){
  # Create node list
  proteinIds <- unique(c(interactions$from, interactions$to))
  
  nodes <- data.frame(id = proteinIds, stringsAsFactors = FALSE)
  
  nodes <- nodes %>%
    mutate(
      cleanId = str_replace(id, "10090\\.", ""),
      label = mapping$external_gene_name[match(cleanId, mapping$ensembl_peptide_id)],
      group = mapping$cell_component_group[match(cleanId, mapping$ensembl_peptide_id)]
    )
  
  edges <- data.frame(
    from = interactions$from,
    to = interactions$to,
    width = interactions$combined_score / 200
  )
  
  # visNetwork plot with nodes grouped by functional group
  visNetwork(nodes, edges) %>%
    visNodes() %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLegend() %>%
    visPhysics(enabled = FALSE)
}
```

```{r}
forString <- del34LiverGo %>%
  dplyr::filter(Description == "regulation of epithelial cell proliferation")%>%
  dplyr::select(geneID) %>%
  str_split("/") %>%
  unlist

string <- stringDBHelperFunction(geneSymbols = forString)
plotVisNet(mapping = string$mapping,
           interactions = string$interactions)








gene_group_map <- group_annotations %>%
  tidyr::separate_rows(genes, sep = "/") %>%
  dplyr::rename(gene = genes)

string$mapping$annotation <- string$mapping %>%
  dplyr::left_join(gene_group_map, by = c("external_gene_name" = "gene")) %>%
  mutate(cell_component_group = coalesce(cell_component_group,
                                         "other"))
string$mapping <- string$mapping %>%
  left_join(gene_group_map, by = c("external_gene_name" = "gene")) %>%
  mutate(cell_component_group = coalesce(cell_component_group, "other"))

string$mapping <- string$mapping %>%
  mutate(color = case_when(
    cell_component_group == "lipid transmembrane transporter activity" ~ "red4",
    cell_component_group == "mRNA binding" ~ "blue4",
    cell_component_group == "G protein-coupled receptor binding" ~ "green4",
    cell_component_group == "beta-catenin binding"  ~ "orange",
    cell_component_group == "Notch binding"  ~ "purple4",
    cell_component_group == "Notch binding"  ~ "grey"
    
  ))
toPlot <- string$mapping %>%
  dplyr::filter(cell_component_group != "other")


plotVisNet(mapping = string$mapping,
           interactions = string$interactions)
```
```{r}
target_gene <- "Vegfa"

target_ids <- string$mapping %>%
  filter(external_gene_name == target_gene) %>%
  pull(ensembl_peptide_id)

filtered_interactions <- string$interactions %>%
  mutate(from_clean = str_split(from, "\\.", simplify = TRUE)[,2]) %>%
  filter(from_clean %in% target_ids | to %in% target_ids)
proteins_in_network <- unique(c(filtered_interactions$from, filtered_interactions$to))

filtered_mapping <- string$mapping %>%
  filter(ensembl_peptide_id %in% str_split(proteins_in_network,
                                           "\\.",simplify = TRUE)[,2])
plotVisNet(mapping = filtered_mapping,
           interactions = filtered_interactions)
```

```{r}
see <- read_csv("../RegEpCellProlifDel34Liver.csv")
see <- see[,-1]
colnames(see) <- c("geneName","RefFromUniprot")
write_csv(see,"../RegEpCellProlifDel34Liver.csv")









see <- see %>%
  mutate(
    RefFromUniprot = str_remove_all(RefFromUniprot, "\\[|\\]|'"),        # remove brackets/quotes
    RefFromUniprot = str_split(RefFromUniprot, ",\\s*"),
    n_refs   = sapply(RefFromUniprot, length) 
  )

ggplot(see, aes(x = geneName, y = n_refs)) +
  geom_col(fill = "steelblue") +
  labs(x = "Gene", y = "Number of references") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
