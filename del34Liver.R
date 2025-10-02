source("~/omicsDataSets/rnaSeq/cnric/scripts/r_scripts/functions.R")
source("../../../singleCell/embryo/helperFunctions.R")
library(org.Mm.eg.db)
del34Liver <- read_csv("../../cnric/deg/D1_WT_malelivDEG.csv")
sigDel34Liver <- del34Liver %>%
  dplyr::filter(padj <= 0.05)

del34LiverGo <- enrichGOHelperFunction(geneSymbols = sigDel34Liver$gene_name,
                                       orgDb = org.Mm.eg.db,ont = "BP")

topGo <- del34LiverGo %>%
  dplyr::filter(Description %in% c("mesenchyme development",
                                   "cell-substrate adhesion",
                                   "epithelial cell development",
                                   "regulation of epithelial cell proliferation",
                                   "wound healing",
                                   "connective tissue development"
                                   ))

topGo <- topGo %>%
  dplyr::mutate(counts = geneID %>%
                  str_split("/") %>%
                  unlist %>%
                  length)
ggplot(topGo, aes(x = reorder(Description, ),
                   y = geneID %>%
                    str_split("/") %>%
                    unlist %>%
                    length,
                   fill = geneID %>%
                    str_split("/") %>%
                    unlist %>%
                    length)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red4", midpoint = 0) 











pathwayGenes <- del34LiverGo %>%
  dplyr::filter(Description == "mesenchyme development") %>%
  dplyr::select(geneID) %>%
  str_split("/")%>%
  unlist

pathwayDeg <- sigDel34Liver %>%
  dplyr::filter(gene_name %in% pathwayGenes)

pathwayGo_CC <- enrichGOHelperFunction(geneSymbols = pathwayDeg$gene_name,
                                            orgDb = org.Mm.eg.db,
                                            ont = "CC")

pathwayGo_MF <- enrichGOHelperFunction(geneSymbols = pathwayDeg$gene_name,
                                       orgDb = org.Mm.eg.db,
                                       ont = "MF")
group_annotations <- tribble(
  ~cell_component_group, ~genes,
  "SMAD protein complex", pathwayGo_CC %>%
    dplyr::filter(Description == "SMAD protein complex") %>%
    dplyr::select(geneID)%>%unlist,
  "cell leading edge", pathwayGo_CC %>%
    dplyr::filter(Description == "cell leading edge") %>%
    dplyr::select(geneID)%>%unlist,
  "beta-catenin destruction complex", pathwayGo_CC %>%
    dplyr::filter(Description == "beta-catenin destruction complex") %>%
    dplyr::select(geneID)%>%unlist,
  "retinoid binding", pathwayGo_MF %>%
    dplyr::filter(Description == "retinoid binding") %>%
    dplyr::select(geneID)%>%unlist,
  "site of polarized growth", pathwayGo_CC %>%
    dplyr::filter(Description == "site of polarized growth") %>%
    dplyr::select(geneID)%>%unlist,
  "collagen-containing extracellular matrix", pathwayGo_CC %>%
    dplyr::filter(Description == "collagen-containing extracellular matrix") %>%
    dplyr::select(geneID)%>%unlist,
  "frizzled binding",pathwayGo_MF %>%
    dplyr::filter(Description == "frizzled binding") %>%
    dplyr::select(geneID)%>%unlist,
  "protein serine/threonine kinase activity",pathwayGo_MF %>%
    dplyr::filter(Description == "protein serine/threonine kinase activity") %>%
    dplyr::select(geneID)%>%unlist,
  "semaphorin receptor binding",pathwayGo_MF %>%
    dplyr::filter(Description == "semaphorin receptor binding") %>%
    dplyr::select(geneID)%>%unlist,
  "peptidase inhibitor activity",pathwayGo_MF %>%
    dplyr::filter(Description == "peptidase inhibitor activity") %>%
    dplyr::select(geneID)%>%unlist
)

gene_group_map <- group_annotations %>%
  separate_rows(genes, sep = "/") %>%
  dplyr::rename(gene = genes)

pathwayDeg_annotated <- pathwayDeg %>%
  dplyr::left_join(gene_group_map, by = c("gene_name" = "gene")) %>%
  mutate(cell_component_group = coalesce(cell_component_group,
                                         "other"))


toPlot <- pathwayDeg_annotated %>%
  dplyr::filter(cell_component_group != "other")








ggplot(toPlot, aes(x = reorder(gene_name, log2FoldChange),
                   y = log2FoldChange,
                   fill = log2FoldChange)) +
  geom_col() +
  facet_wrap(~ cell_component_group, scales = "free_y") +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red4", midpoint = 0) +
  labs(title = "mesenchyme development DEGs",
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


forString <- del34LiverGo %>%
  dplyr::filter(Description == "cell-substrate adhesion")%>%
  dplyr::select(geneID) %>%
  str_split("/") %>%
  unlist

wntString <- stringDBHelperFunction(geneSymbols = forString)
plotVisNet(mapping = wntString$mapping,
           interactions = wntString$interactions)

target_gene <- "Fam107a"

target_ids <- wntString$mapping %>%
  filter(external_gene_name == target_gene) %>%
  pull(ensembl_peptide_id)

filtered_interactions <- wntString$interactions %>%
  mutate(from_clean = str_split(from, "\\.", simplify = TRUE)[,2]) %>%
  filter(from_clean %in% target_ids | to %in% target_ids)
proteins_in_network <- unique(c(filtered_interactions$from, filtered_interactions$to))

filtered_mapping <- wntString$mapping %>%
  dplyr::filter(ensembl_peptide_id %in% str_split(proteins_in_network,
                                           "\\.",simplify = TRUE)[,2])
plotVisNet(mapping = filtered_mapping,
           interactions = filtered_interactions)



library(ReactomePA)
pathwayGenes <- del34LiverGo %>%
  dplyr::filter(Description == "regulation of epithelial cell proliferation") %>%
  dplyr::select(geneID) %>%
  str_split("/")%>%
  unlist

pathwayDeg <- sigDel34Liver %>%
  dplyr::filter(gene_name %in% pathwayGenes)
# Convert to Entrez IDs
entrezIds <- bitr(sigDel34Liver$gene_name,
                  fromType="SYMBOL",
                  toType="ENTREZID", OrgDb="org.Mm.eg.db")

reactome_results <- enrichPathway(gene = entrezIds$ENTREZID,
                                  organism = "mouse",
                                  readable = TRUE)

barplot(reactome_results)
reactomeRes <- as.data.frame(reactome_results)
pathway_ids <- reactomeRes$ID
paste0("https://reactome.org/PathwayBrowser/#/R-HSA-", pathway_ids)












goTermDescriptions <- unique(del34LiverGo$Description)
all_terms <- as.list(GOTERM)
all_go_ids <- names(all_terms)
all_go_names <- sapply(all_terms, Term)

go_lookup <- setNames(all_go_ids, all_go_names)
matched_ids <- go_lookup[goTermDescriptions]
matched_ids <- matched_ids[!is.na(matched_ids)]
length(matched_ids) 

getAllParents <- function(go_ids, levels = 2, ontology = "BP") {
  if (ontology == "BP") {
    parent_db <- GOBPPARENTS
  } else if (ontology == "CC") {
    parent_db <- GOCCPARENTS
  } else if (ontology == "MF") {
    parent_db <- GOMFPARENTS
  } else {
    stop("Invalid ontology")
  }
  
  current <- go_ids
  for (i in 1:levels) {
    parents <- AnnotationDbi::mget(current, parent_db, ifnotfound = NA)
    current <- unique(unlist(parents))
    current <- current[!is.na(current)]
  }
  return(current)
}
high_level_ids <- getAllParents(matched_ids, levels = 5)
high_level_terms <- as.list(GOTERM[high_level_ids])
high_level_descriptions <- sapply(high_level_terms, Term)
high_level_descriptions





all_terms <- as.list(GOTERM)
go_id_to_desc <- sapply(all_terms, Term)  # GO ID -> Description
desc_to_go_id <- setNames(names(go_id_to_desc), go_id_to_desc)  # Description -> GO ID

# Match descriptions in your table to GO IDs
testGo <- del34LiverGo
testGo$go_id <- desc_to_go_id[testGo$Description]

# Remove rows that didn't match
testGo <- testGo[!is.na(testGo$go_id), ]
parents_list <- AnnotationDbi::mget(testGo$go_id,
                                    GOBPPARENTS, ifnotfound = NA)

# Extract first parent (some terms have multiple; you can handle differently if needed)
first_parents <- sapply(parents_list, function(x) if (!is.null(x)) x[3] else NA)

# Add parent GO ID column
testGo$parent_go_id <- first_parents

# Remove NAs (terms without parents)
noParentsGo <- testGo %>%
  dplyr::filter(is.na(parent_go_id))
testGo <- testGo[!is.na(testGo$parent_go_id), ]


parent_descs <- sapply(GOTERM[testGo$parent_go_id], Term)
testGo$parent_description <- parent_descs

final_df <- testGo[, c("parent_description", "geneID")]



head(final_df)
final_df <- final_df %>%
  distinct(parent_description,.keep_all = TRUE)



library(purrr)

final_df <- final_df %>%
  mutate(count = map_int(str_split(geneID, "/"), length))

ggplot(final_df, aes(x = reorder(parent_description,count),
                  y = count,
                  fill = count)) +
  geom_point() +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "purple4", high = "red4", midpoint = 10)+
  xlab(label = "Gene Ontology Term")+
  ylab(label = "Number of genes")+
  ggtitle(label = "GO Term Summarized of Pkhd1 del 3-4/del 3-4 Liver")



del34LiverGo <- del34LiverGo %>%
  mutate(count = map_int(str_split(geneID, "/"), length))



sorted <- del34LiverGo %>%
  dplyr::filter(Description %in% c("steroid metabolic process",
                                   "lipid transport",
                                   "organic acid biosynthetic process",
                                   "cell-matrix adhesion",
                                   "wound healing",
                                   "regulation of epithelial cell proliferation",
                                   "taxis"))
ggplot(sorted, aes(x = reorder(Description,p.adjust),
                     y = count,
                     fill = p.adjust,
                   size= count)) +
  geom_point() +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "purple4", high = "red4", midpoint = 3.125978e-06)+
  xlab(label = "Gene Ontology Term")+
  ylab(label = "Number of genes")+
  ggtitle(label = "GO Term Summarized of Pkhd1 del 3-4/del 3-4 Liver")
