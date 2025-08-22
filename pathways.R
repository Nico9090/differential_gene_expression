#Gene ontologies: BP, MF, CC
if(!"tidyverse" %in% loadedNamespaces()){
  library(tidyverse)
}
goBP<-readr::read_csv("pulmonary_fiborosis/go_enrichment_bp.csv")
goCC<-readr::read_csv("pulmonary_fiborosis/go_enrichment_cc.csv")
goMF<-readr::read_csv("pulmonary_fiborosis/go_enrichment_mf.csv")

#function
goTERM_plot <- function(goTERM_table, head, tail) {
  terms_toPLOT <- goTERM_table %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice(head:tail)
  
  ggplot(terms_toPLOT, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
    geom_col() +
    xlab("Count") +
    ylab("Ontology") +
    coord_flip() +
    theme_minimal()
}
SYMBOL_fromgoTERM<-function(goTERM_table, goID){
  geneIDS<-goTERM_table %>% 
    dplyr::filter(ID == goID) %>%
    dplyr::select(geneID)
  
  geneIDS<-unlist(stringr::str_split(geneIDS,"/"))
  if(!"clusterProfiler" %in% loadedNamespaces()){
  library(clusterProfiler)}
  gene_names<-clusterProfiler::bitr(geneID = geneIDS,
                                    fromType = "ENTREZID",
                                    toType = "SYMBOL",
                                    OrgDb = org.Hs.eg.db)
  return(gene_names$SYMBOL)
  
}

clusterGO<-function(goTERM_table){
  goTERM_table <- goTERM_table %>%
  mutate(gene_list = strsplit(geneID, "/"))
  all_genes <- unique(unlist(goTERM_table$gene_list))
  binary_matrix <- map(goTERM_table$gene_list, ~ all_genes %in% .x) %>%
  do.call(rbind, .)
  rownames(binary_matrix) <- goTERM_table$ID
  colnames(binary_matrix) <- all_genes
  
  library(umap)
  
  umap_result <- umap(binary_matrix)
  embedding <- as.data.frame(umap_result$layout)
  embedding$GO_ID <- rownames(binary_matrix)
  embedding <- embedding %>%
  left_join(goTERM_table %>% dplyr::select(ID, Description), by = c("GO_ID" = "ID"))
  
  
  plot<-plot_ly(data = embedding,x = ~V1,y = ~V2,type = 'scatter',
                mode = 'markers',  # Only markers, no static text
                hoverinfo = 'text',
                hovertext = ~paste("Description:", Description),
                marker = list(
                  symbol = 'x',
                  color = 'steelblue',
                  size = 10))
  return(list(plot=plot,clusterDATA=embedding))

}

#Visualize enrichment
if(!"clusterProfiler" %in% loadedNamespaces()){
  library(clusterProfiler)
}
goTERM_plot(goTERM_table = goBP,head = 1,tail = 10)
goTERM_plot(goTERM_table = goCC,head = 1,tail = 10)
goTERM_plot(goTERM_table = goMF,head = 1,tail = 10)

umap_ofCC<-clusterGO(goTERM_table = goCC)
umap_ofCC$plot

goCC %>% dplyr::filter(Description=="plasma membrane bounded cell projection cytoplasm") %>% dplyr::select(ID)

#Select a gene set from the go terms
library(org.Hs.eg.db)
collagen_genes<-SYMBOL_fromgoTERM(goTERM_table = goCC,goID = "GO:0062023") 
#for collagen-containing extracellular matrix
collagen_genes 

axoneme<-SYMBOL_fromgoTERM(goTERM_table = goCC,goID = "GO:0005930") 
#for collagen-containing extracellular matrix
axoneme 

actin_filament<-SYMBOL_fromgoTERM(goTERM_table = goCC,goID = "GO:0005884") 
#actin filament proteins are interesting, one includes a semaphorin regulating 
#(DPYSL3)
#PDZ domain also interesting
#TROPONIN complex
#IFT57 is also interesting
actin_filament 

synaptic_ves<-SYMBOL_fromgoTERM(goTERM_table = goCC,goID = "GO:0030672") 
synaptic_ves


pls_memb_proj<-SYMBOL_fromgoTERM(goTERM_table = goCC,goID = "GO:0032838")
pls_memb_proj
