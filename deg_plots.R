#Plots
clusterProfiler::dotplot(deg_go,
                         showCategory=5,
                         split=".sign") + 
  facet_grid(.~.sign) + 
  labs(caption = "Enriched pathways in human ARPKD")


plot_ly(
  data = growth_factor_act,
  x = ~log2FoldChange,
  y = ~-log10(padj),
  text = ~gene_name,
  type = 'scatter',
  mode = 'markers',
  marker = list(color = 'rgba(255, 0, 0, .5)', size = 10)
)
#<--- see genes in up vs down regulated
ggplot2::ggplot(deg_table[deg_table$gene_name %in%
                            deg_go@geneSets[["GO:0002181"]],], 
                aes(x = reorder(gene_name, log2FoldChange),
                    y = log2FoldChange,
                    fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  scale_fill_manual(values = c("red", "blue"), labels = c("Suppressed", "Activated")) +
  labs(x = "protein",
       y = "log Fold Change",
       fill = "Expression",
       caption = "Significant differentially expressed:
      Cytoplasmic translation.
       ") +
  theme_minimal()

#Reactome 
entrez_ids <- bitr(hs_deg$gene.symbol,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)


reactome_enrich <- enrichPathway(gene = entrez_ids$ENTREZID,
                                 organism = "human",  # or "mouse", etc.
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)
barplot(reactome_enrich, showCategory = 5)


