go_table<-tibble(ontology=cluster0_0_go@result$ONTOLOGY,
                 id=cluster0_0_go@result$ID,
                 desc=cluster0_0_go@result$Description,
                 genes=cluster0_0_go@result$core_enrichment)

bp_go_terms<-go_table[go_table$ontology=="BP",]

bp_go_terms <- bp_go_terms %>%
  mutate(gene_list = strsplit(genes, "/"))
all_genes <- unique(unlist(bp_go_terms$gene_list))
binary_matrix <- map(bp_go_terms$gene_list, ~ all_genes %in% .x) %>%
  do.call(rbind, .)
rownames(binary_matrix) <- bp_go_terms$id
colnames(binary_matrix) <- all_genes
library(umap)
umap_result <- umap(binary_matrix)
embedding <- as.data.frame(umap_result$layout)
embedding$GO_ID <- rownames(binary_matrix)
embedding <- embedding %>%
  left_join(bp_go_terms %>% dplyr::select(id, desc), by = c("GO_ID" = "id"))

ggplot(embedding, aes(x = V1, y = V2)) +
  geom_point(shape = "x",color = "steelblue", size = 3) +
  geom_text(aes(label = desc), size = 3, hjust = 1.2, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "UMAP of GO BP Terms", x = "UMAP1", y = "UMAP2")


plot_ly(
  data = embedding,
  x = ~V1,
  y = ~V2,
  type = 'scatter',
  mode = 'markers',  # Only markers, no static text
  hoverinfo = 'text',
  hovertext = ~paste(
    "Description:", desc,
    "<br>GO ID:", GO_ID
    # Add more fields if needed (e.g. p-value, NES)
  ),
  marker = list(
    symbol = 'x',
    color = 'steelblue',
    size = 10
  )
) %>%
  layout(
    title = "UMAP of GO BP Terms",
    xaxis = list(title = "UMAP1"),
    yaxis = list(title = "UMAP2")
  )
