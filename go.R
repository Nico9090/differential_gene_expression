#differentially expressed genes
gene_info <- bitr(for_clusterof16_cells_with_fpcL$gene_symbol,
                  fromType="SYMBOL",
                  toType="ENTREZID", 
                  OrgDb="org.Mm.eg.db")
go_results <- enrichGO(gene         = gene_info$ENTREZID,
                       OrgDb        = org.Mm.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "BP",     # "MF", "CC" also valid
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)
barplot(go_results, showCategory=16)

write.csv(as.data.frame(go_results),
          "GO_enrichment_chronic_stress.csv",
          row.names = FALSE)
go_results_df<-readr::read_csv("GO_enrichment_chronic_stress.csv")
