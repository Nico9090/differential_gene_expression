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

extract_go<-function(deg_table,org_db){
  go_results <- enrichGO(gene         = deg_table$gene_symbol,
                         OrgDb        = org_db,
                         keyType      = "SYMBOL",
                         ont          = "ALL",     # "MF", "CC" also valid
                         pAdjustMethod = "BH",pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2)
  return(go_results)
}

library(org.Mm.eg.db)
deg_list <- list(
  del_3_67 = mouse_del_3_67,
  del_3_4 = mouse_del_3_4,
  olson_pkhd1 = mouse_olson_pkhd1
)

go_results_list <- map(deg_list, ~ extract_go(.x, org.Mm.eg.db))
