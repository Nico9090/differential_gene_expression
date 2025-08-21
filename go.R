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

packages<-c("tidyverse","limma","clusterProfiler",
            "plotly","DOSE","enrichplot","maEndToEnd","DESeq2",
            "EnhancedVolcano")
lapply(packages,base::library,character.only=TRUE)

#kidney differentially expressed genes
#load data
mouse_del_3_67<-readr::read_csv("") #data set 1
mouse_del_3_67<-mouse_del_3_67 %>% dplyr::filter(cluster==0)
#data set 2
mouse_del_3_4<-readr::read_csv("")

mouse_hybrid_pkd_pkhd1<-readr::read_delim("") 

mouse_olson_pkhd1<-mouse_hybrid_pkd_pkhd1[2:nrow(mouse_hybrid_pkd_pkhd1),
                                          c(1,2,9,10)] #data set 3, day 0
ols_p3<-mouse_hybrid_pkd_pkhd1[2:nrow(mouse_hybrid_pkd_pkhd1),
                               c(1,2,18,19)] #day 3
ols_p12<-mouse_hybrid_pkd_pkhd1[2:nrow(mouse_hybrid_pkd_pkhd1),
                               c(1,2,27,28)] #day 12
colnames(ols_p12)<-NULL
colnames(ols_p12) <- as.character(unlist(ols_p12[1, ]))
colnames(ols_p3)<-NULL
colnames(ols_p3)<-as.character(unlist(ols_p3[1, ]))
human_fpc_lof<-reahuman_fpc_lof<-reahuman_fpc_lof<-readr::read_csv("")
```

```{r}
#functions
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
library(org.Hs.eg.db)
library(enrichplot)
deg_list <- list(
  del_3_67 = mouse_del_3_67%>%dplyr::filter(p_val_adj<=0.05),
  del_3_4 = mouse_del_3_4%>%dplyr::filter(padj<=0.05),
  olson_pkhd1 = mouse_olson_pkhd1%>%dplyr::filter(`P-Value`<=0.05)
)

go_results_list <- map(deg_list, ~ extract_go(.x, org.Mm.eg.db))
del367_GO<-go_results_list$del_3_67
del34_GO<-go_results_list$del_3_4
ols_GO<-go_results_list$olson_pkhd1

human_GO<-extract_go(deg_table = human_fpc_lof%>%dplyr::filter(padj<=0.05),
                     org_db = org.Hs.eg.db)
ols_GO_p3<-extract_go(deg_table = ols_p3 %>% dplyr::filter(`P-Value`<=0.05),
                      org_db = org.Mm.eg.db)
ols_GO_p12<-extract_go(deg_table = ols_p12 %>% dplyr::filter(`P-Value`<=0.05),
                      org_db = org.Mm.eg.db)

barplot(del367_GO, showCategory = 10)+labs(caption="")
barplot(del34_GO, showCategory = 10)+labs(caption="")
barplot(ols_GO, showCategory = 10)+labs(caption="")
barplot(human_GO, showCategory = 20)+
  theme(axis.text.y = element_text(size = 8)) +labs(caption="")

barplot(ols_GO_p3, showCategory = 10)+labs(caption="")

barplot(ols_GO_p12, showCategory = 10)+labs(caption="")

human_GO <- pairwise_termsim(human_GO)
treeplot(human_GO, showCategory = 20, color = "p.adjust", label_format = 40)

















