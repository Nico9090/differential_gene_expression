#polycystsic organoid dep plots
clusterProfiler::dotplot(interested_cyto_skel_gsea_go,
                         showCategory=10,
                         split=".sign") + 
  facet_grid(.~.sign) + 
  labs(caption = "Enriched cytoskeletal protein pathways in polycystic organoid")


cyt.motor.volcano<-EnhancedVolcano::EnhancedVolcano(toptable = cyt.motor,
                    lab =cyt.motor[["protein"]],
                    x ="logFC",y ="adj.P.Val",xlim =c(-4,4),ylim =c(0.0,8),
                    pCutoff =0.05,FCcutoff = 1,pointSize = 2.5,
                    labSize = 3.5,
                    caption ="Smooth muscle cell protein enrichment in polycystic disease kidney organoid",
                    legendLabels =c("NS","logFC","p-adj","p-adj & logFC"),
                    legendPosition = "right",legendLabSize = 12,legendIconSize = 4.0,
                    drawConnectors = FALSE)
#<--interactive volcano plot--->
plot_ly(
  data = cyt.motor,
  x = ~logFC,
  y = ~-log10(adj.P.Val),
  text = ~protein,
  type = 'scatter',
  mode = 'markers',
  marker = list(color = 'rgba(255, 0, 0, .5)', size = 10)
)

#<---gsea emap plot---->
enrichplot::emapplot(interested_cyto_skel_gsea_go,
         showCategory = 17,
         cex_label_category = 1)




enrichplot::cnetplot(dis_normal_organ.go,
         foldChange=dis_normal_organ.go@geneList)

#<---quick visual of most interactive go terms
ggplot2::ggplot(most_enriched_go,
                aes(x=stats::reorder(Description,enrichmentScore),
                    y=enrichmentScore))+
  #geom_bar(stat="identity",fill="steelblue")+
  geom_point(shape=4,size=4,stroke=2,color="steelblue")+
  coord_flip()+
  labs(title = "Most enriched Go terms",
       x="Go terms",
       y="Enrichment Score")

#<--- see genes in up vs downregulated
ggplot2::ggplot(prot_to_organelle, 
       aes(x = reorder(protein, logFC),
           y = logFC,
           fill = logFC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  scale_fill_manual(values = c("red", "blue"), labels = c("Suppressed", "Activated")) +
  labs(x = "protein",
       y = "log Fold Change",
       fill = "Expression",
       caption = "Differentially expressed proteins in protein to organelle transport") +
  theme_minimal()
