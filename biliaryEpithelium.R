source("helperFunctions.R")
be #biliary epithelial cells
day <- be %>%
  dplyr::filter(day == "E10.75") 
keep <- day %>%
  dplyr::filter(exp!=0) #genes that have expression

day10_75BP <- enrichGOHelperFunction(geneSymbols = keep$gene,
                                   orgDb = org.Mm.eg.db,
                                   ont = "BP")
day10_75MF <- enrichGOHelperFunction(geneSymbols = keep$gene,
                                     orgDb = org.Mm.eg.db,
                                     ont = "MF")
forString <- day10_75BP %>%
  dplyr::filter(Description == "epithelial tube morphogenesis")%>%
  dplyr::select(geneID) %>%
  str_split("/") %>%
  unlist

day10_75STRING <- stringDBHelperFunction(geneSymbols = forString)
plotVisNet(mapping = day10_75STRING$mapping,
           interactions = day10_75STRING$interactions)


em <- be %>%
  dplyr::filter(gene %in% forString)

where <- enrichGOHelperFunction(geneSymbols = forString,
                                orgDb = org.Mm.eg.db,
                                ont = "CC")

binding <- enrichGOHelperFunction(geneSymbols = forString,
                                orgDb = org.Mm.eg.db,
                                ont = "MF")
group_annotations <- tribble(
  ~cell_component_group, ~genes,
  "collagen-containing ECM", "Bmp7/Fgfr2/Col4a1/Dlg1/Vegfa/Gpc3",
  "RNA pol II transcription regulator complex", "Pbx1/Ncoa3/Tead2/Tead1/Hif1a",
  "mitochondrial crista", "Bcl2/Opa1",
  "growth cone","Abl1/Smo/Pak1/Hdac5",
  "site of polarized growth","Abl1/Smo/Pak1/Hdac5",
  "basement membrane","Col4a1/Dlg1/Vegfa"
)


bindingAnnotations <- tribble(
  ~binding, ~genes,
  "transcription coregulator binding","Pbx1/Map3k7/Foxp1/Tead2/Hdac5/Hif1a",
  "RNA polymerase II-specific DNA-binding transcription factor","Ncoa3/Arid1a/Foxp1/Hmga2/Hdac5/Hif1a/Med12",
  "histone modifying activity","Kdm5b/Ncoa3/Map3k7/Hdac5/Kdm2a",
  "glycosaminoglycan binding","Lgr4/Bmp7/Fgfr2/Sema5a/Vegfa",
  "ubiquitin-like protein ligase binding","Bcl2/Map3k7/Glmn/Hif1a"
)
gene_group_map <- bindingAnnotations %>%
  separate_rows(genes, sep = "/") %>%
  rename(gene = genes)

em_annotated <- em %>%
  left_join(gene_group_map, by = "gene") %>%
  mutate(binding = coalesce(binding, "other"))


toPlot <- em_annotated %>%
  dplyr::filter(binding != "other")

ggplot(toPlot, aes(x = dayNum, y = exp, color = gene, group = interaction(celltype, gene))) +
  geom_line(size = 1) +
  facet_wrap(~ binding) +
  theme_minimal(base_size = 14)+
  ggtitle("Liver : epithelial tube morphogenesis")+
  xlab("Developmental stage")+
  ylab("Log Normalized Expression")
