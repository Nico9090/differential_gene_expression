humanKidney<- humanKidney %>%
  mutate(color = case_when(
    baseMean>=1000 ~"green4",
    padj <= 0.05 & log2FoldChange >= 1.5  ~ "blue4",
    padj <= 0.05 & log2FoldChange <= -1.5 ~ "red4",
    
    TRUE                                  ~ "black" 
  ))
humanKidney <- humanKidney %>%
  mutate(group = case_when(
    baseMean>=1000 ~"Mean counts >= 1000",
    padj <= 0.05 & log2FoldChange >= 1.5  ~ "Up-regulated",
    padj <= 0.05 & log2FoldChange <= -1.5 ~ "Down-regulated",
    
    TRUE                                  ~ "Not significantly expressed"
  ))


ggplot2::ggplot(data = del34_LIV,
                mapping = aes(x = log2FoldChange,
                              y = -log10(padj),
                              color = group)) +
  xlim(-15,15)+
  geom_point(shape=3,size=2.5,stroke = 1.5)+ 
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = "dashed",  # or "dashed"
             color = "black",
             linewidth=0.8)+      # color of the line
  geom_hline(yintercept = 1.3,      # e.g., -log10(0.05) ≈ 1.3
             linetype = "dashed",
             color = "black",
             linewidth=0.8)+
  scale_color_manual(values = c("Up-regulated" = "blue4",
                                "Down-regulated" = "red4",
                                "Mean counts >= 1000" = "green4",
                                "Not significantly expressed"="black"))+
  labs(color="Expression",title="Pkhd1 del 3-4/del 3-4 Liver")
  
  
