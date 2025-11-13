
long_gene <- counts_long %>%
  mutate(Group = ifelse(grepl("WT", Sample), "WT", "KO"))





ggplot(long_gene, aes(x = external_gene_name,
                      y = Count,
                      fill = Group)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_errorbar(stat = "summary", 
                fun.data = mean_se, 
                position = position_dodge(width = 0.8),
                width = 0.3) +
  ggtitle(label = " kidney",
          subtitle = "Cilium organization genes")+
  labs(
    x = "gene",
    y = "mRNA counts",
    fill = "Condition"
  ) +
  scale_fill_manual(
    name = "Condition",
    values = c("WT" = "#4E79A7", "KO" = "#E15759")
  ) +
  theme(
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15,),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "right"
  )
