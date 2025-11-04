#ADPKD data
norm_counts <- read_csv("GSE289843_global.normalized_counts.csv")
p1 <- norm_counts %>%
  dplyr::select(...1,GC126323,GC126324,GC126325,GC126326,
                GC117831,GC117832,GC117833,GC117833) #cyst A to D

p1 <- p1 %>%
  dplyr::rename(gene_name = ...1)
counts_long <- p1 %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "Sample",
    values_to = "Count"
  )

immunglobulin_genes <- counts_long %>%
  dplyr::filter(str_starts(gene_name,"IGKV1-12|IGKV1-16|IGKV1-27"))
long_gene <- immunglobulin_genes %>%
  mutate(Group = ifelse(grepl("GC12", Sample), "WT", "KO"))

ggplot(long_gene, aes(x = Group, y = Count, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~ gene_name, scales = "free_y") +
  theme_bw() +
  labs(title = "Raw Counts per Gene",
       y = "Read Counts",
       x = "Group")
