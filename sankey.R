df <- tibble::tibble(
  gene   = mtMotorHuman$gene.symbol,
  log2FC_mouse = mtMotorMouse$log2FoldChange,
  log2FC_human = mtMotorHuman$log2FoldChange
)

edges <- df %>%
  mutate(x_start = 1, x_end = 2,
         y_start = log2FC_mouse,
         y_end   = log2FC_human)

ggplot() +
  geom_curve(
    data = edges,
    aes(x = x_start, y = y_start,
        xend = x_end, yend = y_end,
        linewidth = abs(log2FC_mouse),   # thickness
        color = gene),
    curvature = 0.2, alpha = 0.6
  ) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Mouse", "Human")) +
  theme_minimal() +
  labs(y = "log2 Fold Change", x = "")


