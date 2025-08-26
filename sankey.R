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















humanMtMotor<-c("KIF4A", "KIF20A", "DYNLRB1", "KIF1A",
                "CENPE", "KIF2C", "KIF26B", "KIF15", "KIF7",
                "DNAH2", "KIF18B","KIFC1")
mtMotor<- humanKidney[humanKidney$gene.symbol %in%
                        humanMtMotor,]
mtMotorMouse <- mouseKidney[tolower(mouseKidney$gene.symbol) %in%
                              tolower(humanMtMotor),] 
mtMotorHuman <- mtMotor[tolower(mtMotor$gene.symbol) %in%
                        tolower(mtMotorMouse$gene.symbol),]


df <- tibble::tibble(
  human        = mtMotorHuman$gene.symbol,
  mouse        = mtMotorMouse$gene.symbol,
  log2FC_human = mtMotorHuman$log2FoldChange,
  log2FC_mouse = mtMotorMouse$log2FoldChange,
  Direction_human = ifelse(mtMotorHuman$log2FoldChange > 0, "Up", "Down")
)

ggplot(df,
       aes(axis1 = human, axis2 = mouse, y = abs(log2FC_human))) +
  geom_alluvium(aes(fill = Direction_human)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()



df_long <- df %>%
  tidyr::pivot_longer(
    cols = c(log2FC_human, log2FC_mouse),
    names_to = "Species", values_to = "log2FC"
  )

ggplot(df_long,
       aes(axis1 = human, axis2 = Species, axis3 = mouse, y = abs(log2FC))) +
  geom_alluvium(aes(fill = Species)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()























library(ggplot2)
library(dplyr)

# Example dataframe
df <- tibble::tibble(
  gene   = mtMotorHuman$gene.symbol,
  log2FC_mouse = mtMotorMouse$log2FoldChange,
  log2FC_human = mtMotorHuman$log2FoldChange
)

# Reshape to long format
df_long <- df %>%
  tidyr::pivot_longer(cols = c(log2FC_mouse, log2FC_human),
                      names_to = "Species", values_to = "log2FC")

# Order the x-axis
df_long$Species <- factor(df_long$Species,
                          levels = c("log2FC_mouse", "log2FC_human"),
                          labels = c("Mouse", "Human"))

# Plot
ggplot(df_long, aes(x = Species, y = log2FC, group = gene, color = gene)) +
  geom_line(linewidth = 1.2, alpha = 0.7) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal() +
  labs(y = "log2 Fold Change", x = "")











edges <- df %>%
  mutate(x_start = 1, x_end = 2,
         y_start = log2FC_mouse,
         y_end   = log2FC_human)

ggplot() +
  geom_curve(
    data = edges,
    aes(x = x_start, y = y_start,
        xend = x_end, yend = y_end,
        linewidth = abs(log2FC_human),   # thickness
        color = gene),
    curvature = 0.2, alpha = 0.6
  ) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Mouse", "Human")) +
  theme_minimal() +
  labs(y = "log2 Fold Change", x = "",
       caption = "Human vs mouse DEG in microtubule motor")


