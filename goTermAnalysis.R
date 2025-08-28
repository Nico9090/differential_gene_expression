Pkd1KoGo <- read_csv("Pkd1KoMouseKidneyGoTerms.csv")
glis3KoGo <- read_csv("glis3KoMouseKidneyGoTerms.csv")
doubleKoGo <-read_csv("DoubleKoMouseKidneyGoTerms.csv")

#so they have same number of GO terms
nrow(Pkd1KoGo)
nrow(glis3KoGo)
nrow(doubleKoGo) #least

Pkd1KoGo <- Pkd1KoGo[Pkd1KoGo$category %in%
                       doubleKoGo$category,]
glis3KoGo <- glis3KoGo[glis3KoGo$category %in%
                         doubleKoGo$category,]

goTibble <- tibble(
  goId = Pkd1KoGo$category,
  goDesc= Pkd1KoGo$term,
  pvalPkd1 = Pkd1KoGo$over_represented_pvalue,
  pvalglis3 = glis3KoGo$over_represented_pvalue,
  pvalDoubleKo = doubleKoGo$over_represented_pvalue,
  logFcPkd1 = Pkd1KoGo$avgLog2FC,
  logFcglis3 = glis3KoGo$avgLog2FC,
  logFcDoubleKo = doubleKoGo$avgLog2FC
)

goOfInterest <- goTibble %>%
  dplyr::filter(goId %in% c("GO:0043410",
                            "GO:0043408",
                            "GO:0000165",
                            "GO:1900745",
                            "GO:1900744",
                            "GO:0038066"))
edges <- goOfInterest %>%
  mutate(xStart = 1,xP1=2, xEnd = 3,
         yStart = 2**logFcPkd1,yP1=2**logFcglis3,
         yEnd   = 2**logFcDoubleKo)

ggplot() +
  geom_curve(
    data = edges,
    aes(x = xStart, y = yStart,
        xend = xEnd, yend = yEnd,
        linewidth = 0.5,   # thickness
        color = goDesc),
    curvature = 0.2, alpha = 0.6
  ) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("Pkd1 KO",
                                                     "Glis3 KO",
                                                     "Double KO")) +
  theme_minimal() +
  labs(y = "LogFC", x = "",
       caption = "Go Terms")



goPoints <- goOfInterest %>%
  dplyr::select(goDesc, logFcPkd1, logFcglis3, logFcDoubleKo) %>%
  pivot_longer(
    cols = starts_with("logFc"),
    names_to = "condition",
    values_to = "logFC"
  ) %>%
  mutate(
    x = case_when(
      condition == "logFcPkd1" ~ 1,
      condition == "logFcglis3" ~ 2,
      condition == "logFcDoubleKo" ~ 3
    ),
    condition = factor(condition, levels = c("logFcPkd1", "logFcglis3", "logFcDoubleKo"))
  )
ggplot(goPoints, aes(x = x, y = logFC, color = goDesc)) +
  geom_point(shape=4,stroke=4,
    position = position_dodge(width = 0.4),
    size = 3, alpha = 4
  ) +
  geom_line(aes(group = goDesc), position = position_dodge(width = 0.4), alpha = 0.5)+
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Pkd1 KO", "Glis3 KO", "Double KO")
  ) +
  theme_minimal() +
  labs(x = "", y = "log2 Fold Change", color = "GO Term") +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "right"
  )




goPoints <- goOfInterest %>%
  dplyr::select(goDesc,
         logFcPkd1, logFcglis3, logFcDoubleKo,
         pvalPkd1, pvalglis3, pvalDoubleKo) %>%
  pivot_longer(
    cols = c(logFcPkd1, logFcglis3, logFcDoubleKo,
             pvalPkd1, pvalglis3, pvalDoubleKo),
    names_to = c(".value", "condition"),
    names_pattern = "(logFc|pval)(.*)"
  ) %>%
  mutate(
    x = case_when(
      condition == "Pkd1" ~ 1,
      condition == "glis3" ~ 2,
      condition == "DoubleKo" ~ 3
    ),
    condition = factor(condition, levels = c("Pkd1", "glis3", "DoubleKo"))
  )
ggplot(goPoints, aes(x = x, y = logFc, color = goDesc)) +
  geom_point(
    aes(size = -log10(pval)),
    position = position_dodge(width = 0.4),
    alpha = 0.8
  ) +
  geom_line(aes(group = goDesc), position = position_dodge(width = 0.4), alpha = 0.5)+
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Pkd1 KO", "Glis3 KO", "Double KO")
  ) +
  scale_size_continuous(name = "-log10(p-value)") +
  theme_minimal() +
  labs(x = "", y = "log2 Fold Change", color = "GO Term") +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "right"
  )
