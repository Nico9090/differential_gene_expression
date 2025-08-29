#Sankey Plot
library(ggplot2)
library(ggalluvial)
library(tidyverse)

goOfInterest <- c("GO:0042110",
                  "GO:0046649",
                  "GO:0050870",
                  'GO:0000165')
goOfInterest <- c("GO:0035329","GO:0000165")

humanArpkdGo <- read_csv("../")
del34KidneyGo <- read_csv("..")
del34LiverGo <-read_csv("../csv")
l1KidneyGo <- read_csv("../v")
l1LiverGo <- read_csv("../")
p0pkhd1Go <- read_csv("../") 
p3pkhd1Go <- read_csv("../")
p12pkhd1Go <- read_csv("../")
del67KidneyGo <- read_csv(".")


goOfInterest <- c(
    "GO:0050870",
    "GO:0046649",
    'GO:0042110')

humanInterest<-humanArpkdGo %>% 
  dplyr::filter(category %in% goOfInterest)
del34KidInterest<-del34KidneyGo %>% 
  dplyr::filter(category %in% goOfInterest)
del34LivInterest<-del34LiverGo %>% 
  dplyr::filter(category %in% goOfInterest)
l1KidInterest<-l1KidneyGo %>% 
  dplyr::filter(category %in% goOfInterest)
l1LivInterest<-l1LiverGo %>% 
  dplyr::filter(category %in% goOfInterest)
del67Interest<-del67KidneyGo %>% 
  dplyr::filter(category %in% goOfInterest)
p0Interest<-p0pkhd1Go %>% 
  dplyr::filter(category %in% goOfInterest)
p3Interest<-p3pkhd1Go %>% 
  dplyr::filter(category %in% goOfInterest)
p12Interest<-p12pkhd1Go %>% 
  dplyr::filter(category %in% goOfInterest)


humanInterest <- dplyr::rename(humanInterest, goDesc = term)
del34KidInterest <- dplyr::rename(del34KidInterest, goDesc = term)
del34LivInterest <- dplyr::rename(del34LivInterest, goDesc = term)

l1KidInterest <- dplyr::rename(l1KidInterest, goDesc = term)
l1LivInterest <- dplyr::rename(l1LivInterest, goDesc = term)
del67Interest <- dplyr::rename(del67Interest, goDesc = term)

p0Interest <- dplyr::rename(p0Interest, goDesc = term)
p3Interest <- dplyr::rename(p3Interest, goDesc = term)
p12Interest <- dplyr::rename(p12Interest, goDesc = term)





goTibble <- humanInterest %>%
  dplyr::select(category, goDesc, pvalHuman = over_represented_pvalue, logFcHuman = avgLog2FC) %>%
  dplyr::left_join(
    del34KidInterest %>%
      dplyr::select(category, goDesc, pvaldel34Kidney = over_represented_pvalue, logFcdel34Kidney = avgLog2FC),
    by = c("category", "goDesc")
  ) %>%
  dplyr::left_join(
    del34LivInterest %>%
      dplyr::select(category, goDesc, pvaldel34Liver = over_represented_pvalue, logFcdel34Liver = avgLog2FC),
    by = c("category", "goDesc")
  )%>%
  dplyr::left_join(
    l1KidInterest %>%
      dplyr::select(category, goDesc, pvall1Kidney = over_represented_pvalue, logFcl1Kidney = avgLog2FC),
    by = c("category", "goDesc")
  )%>%
  dplyr::left_join(
    l1LivInterest %>%
      dplyr::select(category, goDesc, pvall1Liver = over_represented_pvalue, logFcl1Liver = avgLog2FC),
    by = c("category", "goDesc")
  )%>%
  dplyr::left_join(
    del67Interest %>%
      dplyr::select(category, goDesc, pvaldel67 = over_represented_pvalue, logFcdel67 = avgLog2FC),
    by = c("category", "goDesc")
  )%>%
  dplyr::left_join(
    p0Interest %>%
      dplyr::select(category, goDesc, pvalp0 = over_represented_pvalue, logFcp0 = avgLog2FC),
    by = c("category", "goDesc")
  )%>%
  dplyr::left_join(
    p3Interest %>%
      dplyr::select(category, goDesc, pvalp3 = over_represented_pvalue, logFcp3 = avgLog2FC),
    by = c("category", "goDesc")
  )%>%
  dplyr::left_join(
    p12Interest %>%
      dplyr::select(category, goDesc, pvalp12 = over_represented_pvalue, logFcp12 = avgLog2FC),
    by = c("category", "goDesc")
  )

goPoints <- goTibble %>%
  pivot_longer(
    cols = c(logFcHuman, logFcdel34Kidney, logFcdel34Liver,logFcl1Kidney,logFcl1Liver,
             logFcdel67,logFcp0,logFcp3,logFcp12,
             pvalHuman, pvaldel34Kidney, pvaldel34Liver,pvall1Kidney,pvall1Liver,pvaldel67,
             pvalp0,pvalp3,pvalp12
    ),
    names_to = c(".value", "condition"),
    names_pattern = "(logFc|pval)(.*)"
  ) %>%
  mutate(
    x = case_when(
      condition == "Human" ~ 1,
      condition == "del34Kidney" ~ 2,
      condition == "del34Liver" ~ 3,
      condition == "l1Kidney" ~ 4,
      condition == "l1Liver" ~ 5,
      condition == "del67" ~ 6,
      condition == "p0" ~ 7,
      condition == "p3" ~ 8,
      condition == "p12" ~ 9
      
    ),
    condition = factor(condition, levels = c("Human", "del34Kidney",
                                             "del34Liver","l1Kidney",
                                             "l1Liver","del67","p0","p3",
                                             "p12"))
  )
ggplot(goPoints, aes(x = x, y = logFc, color = goDesc)) +
  geom_point(
    aes(size = -log2(pval)),
    position = position_dodge(width = 0.4),
    alpha = 0.8,#stroke = 4
  ) +
  geom_line(aes(group = goDesc), position = position_dodge(width = 0.4), alpha = 0.5)+
  scale_color_manual(
    name = "GO Term",
    values = c(
      "MAPK cascade"="blue4",
      "hippo signaling"="red"
    )
  ) +
  scale_x_continuous(breaks = 1:9, labels = c("Human ARPKD Kidney",
                                              "Pkhd1 del 3-4/del 3-4 Kidney",
                                              "Pkhd1 del 3-4/del 3-4 Liver",
                                              "Pkhd1l1 KO/KO Kidney",
                                              "Pkhd1l1 KO/KO Liver",
                                              "Pkhd1 del 3-67 Kidney",
                                              "Pkhd1 KO Olson P0",
                                              "Pkhd1 KO Olson P3",
                                              "Pkhd1 KO Olson P12"
  ))  +
  scale_size_continuous(name = "-log2 p-value") +
  theme_minimal() +
  labs(x = "", y = "Log2 FC", color = "GO Term") +
  geom_hline(yintercept = 0,      
             linetype = "dashed",
             color = "black",
             linewidth=0.8)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,face = "bold"),
    axis.text.y = element_text(vjust = 0.5, hjust = 1,face = "bold",size = 10),
    legend.position = "right"
  )
