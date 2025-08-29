library(tidyverse)
humanArpkdGo <- read_csv("HumanArpkdKidneyGoTerms.csv")
del34KidneyGo <- read_csv("Pkhd1Del34KidneyGoTerms.csv")
del34LiverGo <-read_csv("Pkhd1Del34LiverGoTerms.csv")
l1KidneyGo <- read_csv("Pkhd1l1KidneyGoTerms.csv")
l1LiverGo <- read_csv("Pkhd1l1LiverGoTerms.csv")
p0pkhd1Go <- read_csv("pkhd1KoKidneyOlsonGoTerms.csv") #p0
p3pkhd1Go <- read_csv("pkhd1KoP3KidneyOlsonGoTerms.csv")
p12pkhd1Go <- read_csv("pkhd1KoP12KidneyOlsonGoTerms.csv")
del67KidneyGo <- read_csv("pkhd1Del367KidneyGoTerms.csv")




#for this part, keep checking until all tables have equal number of rows
#so they have same number of GO terms
nrow(humanArpkdGo) 
nrow(del34KidneyGo)
nrow(del34LiverGo) 
nrow(l1KidneyGo)
nrow(l1LiverGo) 
nrow(p0pkhd1Go)#least
nrow(p3pkhd1Go)#least
nrow(p12pkhd1Go)#least
nrow(del67KidneyGo) 


humanArpkdGo <- humanArpkdGo[humanArpkdGo$category %in%
                               l1LiverGo$category,]
del34KidneyGo <- del34KidneyGo[del34KidneyGo$category %in%
                                 l1LiverGo$category,]

del34LiverGo <- del34LiverGo[del34LiverGo$category %in%
                               l1LiverGo$category,]

l1KidneyGo <- l1KidneyGo[l1KidneyGo$category %in%
                           l1LiverGo$category,]
l1LiverGo <- l1LiverGo[l1LiverGo$category %in%
                         humanArpkdGo$category,]

p0pkhd1Go <- p0pkhd1Go[p0pkhd1Go$category %in%
                         l1LiverGo$category,]
p3pkhd1Go <- p3pkhd1Go[p3pkhd1Go$category %in%
                         l1LiverGo$category,]
p12pkhd1Go <- p12pkhd1Go[p12pkhd1Go$category %in%
                           l1LiverGo$category,]

del67KidneyGo <- del67KidneyGo[del67KidneyGo$category %in%
                                 l1LiverGo$category,]

#<--same number of rows so proceed ------->


goTibble <- tibble(
  goId = humanArpkdGo$category,
  goDesc= humanArpkdGo$term,
  pvalHuman = humanArpkdGo$over_represented_pvalue,
  pvaldel34Kidney = del34KidneyGo$over_represented_pvalue,
  pvaldel34Liver = del34LiverGo$over_represented_pvalue,
  pvall1Liver = l1LiverGo$over_represented_pvalue,
  pvall1Kidney = l1KidneyGo$over_represented_pvalue,
  pvalp0 = p0pkhd1Go$over_represented_pvalue,
  pvalp3 = p3pkhd1Go$over_represented_pvalue,
  pvalp12 = p12pkhd1Go$over_represented_pvalue,
  pvaldel67 = del67KidneyGo$over_represented_pvalue,
  logFcHuman = humanArpkdGo$avgLog2FC,
  logFcdel34Kidney = del34KidneyGo$avgLog2FC,
  logFcdel34Liver = del34LiverGo$avgLog2FC,
  logFcl1Kidney = l1KidneyGo$avgLog2FC,
  logFcl1Liver = l1LiverGo$avgLog2FC,
  logFcp0 = p0pkhd1Go$avgLog2FC,
  logFcp3 = p3pkhd1Go$avgLog2FC,
  logFcp12 = p12pkhd1Go$avgLog2FC,
  logFcdel67 = del67KidneyGo$avgLog2FC
)

goOfInterest <- goTibble %>%
  dplyr::filter(goId %in% c("GO:0046649",
                            "GO:0051249",
                            "GO:0002757",
                            "GO:0050853",
                            "GO:0042110",
                            "GO:0031294",
                            "GO:0042110",
                            "GO:0050870",
                            "GO:0051251"
                            ))

goOfInterest <- goTibble %>%
  dplyr::filter(goId %in% c(
                            "GO:0035329",
                            "GO:0000165"))


goPoints <- goOfInterest %>%
  dplyr::select(goDesc,
                logFcHuman, logFcdel34Kidney, logFcdel34Liver,
                logFcl1Kidney,logFcl1Liver,logFcp0,logFcp3,logFcp12,logFcdel67,
                pvalHuman, pvaldel34Kidney, pvaldel34Liver,
                pvall1Liver,pvall1Kidney,pvalp0,pvalp3,pvalp12,pvaldel67) %>%
  pivot_longer(
    cols = c(logFcHuman, logFcdel34Kidney, logFcdel34Liver,
             logFcl1Kidney,logFcl1Liver,logFcp0,logFcp3,logFcp12,logFcdel67,
             pvalHuman, pvaldel34Kidney, pvaldel34Liver,
             pvall1Liver,pvall1Kidney,pvalp0,pvalp3,pvalp12,pvaldel67),
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
      condition == "p0" ~ 6,
      condition == "p3" ~ 7,
      condition == "p12" ~ 8,
      condition == "del67" ~ 9
    ),
    condition = factor(condition, levels = c("Human", "del34Kidney",
                                             "del34Liver","l1Kidney",
                                             "l1Liver","p0","p3","p12","del67"))
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
      "MAPK cascade" = "#1f77b4", 
      "hippo signaling" = "red"
    )
  ) +
  scale_x_continuous(breaks = 1:9, labels = c("Human ARPKD Kidney",
                                                       "Pkhd1 del 3-4/del 3-4 Kidney",
                                                       "Pkhd1 del 3-4/del 3-4 Liver",
                                                       "Pkhd1l1 null Kidney",
                                                       "Pkhd1l1 null Liver",
                                                       "Pkhd1 KO P0 Kidney",
                                                       "Pkhd1 KO P3 Kidney",
                                                       "Pkhd1 KO P12 Kidney",
                                                       "Pkhd1 del 3-67 Kidney"))  +
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
library(knitr)
kable(goPoints, caption = "", digits = 3)
































