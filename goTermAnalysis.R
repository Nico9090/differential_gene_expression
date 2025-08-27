del34KidGo <- read_csv("name1GoTerms.csv")
del34LivGo <- read_csv("name2GoTerms.csv")

l1KidGo <- read_csv("name3GoTerms.csv")
l1LivGo <- read_csv("name4GoTerms.csv")
#so they have same number of GO terms
nrow(del34KidGo)
nrow(del34LivGo)
nrow(l1KidGo)
nrow(l1LivGo) #least

del34KidGo <- del34KidGo[del34KidGo$category %in%
                           l1LivGo$category,]
del34LivGo <- del34LivGo[del34LivGo$category %in%
                           l1LivGo$category,]
l1KidGo <- l1KidGo[l1KidGo$category %in%
                           l1LivGo$category,]
goTibble <- tibble(
  goId = del34KidGo$category,
  goDesc= del34KidGo$term,
  pvalDel34Kid = del34KidGo$over_represented_pvalue,
  pvalDel34Liv = del34LivGo$over_represented_pvalue,
  pvalL1Kid = l1KidGo$over_represented_pvalue,
  pvalL1Liv = l1LivGo$over_represented_pvalue,
  logFcDel34Kid = del34KidGo$avgLog2FC,
  logFcDel34Liv = del34LivGo$avgLog2FC,
  logFcL1Kid = l1KidGo$avgLog2FC,
  logFcL1Liv = l1LivGo$avgLog2FC
)

goOfInterest <- goTibble %>%
  dplyr::filter(goId %in% c("GO:0046649",
                            "GO:0051249",
                            "GO:0002757",
                            "GO:0050853",
                            "GO:0042110",
                            "GO:0031294",
                            "GO:0042110",
                            "GO:0050870"))
edges <- goOfInterest %>%
  mutate(xStart = 1,xP1=2,xP2=3, xEnd = 4,
         yStart = 2**logFcDel34Kid,yP1=2**logFcL1Kid,
         yP2=2**logFcL1Liv,
         yEnd   = 2**logFcDel34Liv)

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
  scale_x_continuous(breaks = c(1,2,3,4), labels = c("Del 3-4 Kidney",
                                                     "L1 Kidney",
                                                     "L1 Liver",
                                                     "Del 3-4 Liver")) +
  theme_minimal() +
  labs(y = "Significance", x = "",
       caption = "Go Terms")










library(dplyr)
library(tidyr)
library(ggplot2)

edges_long <- edges %>%
  left_join(goLong, by = "goDesc")

edges_points <- edges_long %>%
  # build start points
  transmute(goDesc, pval, condition,
            x = xStart, y = yStart, end = "Start") %>%
  bind_rows(
    edges_long %>%
      transmute(goDesc, pval, condition,
                x = xEnd, y = yEnd, end = "End")
  )
# Now you have columns: x, y, end (Start/End), plus goDesc, pval, condition

# 2. Plot with geom_point()
ggplot(edges_points, aes(x = x, y = y,
                         size = -log10(pval),
                         color = goDesc)) +
  geom_point(alpha = 0.7,shape = 3, stroke=4) +
  scale_size(range = c(1, 6)) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c("Del 3-4 Kidney",
                                "L1 Kidney",
                                "L1 Liver",
                                "Del 3-4 Liver")) +
  theme_minimal() +
  labs(y = "Fold Change", x = "",
       caption = "GO Terms")
