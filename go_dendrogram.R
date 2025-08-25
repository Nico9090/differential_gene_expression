del34LIV_GO<-pairwise_termsim(del34LIV_GO)
goPlot_dendrogram <- function(simMatrix, 
                              numClusters,
                              predictedClusters){
  library(ape)
  library(ggtree)
  library(dplyr)
  
  # distance + clustering
  distMatrix <- as.dist(1 - simMatrix)
  tierCluster <- hclust(distMatrix, method = "average")
  phyloTree <- as.phylo(tierCluster)
  
  # assign clusters
  clusters <- cutree(tierCluster, k = numClusters)
  clusterTable <- tibble(
    label = names(clusters),   # ggtree expects 'label'
    cluster = as.factor(clusters[names(clusters)])
  )
  
  # map numeric clusters → predicted cluster names
  clusterTable$cluster <- factor(
    clusterTable$cluster,
    levels = seq_len(numClusters),
    labels = names(predictedClusters)
  )
  
  # compute MRCA nodes, skipping singletons
  nodes <- clusterTable %>%
    group_by(cluster) %>%
    summarise(node = {
      mrca <- getMRCA(phyloTree, label)
      if (is.null(mrca)) NA_integer_ else mrca
    }, .groups="drop") %>%
    filter(!is.na(node))
  
  nodes$cluster_node <- nodes$node
  
  # plot
  treePlot <- ggtree(phyloTree) %<+% clusterTable +
    geom_tiplab(size = 3, family = "Consolas") +   # GO terms on tips
    theme(
      legend.position = "right",
      text = element_text(family = "Consolas"),
      axis.text = element_text(family = "Consolas"),
      legend.text = element_text(family = "Consolas"),
      legend.title = element_text(family = "Consolas")
    ) +
    scale_color_manual(values = predictedClusters, name = "Clusters") +
    scale_fill_manual(values = predictedClusters, name = "Clusters") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.6))) +
    geom_hilight(data = nodes,
                 aes(node = cluster_node, fill = cluster),
                 alpha = 0.2)
  
  return(treePlot)
}

del34Kidney_clusters <- c(
  "primary cilia component" = "#E41A1C",  # red
  "abnormal protein folding"             = "#377EB8",  # blue
  "immune response"                     = "#4DAF4A"   # green
)

goPlot_dendrogram(simMatrix = del34LIV_GO@termsim,
                  numClusters = 3,
                  predictedClusters = del34Kidney_clusters) +
  labs(caption = "liver")
