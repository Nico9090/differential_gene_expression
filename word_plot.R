go_results_df<-readr::read_csv("GO_enrichment_chronic_stress.csv")

go_results_df$short_desc<-ifelse(
  nchar(go_results_df$Description) > 30,
  paste0(substr(go_results_df$Description,1,27),"..."),
  go_results_df$Description
)
png("wordcloud.png", width = 1200, height = 800, res = 150)
wordcloud(words = go_results_df$short_desc,
          freq = go_results_df$Count,
          min.freq = 1,
          scale = c(4,0.8),
          random.order = FALSE,
          colors=c("red4","blue4","yellow4"))
dev.off() 
