#<---Functions--->
packages<-c("tidyverse","limma","clusterProfiler")
base::lapply(packages,base::library,character.only=TRUE)
#<---limma
limma_fit_for_deg<-function(counts,
                            design,
                            meta_data,
                            sample_identifier_by){
  #<---sample_identifier_by= condition or column 2 of meta data
  
  colnames(design)<-levels(meta_data[[sample_identifier_by]])
  
  fit<-limma::lmFit(as.matrix(log2(counts+1)),design)
  
  #<--change PLDF - WTF when needed --->
  contrast_matrix <- limma::makeContrasts(PLDF - WTF,
                                   levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  deg_table <- limma::topTable(fit2, adjust = "fdr", number = Inf)
  
  return(deg_table)
}

#<--gene set enrichment analysis
gseago_fn<-function(deg_table,gene_col_name,fc_col,organism){
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(org.Mm.eg.db)
  geneList <- deg_table[[fc_col]]
  names(geneList) <- deg_table[[gene_col_name]]
  gse <- gseGO(geneList=sort(geneList,
                             decreasing = TRUE), 
               ont ="ALL", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  return(list(go_obj=gse))
}

#<---obtain affymetrix gene names
convert_affy_id.fn<-function(gene.ids,affy.org.db){
  library(mouse4302.db)
  library(rat2302.db)
  library(pd.rat230.2)
  df<-select(affy.org.db,
             keys=gene.ids,
             columns=c("SYMBOL","ENTREZID","GENENAME"))
  return(df)
  
}

#<---DEG table using DESeq2
DESeq2.fn<-function(raw.counts,meta.data,design){
  dds<-DESeqDataSetFromMatrix(countData=raw.counts,
                              colData = meta.data,
                              design= as.formula(paste("~", design)))
  dds <-DESeq(dds)
  resultsNames(dds)
  res <- results(dds)
  return(list(dds_res=res,dds=dds))
}
