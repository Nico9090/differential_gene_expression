#<---Functions--->
packages<-c("tidyverse","limma","clusterProfiler",
            "plotly","DOSE","enrichplot","maEndToEnd","DESeq2"
            )
lapply(packages,base::library,character.only=TRUE)
#<---limma
limma_fit_for_deg<-function(counts,
                            design,
                            meta_data,
                            sample_identifier_by){
  #<---sample_identifier_by= condition or column 2 of meta data
  colnames(design)<-levels(meta_data[[sample_identifier_by]])
  fit<-limma::lmFit(as.matrix(log2(counts+1)),design)
  #<--change PLDF - WTF when needed --->
  contrast_matrix <- limma::makeContrasts(KO - WT,
                                   levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  deg_table <- limma::topTable(fit2, adjust = "fdr", number = Inf)
  return(deg_table)
}
gseago_fn<-function(deg_table,gene_col_name,fc_col,organism){
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
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
gene_set_enrichment <- function(deg_table,
                                gene_column,
                                fold_change_column,
                                organism_db,
                                ontology_class = "BP",
                                pvalue_cutoff = 0.05,
                                p_adjust_method = "BH",
                                key_type = "SYMBOL",
                                min_gs_size = 3,
                                max_gs_size = 800) {
  
  # Load required libraries
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("enrichplot", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required packages: clusterProfiler, enrichplot, ggplot2.")
  }
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  
  # Validate columns
  if (!all(c(gene_column, fold_change_column) %in% colnames(deg_table))) {
    stop("Provided column names not found in 'deg_table'.")
  }
  
  # Validate OrgDb
  if (is.character(organism_db)) {
    # Try to load the OrgDb object from the namespace
    if (!requireNamespace(organism_db, quietly = TRUE)) {
      stop("OrgDb package not found: ", organism_db)
    }
    organism_db <- get(organism_db, envir = asNamespace(organism_db))
  }
  if (!inherits(organism_db, "OrgDb")) {
    stop("'organism_db' must be an OrgDb object or a valid package name.")
  }
  
  # Build geneList
  geneList <- deg_table[[fold_change_column]]
  names(geneList) <- deg_table[[gene_column]]
  
  # Remove NAs and sort
  geneList <- geneList[!is.na(names(geneList)) & !is.na(geneList)]
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Run GSEA
  gse_result <- tryCatch({
    gseGO(
      geneList = geneList,
      ont = ontology_class,
      keyType = key_type,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = p_adjust_method,
      verbose = FALSE,
      OrgDb = organism_db
    )
  }, error = function(e) {
    warning("gseGO failed: ", e$message)
    return(NULL)
  })
  
  return(list(go_obj = gse_result))
}

convert_id.fn<-function(gene.ids,affy.org.db){
  library(mouse4302.db)
  library(rat2302.db)
  library(pd.rat230.2)
  library(hgu133plus2.db)
  library(hugene10stv1cdf)
  library(hugene10sttranscriptcluster.db)
  library(org.Hs.eg.db)
  df<-AnnotationDbi::select(affy.org.db,
             keys=gene.ids,
             columns=c("SYMBOL","ENTREZID","GENENAME"))
  return(df)
  
}
DESeq2.fn<-function(raw.counts,meta.data,design){
  dds<-DESeqDataSetFromMatrix(countData=raw.counts,
                              colData = meta.data,
                              design= as.formula(paste("~", design)))
  dds <-DESeq(dds)
  resultsNames(dds)
  res <- results(dds)
  return(list(dds_res=res,dds=dds))
}
