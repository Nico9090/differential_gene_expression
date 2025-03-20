#Functions
library(Rsubread)
library(biomaRt)
library(AnnotationDbi)
inPath = "path/to/data"
s12_SAM<-file.path(inPath,"/S12.sam")
s12_countMatrix<- featureCounts(files=s12_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s12_counts<-s12_countMatrix$counts

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = s12_countMatrix$annotation$GeneID,  # This is a vector of gene IDs
                   mart = mart)
