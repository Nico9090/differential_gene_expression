#LIBRARIES______________________________________________________________________
#ADD GENE NAMES TO RAW DATA
#_______________________________________________________________________________
library(biomaRt)
library(goseq)
#_______________________________________________________________________________
#_______________________________________________________________________________
raw_data<-"" #csv file
counts<-read_csv(raw_data)
genome<-"" #enseml genome. Ex:mmusculus_gene_ensembl

search_gene_names<-function(genome,raw_data){
  ensembl <- useMart("ensembl", 
                   dataset = genome)
  ens_ID<-"" #column name for ensembl IDs in raw data
  IDs<- raw_data$ens_ID
  names<-getBM(attributes = c("ensembl_gene_id", 
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = IDs,
                   mart = ensembl)
  gene_marked_data<-merge(raw_data, 
                          names, 
                          by.x = ens_ID, 
                          by.y = "ensembl_gene_id", 
                          all.x = TRUE)
  write_csv(gene_marked_data,"name.csv") #change to desired name
}
