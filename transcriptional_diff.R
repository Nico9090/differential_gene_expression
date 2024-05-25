#importing needed libraries
library(tidyverse)
library(EnhancedVolcano)
library(DESeq2)
#Adding and viewing data 
data<-read.delim("~/Bioinformatics Projects/GSE232221_raw_counts_GRCh38.p13_NCBI.tsv")
View(data)
d=data[1:37,]
#Add annotation sheet

annot<-read.delim("~/Bioinformatics Projects/Human.GRCh38.p13.annot.tsv")
View(annot)
a=annot[1:37,]

#Creating a new data frame of selected data
table1<-data.frame(d,a)
Genes<-d$GSM7321059
Counts<-a$Symbol
y<-fct_reorder(Counts,Genes)


#Applying levels to the gene symbols so that graph is clear
ggplot(data=table1,mapping=aes(x=Genes,fct_reorder(Counts,Genes)))+
  geom_point()
