#importing needed libraries

#Adding and viewing data 
data<-read.delim("GSE232221_raw_counts_GRCh38.p13_NCBI.tsv")
View(data)

#Add annotation sheet

annot<-read.delim("Human.GRCh38.p13.annot.tsv")
View(annot)

# Viewing the type of arrangement needed momentarily
dataset<-data.frame(x=annot$Symbol,HbA_sample1=data$GSM7321059)
print(head(dataset))

#Preparing volcano plot
dds = DESeqDataSetFromMatrix(dataset,Metadata,~type)
