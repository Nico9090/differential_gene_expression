# Load libraries
packages<-c('tidyverse','EnhancedVolcano','DESeq2')
lapply(packages,library,character.only=TRUE)

#Import Data
neg1<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM5473895_Negative_gRNA_rep1.hg38.tsv.gz')
neg2<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM5473896_Negative_gRNA_rep2.hg38.tsv.gz')
na1<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM5473897_NFIA_gRNA1_rep1.hg38.tsv.gz')
na2<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM5473898_NFIA_gRNA1_rep2.hg38.tsv.gz')
nx1<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM5473901_NFIX_gRNA1_rep1.hg38.tsv.gz')
nx2<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM5473902_NFIX_gRNA1_rep2.hg38.tsv.gz')

#Make data frame of data of choice
dat<-data.frame(neg1$gene_id,neg1$expected_count,neg2$expected_count,na1$expected_count,
       na2$expected_count,nx1$expected_count,nx2$expected_count)


#Generating metadata
l1<-rep('Control',2)
l2<-rep('NFIA',2)
l3<-rep('NFIX',2)
label<-c(l1,l2,l3)

cells<-c('neg1','neg2','na1','na2','nx1','nx2')
metadata<-data.frame(id=cells,type=label)

#Running DESeq2
dds = DESeqDataSetFromMatrix(dat,metadata,~type)
