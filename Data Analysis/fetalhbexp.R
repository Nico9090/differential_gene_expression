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

write.table(dat, file = "dat.csv",
            sep = "\t", row.names = F)

#Read new file created with all the variables I want

dat<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/dat.csv')

genes<-dat$neg1.gene_id
counts<-c(dat$neg1.expected_count,dat$neg2.expected_count,dat$na1.expected_count
          ,dat$na2.expected_count, dat$nx1.expected_count, dat$nx2.expected_count)
new_dat<-data.frame(genes,dat$neg1.expected_count,dat$neg2.expected_count,dat$na1.expected_count
                    ,dat$na2.expected_count, dat$nx1.expected_count, dat$nx2.expected_count)
write.table(new_dat, file = "new_dat.csv",
            sep = "\t", row.names = F)

new_dat<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/new_dat.csv')
#Generating metadata
l1<-rep('Control',2)
l2<-rep('NFIA',2)
l3<-rep('NFIX',2)
label<-c(l1,l2,l3)

cells<-c('neg1','neg2','na1','na2','nx1','nx2')
metadata<-data.frame(id=cells,type=label)

newest_dat<-data.frame(rows=cells,col=c(genes,dat$neg1.expected_count,dat$neg2.expected_count,dat$na1.expected_count
                       ,dat$na2.expected_count, dat$nx1.expected_count, dat$nx2.expected_count))
#Running DESeq2
dds = DESeqDataSetFromMatrix(new_dat,metadata,~type)
