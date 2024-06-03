# Load libraries
some_packages<-c('EnhancedVolcano','tidyverse','DESeq2')
lapply(some_packages, library, character.only=TRUE)

#Import and view data
counts<-read.delim('GSE232221_raw_counts_GRCh38.p13_NCBI.tsv') #Raw counts
atable<-read.delim('Human.GRCh38.p13.annot.tsv') #Human gene annotation table
View(counts)
View(atable)

#Data wrangling

counts<-counts[1:39376,]
atable<-atable[1:39376,]

control<-rep(c('Control'),5)
ex<-rep(c('Experimental'),5)
label<-c(control,ex)
cells<-c('GSM7321059','GSM7321060','GSM7321061','GSM7321062','GSM7321063'
         ,'GSM7321064','GSM7321065','GSM7321066','GSM7321067','GSM7321068')
metadata<- data.frame(id=cells,type=label)
genes=counts[,cells]

#DESeq2
dds = DESeqDataSetFromMatrix(genes,metadata,~type)
dds <- DESeq(dds)
res=results(object = dds, contrast = c('type','Control','Experimental'),
            pAdjustMethod = 'holm', alpha = 0.000001)
row.names(res)=atable$Symbol
summary(res)
EnhancedVolcano(res,
                lab = atable$Symbol,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-6',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)
