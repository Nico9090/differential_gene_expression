# Importing libraries
```{r}
packages<-c('tidyverse','DESeq2','modelr') #Loaded libraries in the console
```

# Importing the data( #RBM disrupted samples vs Control)

#RBM disrupted
```{r}
rbm1<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM6041892_RBM1.tsv.gz')

rbm2<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM6041893_RBM2.tsv.gz')


rbm3<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM6041894_RBM3.tsv.gz')

```

#Controls
```{r}
control1<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM6041889_NT1.tsv.gz')

control2<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM6041890_NT2.tsv.gz')

control3<-read_delim('~/Bioinformatics Projects/SCD_hemoglobin/GSM6041891_NT3.tsv.gz')

```

#Tidy

```{r}
# RBM disrupted
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db) 
library(AnnotationDbi)
ids <- rbm1$gene_id
```
