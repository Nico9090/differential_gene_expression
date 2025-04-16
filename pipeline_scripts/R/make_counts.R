#Libraries#_______________________________________________________________________
library(Rsubread)
#_________________________________________________________________________________
gen_Counts<-function(
  sam_file,
  gtf_file){
        SAM<-sam_file
        countMatrix<-featureCounts(files=SAM,
                                   annot.inbuilt = NULL,
                                   annot.ext= gtf_file
                                )
        counts<-countMatrix$counts
        write.csv(counts,
                  file = paste0(sub("\\..*$",
                                    "",
                                    sam_file[1]),
                                "auto_counts.csv"),
                  row.names = TRUE)
}

files<-list() #add the names of the files
gtf<-"" #add gtf file
for(file in files){
  gen_Counts(file,
             gtf)
}
