#Obtain counts
#==================================================================================
library(Rsubread)
#=================================================================================
gen_Counts<-function(sam_file,gtf_file,csv_file){
        SAM<-sam_file
        countMatrix<-featureCounts(files=SAM,
                                   annot.inbuilt=NULL,
                                   annot.ext=gtf_file,
                                   isGTFAnnotationFile = TRUE,
                                   isPairedEnd = TRUE
                                )
        counts<-countMatrix$counts
        write.csv(counts,
                  file = csv_file,
                  row.names = TRUE)
}
