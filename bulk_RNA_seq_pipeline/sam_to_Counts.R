#Obtain counts
#==================================================================================
library(Rsubread)

#=================================================================================
#SAM files

s11_SAM<-"../output.sam"
s12_SAM<-"S12.sam"
s13_SAM<-"S13.sam"
s14_SAM<-"S14.sam"
s15_SAM<-"S15.sam"
s16_SAM<-"S16.sam"

#=================================================================================
#Generate count matrices
s11_countMatrix<- featureCounts(files=s11_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s12_countMatrix<- featureCounts(files=s12_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s13_countMatrix<- featureCounts(files=s13_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s14_countMatrix<- featureCounts(files=s14_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s15_countMatrix<- featureCounts(files=s15_SAM,
                                annot.inbuilt = "hg38",                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s15_countMatrix<- featureCounts(files=s15_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
s16_countMatrix<- featureCounts(files=s16_SAM,
                                annot.inbuilt = "hg38",
                                isGTFAnnotationFile = FALSE,
                                isPairedEnd = TRUE)
#===================================================================================
#SAVE count matrices
s11_counts<-s11_countMatrix$counts
s12_counts<-s12_countMatrix$counts
s13_counts<-s13_countMatrix$counts
s14_counts<-s14_countMatrix$counts
s15_counts<-s15_countMatrix$counts
s16_counts<-s16_countMatrix$counts

write.csv(s11_counts, file = "s11_counts.csv", row.names = TRUE)
write.csv(s12_counts, file = "s12_counts.csv", row.names = TRUE)
write.csv(s13_counts, file = "s13_counts.csv", row.names = TRUE)
write.csv(s14_counts, file = "s14_counts.csv", row.names = TRUE)
write.csv(s15_counts, file = "s15_counts.csv", row.names = TRUE)
write.csv(s16_counts, file = "s16_counts.csv", row.names = TRUE)
#=================================================================================
csv_files <- list.files(pattern = "*.csv")

# Read and combine all the CSV files into one data frame
merged_data <- do.call(rbind, lapply(csv_files, read.csv))
write.csv(merged_data, "merged_data.csv", row.names = TRUE)

                                
