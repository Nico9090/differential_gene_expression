data<-read.delim('GSM6041889_NT1.tsv')
View(data)
data<-head(data) #Playing around with data transformation
library() #Add library needed. I worked with tidyverse
ggplot(data,mapping=aes(x=length,y=effective_length,size=expected_count))+
  geom_point()
data<-read.delim('GSM6041889_NT1.tsv')
View(data)
data<-select(data,length:FPKM) #Playing around with data transformation
head(data)
