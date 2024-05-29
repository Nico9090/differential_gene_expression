#Importing packages
packages<-c('tidyverse','modelr')
lapply(packages, library, character.only=TRUE)

#Creating sample data
x<-sample(-5:5,30,replace=TRUE)
y<-round(rnorm(30),4)

scores<-data.frame(x,y)

#Plot the data

plot<-ggplot(data=scores,mapping=aes(x,y))+
  geom_point()

#Linear enough so write model
l_model<-function(b,m){
  m*scores$x + b
}
trendline=l_model(0.5,2)

#Plot the trend model

plot +
  geom_abline(data=trendline)

#Generate many models
models<-tibble(x=runif(300,-5,5),y=runif(300,-2,3))
#Plot many models
plot+
  geom_abline(data=models, aes(intercept=x,slope=y,
                               alpha=1/4))
# Root mean square deviation for model we created
rms<-function(model_line,dataset){
  r1<-(abs(scores$y -model_line))^2
  sum=0
  for (i in r1){
    sum=sum + i
  }
  sqrt(sum)
}
rms(trendline,scores)

#Distances for all the bad models using purrr
distances<-function(dataset){
  models %>%
    purrr::map2_dbl(x,y,dataset)
  
}
distances(scores)

