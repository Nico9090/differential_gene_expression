#Creating models for sample data
#Importing packages
packages<-c('tidyverse','modelr')
lapply(packages, library, character.only=TRUE)

#Creating sample data
x<-sample(-5:5,30,replace=TRUE)
y<-round(rnorm(30),4)

data<-data.frame(x,y)

#Plot the data

plot<-ggplot(data=data,mapping=aes(x,y))+
  geom_point()

#Linear enough so write model
l_model<-function(b,m){
  m*data$x + b
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
# Root mean square deviation
rms<-function(model_y_val,actual_data){
  (model_y_val[i]-actual_data[i])^2
}
