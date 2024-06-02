#Import libraries
```{r}
pkg<-c('modelr','tidyverse')
lapply(pkg,library,character.only=TRUE)
```

#Creating sample dataset to generate the model of
```{r}
#Creating a data frame to obtain a random distribution of x and y values
#The plot the data with ggplot
dat<-tibble(
  x=runif(100,-5,5),
  y=runif(100,20,50)
)

plot<-ggplot(dat,aes(x,y))+
  geom_point()
```

#Generating models for the data
```{r}
mdl<-tibble(
  s=runif(100,-10,10),
  i=runif(100,20,50)
)

lines<-ggplot(dat,aes(x,y))+
  geom_abline(data=mdl,aes(intercept=i,slope=s),alpha=1/4)+
  geom_point()
```

#To find the best model of them all, emply the RMSD
```{r}
rmsd<-function(dataset,model){
  sqrt((mean(dataset$y-model$s))^2)
}
```

