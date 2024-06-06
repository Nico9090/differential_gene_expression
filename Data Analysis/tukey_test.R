#tukey test
df<-data.frame(hemoglobin=rep(c("HbA","HbS"),each=15),
               Oxygen=c(94,98,89,120,111,121,101,102,109,108,107,
                        104,100,107,103,82,32,43,54,74,44,43,
                        65,43,66,43,23,60,34,72))
))
atest<-aov(Oxygen~hemoglobin,data=df)
atest_summary<-summary(atest)




