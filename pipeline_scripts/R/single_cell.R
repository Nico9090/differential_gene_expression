#Libraries___________________________________________________________________________________
library(dplyr)
library(Seurat)
library(patchwork)
#_____________________________________________________________________________________________
#Create seurat objects for each condition
ko1.data<-Read10X(data.dir="path/to/knockout/matrix")
wt1.data <-Read10X(data.dir = "path/to/control/matrix")
#______________________________________________________________________________________________
ko1_data<-CreateSeuratObject(counts = ko1.data, 
                             project = "SnRNA-Seq Analysis", 
                             min.cells =3, 
                             min.features = 200)
