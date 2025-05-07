library(Matrix)
values<-c(1,2,3,4,5,6)
#________________________________________________________________________
matrix<-matrix(values,
              nrow=2,
              ncol=3) #product of rows and columns must equal # of values
#This matrix is a dense matrix. When you have a dataset with more non-zero values, use this

#________________________________________________________________________
#Sparse matrix (dcgMatrix)
sparse_matrix<-Matrix(matrix,
                      sparse=TRUE) #use this when there are many 0 values, like in RNA sequencing
