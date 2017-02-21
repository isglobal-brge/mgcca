rm(list=ls())
load("/Users/harryodell-mac/Desktop/GitHub/GCCA/data/cardiovascular.Rdata")
x<-list(X1,X2,X3)
rm(X1,X2,X3)


# input x (list of matrices)

n<-length(x)
for (i in 1:n) assign(paste0("X", i), x[[i]])


n<-length(x)
for (i in 1:n) {
  assign(paste0("X", i), x[[i]])
}





rn<-sort(unique(c(rownames(X1),rownames(X2),rownames(X3)))) # unique elements 
m<-length(rn)  # get the maximum number of individuals 

