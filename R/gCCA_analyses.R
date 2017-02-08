#install.packages("RGCCA")
#library(RGCCA)
#install.packages("pracma")

rm(list=ls())

library(gdata)
library(pracma)



####################################
#### prepare Xi and Ki matrices ####

load("cardiovascular.RData")

# this is performing step 1,2 and 3 (just prepare data)

rn<-sort(unique(c(rownames(X1),rownames(X2),rownames(X3))))
m<-length(rn)  # get the maximum number of individuals 

X1<-merge(data.frame(rn=rn),cbind(rn=rownames(X1),as.data.frame(X1)),by="rn",all.x=TRUE)
rownames(X1)<-X1$rn
X1<-remove.vars(X1,"rn")
ww<-apply(is.na(X1),1,all)
X1[ww,]<-0
K1<-matrix(0,nrow=m,ncol=m)
diag(K1)[!ww]<-1

X2<-merge(data.frame(rn=rn),cbind(rn=rownames(X2),as.data.frame(X2)),by="rn",all.x=TRUE)
rownames(X2)<-X2$rn
X2<-remove.vars(X2,"rn")
X2[ww<-apply(is.na(X2),1,all),]<-0
K2<-matrix(0,nrow=m,ncol=m)
diag(K2)[!ww]<-1

X3<-merge(data.frame(rn=rn),cbind(rn=rownames(X3),as.data.frame(X3)),by="rn",all.x=TRUE)
rownames(X3)<-X3$rn
X3<-remove.vars(X3,"rn")
X3[ww<-apply(is.na(X3),1,all),]<-0
K3<-matrix(0,nrow=m,ncol=m)
diag(K3)[!ww]<-1

X<-lapply(list(X1,X2,X3),as.matrix)
K<-list(K1,K2,K3)

#########################
#### start algorithm ####

n<-length(X) # number of tables
m<-nrow(X[[1]]) # maximum number of different individuals
p<-sapply(X,ncol) # number of variables per table
numvars<-min(p) # minimum number of variables


## compute eigen values and vectors

M<-Ksum<-matrix(0,nrow=m,ncol=m)
for (i in 1:n){
  Xi<-X[[i]]
  Ki<-K[[i]]
  Mi<-Ki%*%Xi%*%pinv(t(Xi)%*%Ki%*%Xi)%*%t(Xi)%*%Ki
  M<-M+Mi
  Ksum<-Ksum+Ki
}
Ksum05<-Ksum
diag(Ksum05)<-diag(Ksum05)^(-.5)
M<-Ksum05%*%M%*%Ksum05


eig<-eigen(M)
Yast<-eig$vectors
lambda<-eig$values
Y<-sqrt(n)*Ksum05%*%Yast

B<-A<-list()
for (i in 1:n){
  Xi<-X[[i]]
  Ki<-K[[i]]
  B[[i]]<-pinv(t(Y)%*%Ki%*%Y)%*%t(Y)%*%Ki%*%Xi
  A[[i]]<-pinv(t(Xi)%*%Ki%*%Xi)%*%t(Xi)%*%Ki%*%Y
}



As<-A    ## var(X*a)=1
for (i in 1:n){
  Ai<-A[[i]]
  Xi<-X[[i]]
  Ki<-K[[i]]
  vv<-1/apply(Ki%*%Xi%*%Ai,2,sd)
  vv<-cbind(rep(1,nrow(Ai)))%*%rbind(vv)
  As[[i]]<-A[[i]]*vv
}


#### finish algorithm ####
#########################


#
# Data visualization
#


## Select CpG (in heuristic way), |As| > 0.1
cutoff<-0.1
ww<-abs(As[[1]][,1])>cutoff | abs(As[[1]][,2])>cutoff
plot(As[[1]][,1],As[[1]][,2],col="white",pch=19,xlab=expression(As[1]),ylab=expression(As[2]))
abline(h=c(-cutoff,cutoff),v=c(-cutoff,cutoff),lty=2,col="grey")
abline(h=0,v=0,lty=2,col="grey")
text(As[[1]][,1],As[[1]][,2],colnames(X[[1]]),col=ifelse(ww,"red","black"),cex=0.8)

selsnps<-colnames(X[[1]])[ww]
length(selsnps)


## Plot original variables on Y factors.

  # only represent the CpG with higher A coefficients
par(xpd=NA)
plot(c(B[[1]][1,selsnps],B[[2]][1,],B[[3]][1,]),c(B[[1]][2,selsnps],B[[2]][2,],B[[3]][1,]),col="white",pch=19,xlab=expression(Y[1]),ylab=expression(Y[2]))
arrows(0,0,B[[2]][1,],B[[2]][2,],col="blue")
arrows(0,0,B[[3]][1,],B[[3]][2,],col="red",len=0.1)
text(B[[1]][1,selsnps],B[[1]][2,selsnps],selsnps,col="black",cex=0.8)
text(B[[2]][1,],B[[2]][2,],colnames(B[[2]]),col="blue")
text(B[[3]][1,],B[[3]][2,],colnames(B[[3]]),col="red")
par(xpd=FALSE)
abline(v=0,h=0,col="grey")
title("Original variables on first two Y factors")

## First canonical variables vs First factor Y

par(mar=c(5,4,3,2),oma=c(0,0,1,0))
layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = TRUE), respect = FALSE)
plot(Yast[,1],K[[1]]%*%X[[1]]%*%A[[1]][,1],xlab=expression(paste(Y[1]^"*")),ylab=expression(U[1]),pch=19)
abline(mod<-lm(I(K[[1]]%*%X[[1]]%*%A[[1]][,1])~I(Yast[,1]))); mtext(bquote(R^2==.(round(summary(mod)$r.squared,4))))
plot(Yast[,1],K[[2]]%*%X[[2]]%*%A[[2]][,1],xlab=expression(paste(Y[1]^"*")),ylab=expression(U[2]),pch=19)
abline(mod<-lm(I(K[[2]]%*%X[[2]]%*%A[[2]][,1])~I(Yast[,1]))); mtext(bquote(R^2==.(round(summary(mod)$r.squared,4))))
plot(Yast[,1],K[[3]]%*%X[[3]]%*%A[[3]][,1],xlab=expression(paste(Y[1]^"*")),ylab=expression(U[3]),pch=19)
abline(mod<-lm(I(K[[3]]%*%X[[3]]%*%A[[3]][,1])~I(Yast[,1]))); mtext(bquote(R^2==.(round(summary(mod)$r.squared,4))))
title("First canonical variables vs First factor Y*",out=TRUE)                                                           


# Redundancy Index

RI<-NULL
numvars<-6
for (i in 1:n){
  Xi<-X[[i]]
  Ki<-K[[i]]
  RI.i<-sapply(1:numvars,function(k){
    Yk<-Y[,1:k,drop=FALSE]
    sum(diag(t(Xi)%*%Ki%*%Yk%*%pinv(t(Yk)%*%Ki%*%Yk)%*%t(Yk)%*%Ki%*%Xi))/sum(diag(t(Xi)%*%Ki%*%Xi))
  })
  RI<-cbind(RI,RI.i)
}
colnames(RI)<-paste("X",1:n,sep="")
rownames(RI)<-paste("Y",1:numvars,sep="")
RI<-cbind(RI,total=rowMeans(RI))
RI



# correlation between Yj and X[[i]]

r2ij<-sapply(1:n,function(i){
    Xi<-X[[i]]
    Ki<-K[[i]]
    diag(t(Y)%*%Ki%*%Xi%*%pinv(t(Xi)%*%Ki%*%Xi)%*%t(Xi)%*%Ki%*%Y)/diag(t(Y)%*%Ki%*%Y)
  }
)

rowMeans(r2ij)[1:6]
r2ij<-cbind(r2ij,rowMeans(r2ij))
colnames(r2ij)<-c(paste("X",1:n,sep=""),"Total")
rownames(r2ij)<-paste("Y",1:ncol(Y),sep="")
r2ij[1:6,]

lambda[1:6] # this should be the mean of correlations!!





