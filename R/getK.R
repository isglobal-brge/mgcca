
getK <- function(x, ids, m) {
  X <- merge(data.frame(ids=ids), x, by.x="ids", by.y = 0 , all.x=TRUE)[,-1]
  rownames(X)<-ids
  ww<-apply(is.na(X), 1, all)
  X[ww,]<-0
  K <- matrix(0,nrow=m,ncol=m)
  diag(K)[!ww]<-1
  ans <- list(X=X, K=K)
  ans
} 


