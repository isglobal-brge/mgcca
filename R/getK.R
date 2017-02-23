getK <- function(x, ids, m) {
  na.ids <- ids[!ids%in%rownames(x)]
  x.na <- matrix(0, nrow=length(na.ids), ncol=ncol(x), dimnames = list(na.ids, colnames(x)))
  X <- rbind(x, x.na)
  X <- X[ids,]
  K <- matrix(0, nrow=m, ncol=m)
  rownames(K) <- ids
  diag(K) <- 1
  diag(K)[which(rownames(K)%in%na.ids)] <- 0
  ans <- list(X=as.matrix(X), K=K)
  ans
} 
