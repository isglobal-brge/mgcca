# input x (list of matrices)
# NOTE: missing values of some variables should be imputed (to be investigated)

mgcca <- function(x, mc.cores=1, ...) {
  n<-length(x) # number of tables
  
  for (i in 1:n) assign(paste0("X", i), x[[i]])
  
  rn <- NULL
  for (i in 1:n){
  rn<- c(rn, rownames(get(paste0("X", i)))) 
  }
  rn<-sort(unique(rn))
  m<-length(rn)  # get the maximum number of individuals 
  
  X <- K <- vector("list", n)
  for (i in 1:n){
    temp <- getK(get(paste0("X", i)), ids=rn, m=m)
    X[[i]] <- as.matrix(temp$X)
    K[[i]] <- temp$K
  }

  p <- sapply(X, ncol) # number of variables per table
  numvars <- min(p) # minimum number of variables
  
  ans <- list(X=X, K=K)
  ans
}

