# input x (list of matrices - NOTE: each matrix should have the ids in the rownames 
# NOTE: missing values of some variables should be imputed (to be investigated)

mgcca <- function(x, mc.cores=1, ...) {
  n<-length(x) # number of tables
  
  if (!is.list(x))
    stop("x must be a list containing the different matrices")
  
  rn <- NULL
  for (i in 1:n){
   rn<- c(rn, rownames(x[[i]])) 
  }
  rn<-sort(unique(rn))
  m<-length(rn)  # get the maximum number of individuals 
  
  X <- K <- vector("list", n)
  for (i in 1:n){
    temp <- getK(x[[i]], ids=rn, m=m)
    X[[i]] <- as.matrix(temp$X)
    K[[i]] <- temp$K
  }

  p <- sapply(X, ncol) # number of variables per table
  numvars <- min(p) # minimum number of variables
  
  ans <- list(X=X, K=K)
  ans
}

