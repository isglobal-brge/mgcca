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

  a<<-X
  aa<<-K
  
  p <- sapply(X, ncol) # number of variables per table
  numvars <- min(p) # minimum number of variables
  
  Mi <- lapply(1:n, solution, XX=X, K=K)
  M <- Reduce('+', Mi)
  Ksum <- Reduce('+', K)
  
  Ksum05<-Ksum
  diag(Ksum05)<-diag(Ksum05)^(-.5)
  M<-Ksum05%*%M%*%Ksum05
  
  eig<-eigen(M)
  Yast<-eig$vectors
  lambda<-eig$values
  Y<-sqrt(n)*Ksum05%*%Yast
  
  B <- lapply(1:n, productYKX, Y=Y, XX=X, K=K)
  A <- lapply(1:n, productXKY, XX=X, K=K, Y=Y)
  
  As <- lapply(1:n, getWeights, A=A, XX=X, K=K)
  As
}

