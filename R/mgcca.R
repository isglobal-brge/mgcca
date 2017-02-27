# input x (list of matrices - NOTE: each matrix should have the ids in the rownames 
# NOTE: missing values of some variables should be imputed (to be investigated)

mgcca <- function(x, nfac=2, mc.cores=1, ...) {
  
  n<-length(x) # number of tables
  
  if (!is.list(x))
    stop("x must be a list containing the different matrices")
  
  rn <- sort(Reduce('union', lapply(x, rownames)))
  m <- length(rn)  # get the maximum number of individuals 
  
  XK <- lapply(x, getK, ids=rn, m=m)
  X <- lapply(XK, '[[', 1)
  K <- lapply(XK, '[[', 2)
  
  p <- sapply(X, ncol) # number of variables per table
  numvars <- min(p) # minimum number of variables
  
  Mi <- lapply(1:n, solution, XX=X, K=K)
  M <- Reduce('+', Mi)
  Ksum <- Reduce('+', K)
  
  Ksum05<-Ksum
  diag(Ksum05)<-diag(Ksum05)^(-.5)
  M<-Ksum05%*%M%*%Ksum05
  
  eig<-eigen(M, symmetric=TRUE)
  Yast<-eig$vectors[,1:nfac]
  lambda<-eig$values[1:nfac]
  Y<-sqrt(n)*Ksum05%*%Yast
  
  B <- lapply(1:n, productYKX, Y=Y, XX=X, K=K)
  A <- lapply(1:n, productXKY, XX=X, K=K, Y=Y)
  
  As <- lapply(1:n, getWeights, A=A, XX=X, K=K)
  ans <- list(Yast=Yast, B=B, A=A, As=As)
  class(ans) <- "mgcca"
  ans
}

