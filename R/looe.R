# Y = a + bX
# A'A = X
# Solve the linear system (A'A + lambda * I)x = Ab using
# eigen decomposition
#
# example: X data Y outcome
#    M <- tcrossprod(X)
#    LOOE(M, Y)


inversecpp <- function(X, lambda=1, eigen=TRUE,
                       Lambda, Q){
  if (eigen){
   ee <- eigen(X, symmetric = TRUE)
   Lambda <- ee$values
   Lambda[Lambda<0] <- 0
   Q <- ee$vectors
  }
  else
    if(missing(Lambda) | missing(Q))
      stop('SVD results should be provided. \n')

  if (lambda == 1)
    W <-  1/Lambda
  else
    W <- 1/(Lambda + lambda)

  Ginv <- rfunctions::crossprodcpp(t(Q), W)
  # Q%*%diag(W)%*%t(Q)
  Ginv
}

solveEigen <- function(X, Y, lambda){
  XX <- tcrossprod(X)
  Ginv <- inversecpp(XX, lambda=lambda)
  coef <- t(Ginv%*%X)%*%Y
  coef
}


#
# Compute LOOE
#

LOOE <- function(X, Y, nlambdas=100, max.lambda=1, lambdas){

  if (missing(lambdas)){
    lambdas <- seq(0.001, max.lambda, length=nlambdas)
  }
  ee <- eigen(X, symmetric = TRUE)
  Lambda <- ee$values
  Lambda[Lambda<0] <- 0
  Q <- ee$vectors
  looe <- sapply(lambdas, LOOE.i, Lambda=Lambda, Q=Q, Y=Y)

  lambda.min <- lambdas[which.min(looe)]

  Ginv <- inversecpp(X, lambda.min)

  ans <- list(looe=looe, Ginv=Ginv, lambdas=lambdas,
              lambda.min=lambda.min)
  ans
}

LOOE.i <- function(lambda, Lambda, Q, Y){
  Ginv <- inversecpp(lambda=lambda, Lambda=Lambda,
                     Q=Q, eigen=FALSE)
  cte <- Ginv%*%Y
  ans <- sum((cte/diag(Ginv))^2)
  ans
}


LOOE.i <- function(lambda, Lambda, Q, Y){
  xx <- (MM%*%solve(MM + diag(lambda, nrow(MM))))%*%t(MM)
  num <- Y - xx
  out <- sum(num / diag(xx)) ^ 2
  out
}

