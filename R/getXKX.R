getXKX <- function(XX, K, inv, lambda, mc.cores=1){
  ans <- mclapply(1:length(XX), getXKX.i, XX=XX, K=K, inv=inv,
                  lambda=lambda, mc.cores=mc.cores)
  ans
}

getXKX.i <- function(i, XX, K, inv, lambda){
  # xk <- crossprod(XX[[i]], K[[i]])
  xk <- mult_Xw(t(XX[[i]]), diag(K[[i]]))
  M <- xk%*%XX[[i]]
  if (inv==1) # solve
    xkx <- chol2inv(chol(M))
  else if (inv==2) # penalized
    xkx <- chol2inv(chol(M + diag(nrow(M))*lambda[i]))
  else
    stop("need correct method")

  ans <- list(xkx=xkx, xk=xk)
  ans
}
