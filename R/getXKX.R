getXKX <- function(XX, K, inv, lambda, mc.cores=1){
  ans <- mclapply(1:length(XX), getXKX.i, XX=XX, K=K, inv=inv,
                  lambda=lambda, mc.cores=mc.cores)
  ans
}

getXKX.i <- function(i, XX, K, inv, lambda){
  xk <- crossprod(XX[[i]], K[[i]])
  M <- xk%*%XX[[i]]
  if (inv=="solve")
    xkx <- solve(M)
  else if(inv=="geninv")
    xkx <- geninv(M)
  else if(inv=="penalized")
    xkx <- cgls(M, diag(ncol(M)), lambda=lambda[i])$x
  else if(inv=="ginv")
    xkx <- ginv(M)
  else
    stop("need correct method")

  ans <- list(xkx=xkx, xk=xk)
  ans
}
