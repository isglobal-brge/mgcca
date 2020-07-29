getXKX_bd <- function(XX, K, inv, lambda, mc.cores=1){
  ans <- parallel::mclapply(1:length(XX), getXKX_bd.i, XX=XX, K=K, inv=inv,
                  lambda=lambda, mc.cores = mc.cores)
  ans

}

getXKX_bd.i <- function(i, XX, K, inv, lambda)
{

  xk <- BigDataStatMeth::blockmult(t(XX[[i]]),K[[i]])
  M <- BigDataStatMeth::blockmult(xk,XX[[i]])

  if (inv==1) # solve
    xkx <- bdInvCholesky(M)
  else if (inv==2) # penalized
    xkx <- bdInvCholesky(M + bdwproduct(diag(nrow(M)), lambda[i], "wX"))

  else
    stop("need correct method")

  ans <- list(xkx=xkx, xk=xk)
  ans
}
