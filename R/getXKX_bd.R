getXKX_bd <- function(XX, K, inv, lambda, mc.cores=1){
  # ans <- parallel::mclapply(1:length(XX), getXKX_bd.i, XX=XX, K=K, inv=inv,
  #                 lambda=lambda, nthreads=mc.cores)
  ans <- lapply(1:length(XX), getXKX_bd.i, XX=XX, K=K, inv=inv,
                            lambda=lambda, nthreads=mc.cores)
  ans

}

getXKX_bd.i <- function(i, XX, K, inv, lambda, nthreads)
{

  xk <- BigDataStatMeth::bdblockmult(t(XX[[i]]),K[[i]], onmemory = T, threads = nthreads)
  M <- BigDataStatMeth::bdblockmult(xk,XX[[i]], onmemory = T, threads = nthreads)
  if (inv==1) # solve
    xkx <- bdInvCholesky(M)
  else if (inv==2) # penalized
    xkx <- bdInvCholesky(M + bdScalarwproduct(diag(nrow(M)), lambda[i], "wX"))

  else
    stop("need correct method")

  ans <- list(xkx=xkx, xk=xk)
  ans
}

