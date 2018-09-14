getSol <- function(XX, inv, lambda, mc.cores=1){
  ans <- mclapply(1:length(XX), getSol.i, XX=XX, inv=inv,
                  lambda=lambda, mc.cores=mc.cores)
  ans
}

getSol.i <- function(i, XX, inv, lambda){
  M <- crossprod(XX[[i]])
  if (inv==1) # solve
    Minv <- chol2inv(chol(M))
  else if(inv==2) # penalized
    Minv <- chol2inv(chol(M + diag(nrow(M))*lambda[i]))
  else if(inv==3)
    Minv <- geninv(M)
  else if(inv==4)
    Minv <- ginv(M)
  else
    stop("need correct method")

  ans <- XX[[i]]%*%Minv%*%t(XX[[i]])
  ans
}

