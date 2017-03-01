ginv2 <- function (X, tol = sqrt(.Machine$double.eps))
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- irlba(X, nv=5)
  Xsvd$v <- Xsvd$u
  if (is.complex(X))
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) {
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  }
  else if (!any(Positive))
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}


# Returns the Moore-Penrose inverse of the argument"
geninv <- function(G){
  # Transpose if m < n"
  m <- nrow(G)
  n <- ncol(G)
  transpose <- FALSE
   if (m<n){
    transpose <- TRUE
    A <- tcrossprod(G)
    n <- m
   }
   else {
     A <- crossprod(G)
   }
   # Full rank Cholesky factorization of A
   dA <- diag(A)
   tol <- min(dA[dA>0])*1e-9
   L <- matrix(0, nrow=nrow(A), ncol=ncol(A))
   r <- 0
   for (k in 1:n){
     r <- r+1
     L[k:n,r] <- A[(k:n),k,drop=FALSE] - tcrossprod(L[(k:n),(1:(r-1))], L[k,(1:(r-1)), drop=FALSE])
     # Note: for r=1, the substracted vector is zero"
     if (L[k,r]>tol){
       L[k,r] <- sqrt(L[k,r])
       if (k<n){
         L[(k+1):n,r] <- L[(k+1):n,r]/L[k,r]
       }
     }
     else {
       r <- r-1
     }
   }
    L <- L[,1:r]
    # Computation of the generalized inverse of G"
    M <- solve(crossprod(L))
    if (transpose){
       Y <- crossprod(G,L)%*%M%*%tcrossprod(M,L)
    }
    else {
     Y <- L%*%M%*%tcrossprod(tcrossprod(M,L), G)
    }
    Y
}
