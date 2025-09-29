getK_bd <- function(x, ids, m) {
  na.ids <- ids[!ids%in%rownames(x)]
  x.na <- matrix(0, nrow=length(na.ids), ncol=ncol(x), dimnames = list(na.ids, colnames(x)))
  X <- rbind(x, x.na)
  X <- X[ids,]
  K <- diag(1, nrow = m, ncol = m)
  # K <- matrix(0, nrow=m, ncol=m)
  rownames(K) <- ids
  # diag(K) <- 1
  diag(K)[which(rownames(K)%in%na.ids)] <- 0
  ## !!! needed ¿? !!! ## rownames(X) <- ids
  ## !!! needed ¿? !!! ## colnames(X) <- colnames(x)
  ans <- list( X = as.matrix(X),
               K = K)
  return(ans)
}
