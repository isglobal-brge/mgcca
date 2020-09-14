getScores0 <- function(i, Y, XX, Minv){
  minv <<- Minv[[i]]
  # xy <<- crossprod(XX[[i]],Y)
  blockmult(t(XX[[i]]),Y, onmemory = T)
  # A <<- minv%*%xy
  A <- blockmult(minv,xy, onmemory = T)

  # XA <- XX[[i]]%*%A
  XA <- blockmult(XX[[i]],A, onmemory = T)
  vv <- 1/apply(XA, 2, sd)
  # vv <- cbind(rep(1,nrow(A)))%*%rbind(vv)
  vv <- blockmult(cbind(rep(1,nrow(A))),rbind(vv), onmemory = T)
  # As <<- A*vv
  As <<- blockmult(A,vv, onmemory = T)
  # scores <- XX[[i]]%*%As
  scores <- blockmult(XX[[i]],As, onmemory = T)
  scores
}
