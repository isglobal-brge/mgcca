getScores0 <- function(i, Y, XX, Minv){
  minv <<- Minv[[i]]
  xy <<- crossprod(XX[[i]],Y)
  A <<- minv%*%xy
  XA <- XX[[i]]%*%A
  vv <- 1/apply(XA, 2, sd)
  vv <- cbind(rep(1,nrow(A)))%*%rbind(vv)
  As <<- A*vv
  scores <- XX[[i]]%*%As
  scores
}
