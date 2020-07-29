getWeights <- function(i, A, XX, K){

  KXA <- bblockmult(blockmult(K[[i]],XX[[i]]),A[[i]])
  vv <- 1/apply(KXA, 2, sd)
  vv <-  blockmult(cbind(rep(1,nrow(A[[i]]))),rbind(vv) )
  As <- blockmult(A[[i]],vv)
  rownames(As) <- colnames(XX[[i]])
  As

}
