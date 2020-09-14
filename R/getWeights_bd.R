getWeights <- function(i, A, XX, K){

  KXA <- bblockmult(blockmult(K[[i]],XX[[i]]),A[[i]], onmemory = T)
  vv <- 1/apply(KXA, 2, sd)
  vv <-  blockmult(cbind(rep(1,nrow(A[[i]]))),rbind(vv) , onmemory = T)
  As <- blockmult(A[[i]],vv, onmemory = T)
  rownames(As) <- colnames(XX[[i]])
  As

}
