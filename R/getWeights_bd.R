getWeights_bd <- function(i, A, XX, K){

  KXA <- bdblockmult(bdblockmult(K[[i]],XX[[i]], onmemory = T),A[[i]], onmemory = T)
  vv <- 1/apply(KXA, 2, sd)
  vv <-  bdblockmult(cbind(rep(1,nrow(A[[i]]))),rbind(vv) , onmemory = T)
  As <- A[[i]]*vv
  rownames(As) <- colnames(XX[[i]])
  As

}
