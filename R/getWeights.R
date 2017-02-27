getWeights <- function(i, A, XX, K){
  KXA <- K[[i]]%*%XX[[i]]%*%A[[i]]
  vv <- 1/apply(KXA, 2, sd)
  vv <- cbind(rep(1,nrow(A[[i]])))%*%rbind(vv)
  As <- A[[i]]*vv
  rownames(As) <- colnames(XX)
  colnames(As) <- paste0("can", 1:ncol(As))
  As
}