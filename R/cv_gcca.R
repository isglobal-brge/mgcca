cv_gcca <- function(X1, X2, grid1, grid2){
  if (missing(grid1))
    grid1 <- seq(0.001, 1, length = 5)
  if (missing(grid2))
    grid2 <- seq(0.001, 1, length = 5)

  grid <- expand.grid(grid1, grid2)
  res <- apply(grid, 1, function(x) {
    loop(X1, X2, x[1], x[2])})
  opt <- which(res==max(res), arr.ind=TRUE)
  ans <- c(grid1[opt[1]], grid1[opt[2]])
  ans
}


loop <- function(X1, X2, l1, l2){
  ans <- rep(NA, nrow(X1))
  for (i in 1:nrow(X1))
    ans[i] <- mean(mgcca(list(X1[-i,], X2[-i,]),
                         pval = FALSE,
                         method="penalized",
                         lambda=c(l1, l2))$AVE$AVE_outer_model)
  ans
}


