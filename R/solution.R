solution <- function(i, XX, K) {
  xk <- crossprod(XX[[i]],K[[i]])
  xkx <- MASS::ginv(xk%*%XX[[i]])
  ans <- K[[i]]%*%XX[[i]]%*%xkx%*%xk
  ans
}
