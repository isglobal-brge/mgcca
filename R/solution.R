solution0 <- function(i, XX, K) {
  ans <- K[[i]]%*%XX[[i]]%*%pinv(t(XX[[i]])%*%K[[i]]%*%XX[[i]])%*%t(XX[[i]])%*%K[[i]]
  ans
}


solution <- function(i, XX, K) {
  xk <- crossprod(XX[[i]],K[[i]])
  xkx <- pinv(xk%*%XX[[i]])
  ans <- K[[i]]%*%XX[[i]]%*%xkx%*%xk
  ans
}
