productYKX <- function(i, Y, K, XX) {
  ans <- pinv(t(Y)%*%K[[i]]%*%Y)%*%t(Y)%*%K[[i]]%*%XX[[i]]
  ans
}

productXKY <- function(i, XX, K, Y) {
  ans <- pinv(t(XX[[i]])%*%K[[i]]%*%XX[[i]])%*%t(XX[[i]])%*%K[[i]]%*%Y
  ans
}