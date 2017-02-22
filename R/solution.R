solution <- function(i, XX, K) {
  ans <- K[[i]]%*%XX[[i]]%*%pinv(t(XX[[i]])%*%K[[i]]%*%XX[[i]])%*%t(XX[[i]])%*%K[[i]]
  ans
}
