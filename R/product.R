productYKX0 <- function(i, Y, K, XX) {
  ans <- pinv(t(Y)%*%K[[i]]%*%Y)%*%t(Y)%*%K[[i]]%*%XX[[i]]
  ans
}

productXKY0 <- function(i, XX, K, Y) {
  ans <- pinv(t(XX[[i]])%*%K[[i]]%*%XX[[i]])%*%t(XX[[i]])%*%K[[i]]%*%Y
  ans
}


productYKX <- function(i, Y, K, XX) {
  yk <- crossprod(Y, K[[i]])
  yky <- ginv(yk%*%Y)
  ykx <- yk%*%XX[[i]]
  ans <- yky%*%ykx
  ans
}


productXKY <- function(i, XX, K, Y) {
   xk <- crossprod(XX[[i]], K[[i]])
   xkx <- ginv(xk%*%XX[[i]])
   xky <- xk%*%Y
   ans <- xkx%*%xky
   ans
  }