productYKX <- function(i, Y, K, XX) {
  yk <- crossprod(Y, K[[i]])
  yky <- ginv(yk%*%Y)
  ykx <- yk%*%XX[[i]]
  ans <- yky%*%ykx
  rownames(ans) <- colnames(Y)
  colnames(ans) <- rownames(Y)
  ans
}


productXKY <- function(i, XX, K, Y) {
   xk <- crossprod(XX[[i]], K[[i]])
   xkx <- ginv(xk%*%XX[[i]])
   xky <- xk%*%Y
   ans <- xkx%*%xky
   rownames(ans) <- colnames(XX)
   colnames(ans) <- rownames(XX)
   ans
  }