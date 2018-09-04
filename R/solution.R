solution <- function(i, XX, K, XKX) {
  xkx <- XKX[[i]][[1]]
  xk <- XKX[[i]][[2]]
  # ans <- K[[i]]%*%XX[[i]]%*%xkx%*%xk
  ans <- mult_wX(XX[[i]], diag(K[[i]]))%*%xkx%*%xk
  ans
}
