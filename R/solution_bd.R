solution_bd <- function(i, XX, K, XKX) {
  xkx <- XKX[[i]][[1]]
  xk <- XKX[[i]][[2]]
  # ans <- K[[i]]%*%XX[[i]]%*%xkx%*%xk

  res1 <- bdwXw( XX[[i]], diag(K[[i]]), 'wX')
  ans <- blockmult(blockmult(res1,xkx),xk)
  #..# ans <- mult_wX(XX[[i]], diag(K[[i]]))%*%xkx%*%xk

  ans
}

