productXKY <- function(i, Y, XKX, inv="solve") {
  xkx <- XKX[[i]][[1]]
  xk <- XKX[[i]][[2]]
  xky <- xk%*%Y
  ans <- xkx%*%xky
  ans
}

productYKX <- function(i, Y, K, XX) {
  # not called so far ...  think about how to make comparable with XKX
  yk <- crossprod(Y, K[[i]])
  yky <- geninv(yk%*%Y) # improve
  ykx <- yk%*%XX[[i]]
  ans <- yky%*%ykx
  ans
}


