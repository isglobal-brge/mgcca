productXKY_bd <- function(i, Y, XKX) {
  xkx <- XKX[[i]][[1]]
  xk <- XKX[[i]][[2]]
  xky <- blockmult(xk,y)
  ans <- blockmult(xkx,xky)
  ans
}

productYKX_bd <- function(i, Y, K, XX) {
  # not called so far ...  think about how to make comparable with XKX
  yk <- crossprod(Y, K[[i]])
  yky <- geninv(yk%*%Y) # improve
  ykx <- blockmult(yk,XX[[i]])
  ans <- blockmult(yky,ykx)
  ans
}
