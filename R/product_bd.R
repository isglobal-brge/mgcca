productXKY_bd <- function(i, Y, XKX) {
  xkx <- XKX[[i]][[1]]
  xk <- XKX[[i]][[2]]
  xky <- blockmult(xk,y, onmemory = T)
  ans <- blockmult(xkx,xky, onmemory = T)
  ans
}

productYKX_bd <- function(i, Y, K, XX) {
  # not called so far ...  think about how to make comparable with XKX
  yk <- crossprod(Y, K[[i]])
  yky <- geninv(yk%*%Y) # improve
  ykx <- blockmult(yk,XX[[i]], onmemory = T)
  ans <- blockmult(yky,ykx, onmemory = T)
  ans
}
