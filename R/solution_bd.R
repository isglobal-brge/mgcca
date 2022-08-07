solution_bd <- function(i, XX, K, XKX) {
  # XKX[[i]][[1]] : xkx
  # XKX[[i]][[2]] : xk

  res1 <- bdwproduct(XX[[i]], diag(K[[i]]),"wX")
  ans <- bdblockmult(bdblockmult(res1,XKX[[i]][[1]], onmemory=T),XKX[[i]][[2]], onmemory = T)

  return(ans)
}

