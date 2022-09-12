productXKY_hdf5 <- function(filename, Y, XK, XKX, threads) {
    # xkx <- XKX[[i]][[1]]
    # xk <- XKX[[i]][[2]]

    # xky <- blockmult(xk,y, onmemory = T)
    # ans <- blockmult(xkx,xky, onmemory = T)

    bdapply_Function_hdf5(filename = filename,
                          group = "XK", datasets = XK,
                          b_group = "FinalRes/Y", b_datasets = Y,
                          outgroup = "XKY",func = "blockmult",
                          force = T)
    XKY <- bdgetDatasetsList_hdf5(filename = filename, group = "XKY")

    browser()

    bdapply_Function_hdf5(filename = filename,
                          group = "XKX",datasets = XKX,
                          b_group = "XKY", b_datasets = XKY,
                          outgroup = "scores/A",func = "blockmult",
                          force = T)
    A <- bdgetDatasetsList_hdf5(filename = filename, group = "scores/A")

    return(A)
}

# productYKX_bd <- function(i, Y, K, XX) {
#   # not called so far ...  think about how to make comparable with XKX
#   yk <- crossprod(Y, K[[i]])
#   yky <- geninv(yk%*%Y) # improve
#   ykx <- blockmult(yk,XX[[i]], onmemory = T)
#   ans <- blockmult(yky,ykx, onmemory = T)
#   ans
# }
