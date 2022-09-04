solution_hdf5 <- function(filename, XX, XKX, mc.cores=1 ) {

    # XKX[[i]][[1]] : xkx
    # XKX[[i]][[2]] : xk

    # res1 <- bdwproduct(XX[[i]], diag(K[[i]]),"wX") # ==> No aplicat,
    #                                         es correspon directament amb XX[[i]]
    #         -> Falta testejar-ho amb més casos però de moment és la operació
    #         principal que realitzem per calcular X inicial al fitxer hdf5

    #. BD .#  ans <- bdblockmult(bdblockmult(res1,XKX[[i]][[1]], onmemory=T),XKX[[i]][[2]], onmemory = T)
    #
    #. R-functions .# ans <- mult_wX(XX[[i]], diag(K[[i]]))%*%xkx%*%xk

    # Passos a realitzar :
    #         Tenint en compte : res1 = X
    #     res_prov <- X %*% XKX
    #     ans <- res_prov %*%
    #
    #
    # QUE HE DE CALCULAR REALMENT ??
    #         X %*% XKX %*% t(X)

    bdapply_Function_hdf5(filename = filename,
                        group = "X", datasets = XX,
                        b_group = "XKX", b_datasets = XKX,
                        outgroup = "sol",func = "blockmult",
                        transp_dataset = T,transp_bdataset = T,
                        force = T)
    browser()
    solution <- bdgetDatasetsList_hdf5(filename = filename, group = "sol")
    print("Llistat obtingut")
    bdapply_Function_hdf5(filename = filename,
                        group = "sol", datasets = solution,
                        b_group = "X", b_datasets = XX,
                        outgroup = "solans",func = "blockmult",
                        force = T)
    print("Segon blockmult calculat")
    browser()
  # return(ans)
}

