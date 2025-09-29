getXKX_hdf5 <- function(filename, XX, K, inv, lambda, scores, mc.cores=1) {

    # Get XK and store (to be used later if scores == true)
    if( scores ) {
        bdapply_Function_hdf5(filename = filename,
                              group = "K",datasets = K,
                              b_group = "X", b_datasets = XX,
                              outgroup = "KX",func = "blockmult",
                              transp_dataset = T,transp_bdataset = T,
                              overwrite = TRUE)
    }

print("Debug 7.1")

    bdapply_Function_hdf5(filename = filename,
                          group = "X",datasets = XX,
                          b_group = "K", b_datasets = K,
                          outgroup = "XK",func = "blockmult",
                          # transp_dataset = T,transp_bdataset = T,
                          overwrite = TRUE)
    # browser()
print("Debug 7.2")
    bdapply_Function_hdf5(filename = filename,
                          group = "X",datasets = XX,
                          outgroup = "M",func = "tCrossProd",
                          overwrite = TRUE)

print("Debug 7.3")
    M <- bdgetDatasetsList_hdf5(filename = filename, group = "M")

    if (inv==1) { # solve
print("Debug 7.4.1")
        bdapply_Function_hdf5(filename = filename,
                              group = "M",datasets = M,
                              outgroup = "XKX",func = "invChol",
                              overwrite = TRUE, fullMatrix = T)

    } else if (inv==2) { # penalized

print("Debug 7.4.2")
        sapply(1:length(M), function(i) {

print("Debug  7.4.3 i 7.4.4")
            # Update Diagonal to compute cholesky decomposition
                bdDiag_scalar_hdf5( filename = filename,
                                 group = "M", dataset = M[i],
                                 scalar = lambda[i], operation = "+",
                                 target = "input")

                } )
print("Debug 7.5")
        # bdapply_Function_hdf5(filename = filename,
        #                     group = "tmp", datasets = M,
        #                     outgroup = "XKX", func = "invChol",
        #                     overwrite = TRUE, fullMatrix = TRUE)



        bdapply_Function_hdf5(filename = filename,
                      group = "M", datasets = M,
                      outgroup = "XKX", func = "invChol",
                      overwrite = TRUE, fullMatrix = TRUE)

        # Restore Diagonal after cholesky decomposition
browser()
print("Debug  7.6")
        sapply(1:length(M), function(i) {
            bdDiag_scalar_hdf5( filename = filename,
                                group = "M", dataset = M[i],
                                scalar = lambda[i], operation = "-",
                                target = "input")
        } )

print("Debug 7.5 --> ACABAT !!!!")
    } else {
        stop("need correct method")
    }

print("Finalitzem 7")

}
