getXKX_hdf5 <- function(filename, XX, K, inv, lambda, mc.cores=1){

    bdapply_Function_hdf5(filename = filename,
                          group = "X",datasets = XX,
                          outgroup = "M",func = "tCrossProd",
                          force = T)

    M <- bdgetDatasetsList_hdf5(filename = filename, group = "M")

    if (inv==1) { # solve
        bdapply_Function_hdf5(filename = filename,
                              group = "M",datasets = M,
                              outgroup = "XKX",func = "invChol",
                              force = T)

    } else if (inv==2) { # penalized
        tmp <- bdgetDatasetsList_hdf5(filename = filename, group = "M")

        sapply(1:length(tmp), function(i) {
                dims <- BigDataStatMeth::bdgetDim_hdf5(filename, paste0("M/",tmp[i]))
                tmpResult <- bdScalarwproduct( diag(dims[1]) , lambda[i], "wX")
                bdAdd_hdf5_matrix( tmpResult, filename, "tmp", paste0(tmp[i],"scalarx"), force = T)

                print(paste0("Sumatori de ", tmp[i]))

                bdblockSum_hdf5( filename = filename,
                               group = "M", a = tmp[i],
                               groupB = "tmp", b = paste0(tmp[i],"scalarx"),
                               outgroup = "tmp", outdataset = tmp[i],
                               block_size = 7500 )

          } )

        bdapply_Function_hdf5(filename = filename,
                            group = "tmp",datasets = tmp,
                            outgroup = "XKX",func = "invChol",
                            force = T)

    } else {
    stop("need correct method")
    }

}
