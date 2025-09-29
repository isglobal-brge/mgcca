solution_hdf5 <- function(filename, XX, XKX, mc.cores=1 ) {

    bdapply_Function_hdf5(filename = filename,
                          group = "X", datasets = XX,
                          b_group = "XKX", b_datasets = XKX,
                          outgroup = "sol",func = "blockmult",
                          transp_dataset = T,transp_bdataset = T,
                          overwrite = TRUE)

    tmpSol <- bdgetDatasetsList_hdf5(filename = filename, group = "sol")

    bdapply_Function_hdf5(filename = filename,
                        group = "sol", datasets = tmpSol,
                        b_group = "X", b_datasets = XX,
                        outgroup = "Mi",func = "blockmult",
                        overwrite = TRUE)
    # datasetSol <- bdgetDatasetsList_hdf5(filename = filename, group = "Mi")
print("OK -  solution_df5 acabat!!")
    return("Mi")
    # return( list( group = "solans",
    #               datasets = datasetSol) )
}
