solution_hdf5 <- function(filename, XX, XKX, mc.cores=1 ) {

    bdapply_Function_hdf5(filename = filename,
                        group = "X", datasets = XX,
                        b_group = "XKX", b_datasets = XKX,
                        outgroup = "sol",func = "blockmult",
                        transp_dataset = T,transp_bdataset = T,
                        force = T)
    tmpSol <- bdgetDatasetsList_hdf5(filename = filename, group = "sol")

    bdapply_Function_hdf5(filename = filename,
                        group = "sol", datasets = tmpSol,
                        b_group = "X", b_datasets = XX,
                        outgroup = "solans",func = "blockmult",
                        force = T)
    datasetSol <- bdgetDatasetsList_hdf5(filename = filename, group = "solans")

    return("solans")
    # return( list( group = "solans",
    #               datasets = datasetSol) )
}
