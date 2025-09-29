mult_wXw_hdf5 <- function(filename, groupX, X, groupw, w){
  # X matrix
  # w diag of a matrix
    res <- bdcomputeMatrixVector_hdf5(filename, group = groupX, dataset = X,
                                      vectorgroup = groupw, vectordataset = w,
                                      outgroup = "tmp", outdataset = "MKsum05",
                                      func = "*",
                                      byrows = F,overwrite = TRUE)

    bdcomputeMatrixVector_hdf5(filename, group = res$gr, dataset = res$ds,
                               vectorgroup = groupw, vectordataset = w,
                               outgroup = "FinalRes", outdataset = res$ds,
                               func = "*",
                               byrows = F,overwrite = TRUE)
}

