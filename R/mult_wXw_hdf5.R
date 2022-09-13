mult_wXw_hdf5 <- function(filename, groupX, X, groupw, w){
  # X matrix
  # w diag of a matrix

    bdcomputeMatrixVector_hdf5(filename, group = groupX, dataset = X,
                               vectorgroup = groupw, vectordataset = w,
                               outgroup = "tmp", outdataset = "MKsum05",
                               func = "*",
                               byrows = T,force = T)

    bdcomputeMatrixVector_hdf5(filename, group = "tmp", dataset = "MKsum05",
                               vectorgroup = groupw, vectordataset = w,
                               outgroup = "FinalRes", outdataset = "MKsum05",
                               func = "*",
                               byrows = F,force = T)
}

