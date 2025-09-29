getK_hdf5 <- function(x, ids, m, filename, group, ngroup) {

    outgroup = paste0(group, "_out")
    dimnames <- getDimNames_hdf5(filename, group, x)
    na.ids <- ids[ !ids %in% t(dimnames$rownames)]

    # if( !bdExists_hdf5_element(filename, "tmp_mgga") ) {
    # BigDataStatMeth::bdCreate_hdf5_group( filename, "tmp_mgga" )
    # }

    if( length(na.ids) == 0 ) {

        X <- cbind( dimnames[[1]],
                    order = seq(1:(nrow(dimnames[[1]]) + length(na.ids))),
                    diag = c(rep(1,nrow(dimnames[[1]]))) )

        rownames(X) <- X[,1]
        # X <- X[ids,]
        X <- X[ids, -1]
        X$newOrder = seq( 1:nrow(X))
        blocks <- split(X, cumsum(c(TRUE, diff(X$order) != 1)))

        BigDataStatMeth::bdCreate_hdf5_emptyDataset(
              filename = filename, group = "X", dataset = x,
              nrows = nrow(dimnames[[2]]), ncols = nrow(dimnames[[1]]),
              overwriteDataset = TRUE )

        bdSort_hdf5_dataset( filename = filename, group = ngroup, dataset = x,
                             outdataset = x, blockedSortlist = blocks,
                             func = "sortRows", outgroup = "X", overwrite = TRUE )

        # BigDataStatMeth::bdWriteDimnames_hdf5( filename, group = "X" ,
        #                                        dataset = x,
        #                                        rownames = t(as.vector(dimnames[[2]])),
        #                                        colnames = X$chr)
        BigDataStatMeth::bdWrite_hdf5_dimnames( filename, group = "X",
                                        dataset = x,
                                        rownames = dimnames[[2]][,1],
                                        colnames = rownames(X))
        # BigDataStatMeth::bdCreate_hdf5_emptyDataset(filename, "K", x, nrow(X), nrow(X) )
        BigDataStatMeth::bdCreate_hdf5_emptyDataset(
            filename = filename, group = "K", dataset = x,
            nrows = nrow(X), ncols = nrow(X), overwriteDataset = TRUE)
        BigDataStatMeth::bdWriteDiagonal_hdf5( rep(1, nrow(X)), filename, "K", x)
    } else {

        X <- cbind( rbind(dimnames[[1]], na.ids),
                    order = seq(1:(nrow(dimnames[[1]]) + length(na.ids))),
                    diag = c(rep(1,nrow(dimnames[[1]])), rep(0, length(na.ids))) )
        rownames(X) <- X[,1]
        # X <- X[ ids, ]
        X <- X[ ids, -1]

        X$newOrder = seq( 1:nrow(X))
        blocks <- split(X, cumsum(c(TRUE, diff(X$order) != 1)))

        diagonal <- sapply(1:length(blocks), function(x) {
            # bdAdd_hdf5_matrix(blocks[x], filename = filename, group = ngroup,
            #                   dataset = names(blocks[x]))
            # browser()
            bdCreate_hdf5_matrix(object = blocks[[x]], filename = filename, group = ngroup,
                              dataset = names(blocks[x]), overwriteDataset = TRUE)

            if(exists("diagonal")) {
                vect <- c(diagonal, blocks[[x]]$diag)
            } else {
                vect <- blocks[[x]]$diag }
                return(vect)
        }, simplify = F)

        diagonal <- do.call(c, diagonal)

        BigDataStatMeth::bdCreate_hdf5_emptyDataset(
          filename = filename, group = "X", dataset = x,
          nrows = nrow(dimnames[[2]]), ncols = length(diagonal),
          overwriteDataset = TRUE )

        bdSort_hdf5_dataset( filename = filename, group = ngroup, dataset = x,
                             outdataset = x, blockedSortlist = blocks,
                             func = "sortRows", outgroup = "X", overwrite = TRUE )
        # BigDataStatMeth::bdWriteDimnames_hdf5( filename, group = "X" ,
        #                                        dataset = x,
        #                                        rownames = t(as.vector(dimnames[[2]])),
        #                                        colnames = X$chr)

        BigDataStatMeth::bdWrite_hdf5_dimnames( filename, group = "X" ,
                                       dataset = x,
                                       rownames = dimnames[[2]][,1],
                                       colnames = rownames(X))
        # BigDataStatMeth::bdCreate_hdf5_emptyDataset(filename, "K", x, nrow(X), nrow(X) )
        BigDataStatMeth::bdCreate_hdf5_emptyDataset( filename = filename,
                    group = "K", dataset = x, nrows = nrow(X), ncols = nrow(X),
                    overwriteDataset = TRUE )

        BigDataStatMeth::bdWriteDiagonal_hdf5( diagonal = diagonal,
                               filename = filename, group = "K", dataset = x)

    }

}
