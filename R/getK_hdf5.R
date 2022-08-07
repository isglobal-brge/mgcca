getK_hdf5 <- function(x, ids, m, filename, group) {

    outgroup = paste0(group, "_out")
    dimnames <- getDimNames_hdf5(filename, group, x)
    na.ids <- ids[ !ids %in% t(dimnames$rownames) ]

    if( length(na.ids) == 0 ) {

        X <- cbind(dimnames[[1]], order = seq(1:(nrow(dimnames[[1]]) + length(na.ids))))
        rownames(X) <- X[,1]
        X <- X[ids,]
        X$newOrder = seq( 1:nrow(X))
        blocks <- split( X, cumsum(c(TRUE, diff(X$order) != 1)) )

        sapply( 1:length(blocks), function(x) {
          bdAdd_hdf5_matrix(blocks[x], filename = filename, group = group, dataset = names(blocks[x]))
        } )

        browser()

        BigDataStatMeth::bdBind_hdf5( filename, group, bdgetDatasetsList_hdf5(filename, group), "X", outdataset = x, func = "bindRows")
        BigDataStatMeth::bdCreateEmptyDataset_hdf5(filename, "K", x, nrow(X), nrow(X) )
        BigDataStatMeth::bdWriteDiagonal_hdf5( rep(1, nrow(X)), filename, "K", x)

    } else {

        X <- cbind( rbind(dimnames[[1]], na.ids),
                    order = seq(1:(nrow(dimnames[[1]]) + length(na.ids))),
                    diag = c(rep(1,nrow(dimnames[[1]])), rep(0, length(na.ids))) )
        rownames(X) <- X[,1]
        X <- X[ids,]
        X$newOrder = seq( 1:nrow(X))
        blocks <- split(X, cumsum(c(TRUE, diff(X$order) != 1)))

        diagonal <- sapply(1:length(blocks), function(x) {
            bdAdd_hdf5_matrix(blocks[x], filename = filename, group = group, dataset = names(blocks[x]))
            if(exists("diagonal")) {
                vect <- c(diagonal, blocks[[x]]$diag)
            } else {
                vect <- blocks[[x]]$diag }
            return(vect)
        }, simplify = F)

        diagonal <- do.call(c, diagonal)
        BigDataStatMeth::bdBind_hdf5( filename, group, bdgetDatasetsList_hdf5(filename, group), "X", outdataset = x, func = "bindRows")
        BigDataStatMeth::bdCreateEmptyDataset_hdf5(filename, "K", x, nrow(X), nrow(X) )
        BigDataStatMeth::bdWriteDiagonal_hdf5( diagonal, filename, "K", x)
    }
}
