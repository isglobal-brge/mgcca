getCor_hdf5 <- function(filename, Xgroup, x, Ygroup, y, byblocks, threads) {

    maxelements <- 30000

    common <- BigDataStatMeth::bdgetDiagonal_hdf5(filename, group = "K", dataset =  x)
    common_elems <- which(common!=0)

    if (!inherits(y, "matrix") && !inherits(y, "data.frame")) {
        if(inherits(y,"character")){
            Y <- rhdf5::h5read(filename,y, index = list(common_elems, NULL))
        }
    } else {
        if(inherits(y, "data.frame")) {
            Y <- as.matrix(y[common_elems,])
        } else {
            Y <- y[common_elems,]
        }
    }

    if (byblocks) {
            dims <- BigDataStatMeth::bdgetDim_hdf5(filename, paste0(Xgroup,"/", x))

            if(dims[1] < maxelements) { sizeblock <- maxelements / dims[1]
            } else { sizeblock <- 1 }

            xg <- split(1:dims[1], ceiling(seq_along(1:dims[1])/sizeblock))

            res <- lapply(xg, function(dat, mY) {
                X <- rhdf5::h5read(filename, paste0(Xgroup,"/",x), index = list(dat, common_elems))
                cor(t(X),mY)
            }, mY = Y)
            res <- do.call(rbind,res)

    } else {
        X <- rhdf5::h5read(filename, paste0(Xgroup,"/",x),
                           index = list(NULL, common_elems))
        res <- cor(X,Y)
    }

    rownames(res) <- t(getDimNames_hdf5(filename, Xgroup, x)$colnames)
    return(res)
}

