getCor_hdf5 <- function(filename, x, y, byblocks, threads) {

    maxelements <- 300

    if (byblocks) {
            dims <- BigDataStatMeth::bdgetDim_hdf5(filename,x)

            if (!inherits(y, "matrix") && !inherits(y, "data.frame")) {
                if(inherits(Y,"character")){
                    Y <- rhdf5::h5read(filename,y)
                }
            } else {
                if(inherits(y, "data.frame")) {
                    Y <- as.matrix(y)
                } else {
                    Y <- y
                }
            }

            if(dims[1] < maxelements) { sizeblock <- maxelements / dims[1]
            } else { sizeblock <- 1 }

            xg <- split(1:dims[1], ceiling(seq_along(1:dims[1])/sizeblock))

            res <- lapply(xg, function(dat, mY) {
                X <- rhdf5::h5read(filename, x, index = list(dat, NULL))
                cor(t(X),mY)
            }, mY = Y)
            return(do.call(rbind,res))

    } else {
        return(cor(X,Y))
    }

}

