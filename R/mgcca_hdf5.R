#' Generalized canonical correlation with missing individuals for big data
#'
#' @param x string for hdf5 file where data to process will be under MGCCA_IN group .
#' Results will be stored in the same data file under MGCCA_OUT group. Missing is not allowed  (see details).
#' @param filename string for hdf5 file where data to process will be under MGCCA_IN group .
#' Results will be stored in the same data file under MGCCA_OUT group. Missing is not allowed  (see details).
#' @param datanames string array with datasetnames to use with mcgga
#' @param nfac ...
#' @param scale ...
#' @param pval should p-values of correlation between variables and shared canonical variates be computed? Default is TRUE.
#' @param scores should canonical variables be computed for each table?
#'               Default is FALSE. See details
#' @param method ...
#' @param lambda ...
#' @param mc.cores ...
#' @details NOT WORKING
#' @return a list consisting of
#'   \item{Y}{canonical components for the shared space}
#'   \item{corsY}{correlation between variables and shared canonical components}
#'   \item{scores}{canonical componets for each table}
#'   \item{p.values}{p-values of correlation between variables and shared canonical components}
#'   \item{AVE}{indicators of model quality based on the Average Variance Explained (AVE): AVE for each table, AVE_outer (average accross tables), AVE_inner(for the shared component).}
#' @examples
#' see vignette
#'
#' @export
#' @importFrom parallel mclapply
#' @importFrom MASS ginv
#' @importFrom RSpectra eigs eigs_sym

# mgcca_hdf5 <- function(x, datanames, nfac=2, scale=TRUE, pval=TRUE, scores=FALSE,
#                      method="solve", lambda, mc.cores=1, ...) {

mgcca_hdf5 <- function(x, filename, group, datasets, nfac=2, scale=TRUE, pval=TRUE, scores=FALSE,
                       method="solve", lambda, mc.cores=1, ...) {

    # #' @importFrom rfunctions geninv cgls
    inv.type <- c("solve", "penalized")
    inv.method <- charmatch(method, inv.type, nomatch = 0)
    if (inv.method == 0)
        stop("method should be 'solve' or 'penalized' \n")

    if (BigDataStatMeth::bdIsLocked_hdf5(filename))
        stop("HDF5 file is in use by another process; aborting.")

    n <- length(datasets) # number of tables

    if (inv.method == 2) {
      if (missing(lambda)) {
          stop("penalized method requires lambda parameter \n")
      } else {
          if(length(lambda)!=n) {
              stop("lambda must be a vector of length equal to the number of tables \n")
          }
      }
    }

    # Create FINAL RESULTS folder
    BigDataStatMeth::bdCreate_hdf5_group(filename, "FINAL_RESULTS" )

    if (scale) {
        x <- lapply(datasets, bdNormalize_hdf5, filename = filename, group = group, bcenter=TRUE, bscale=TRUE, overwrite = TRUE, byrows = TRUE )
    }

    currentdatasets <- sapply(datasets, function(d, g) {
        if(scale) {
            daasetname <- paste("NORMALIZED", g, d, sep = "/")
        } else {
            daasetname <- paste( g, d, sep = "/")
        }
    }, g = group)

    # Get current groupname
    distancia <- regexpr("\\/[^\\/]*$", currentdatasets[1])[[1]]
    ngroup <-  substr(currentdatasets[1],1,distancia-1)

    if(length(currentdatasets) <= 1)
        stop("we need more than one dataset to perform mgcca analysis")

    ns <- sapply( currentdatasets, function(dataset, file) { return( bdgetDim_hdf5(file, dataset)[2] ) }, file = filename )

    if(max(ns)==min(ns)) { # check whether there are missing individuals
        rn <- Reduce('union', sapply(datasets, getRowNames_hdf5, filename = filename, group = group))
    } else {
        rn <- sort(Reduce('union', sapply(datasets, getRowNames_hdf5, filename = filename, group = group)))
    }

    m <- length(rn)  # get the maximum number of individuals

    mclapply( datasets, getK_hdf5, ids=rn, m=m, mc.cores=mc.cores, filename = filename, group = group, ngroup = ngroup )

    X <- bdgetDatasetsList_hdf5(filename = filename, group = "X")
    K <- bdgetDatasetsList_hdf5(filename = filename, group = "K")

    p <-  sapply( paste0( "X/",X), function(el, filename) {
        res <- BigDataStatMeth::bdgetDim_hdf5(filename, el)
        return(res[2]) # number of variables per table
    }, filename = filename )

    numvars <- min(p) # minimum number of variables
    getXKX_hdf5(filename, X, K, inv.method, lambda=lambda, scores, mc.cores=mc.cores)

    XKX <- bdgetDatasetsList_hdf5(filename = filename, group = "XKX")

print("Mi")
    Mi <- solution_hdf5( filename = filename, X = X, XKX = XKX, mc.cores)

    bdReduce_hdf5_dataset(filename = filename, group = Mi,
              reducefunction = "+", outgroup = "FinalRes", outdataset = "M")
    bdReduce_hdf5_dataset(filename = filename, group = "K",
              reducefunction = "+", outgroup = "FinalRes", outdataset = "Ksum")

print("Debug previ ksum05")

    ksum05 <- bdDiag_scalar_hdf5(filename = filename, group = "FinalRes", dataset = "Ksum",
                                 scalar = -0.5,operation = "pow", target = "new",
                                 outgroup = "FinalRes", outdataset = "Ksum05")

    mult_wXw_hdf5(filename, "FinalRes","M", ksum05$gr, ksum05$ds)

    eigs <- bdEigen_hdf5( filename = filename, group = "FinalRes", dataset = "MKsum05",
                  bcenter = F, bscale = F, k = 2)

print("Debug previ Yast")

    Yast <- Re( rhdf5::h5read(filename, eigs$vectors))[,1:nfac]

    Y <- sqrt(n)* bd_wproduct(Yast,
                      rhdf5::h5read(filename, paste0(ksum05$gr, "/", ksum05$ds)),
                      "wX")
    colnames(Y) <- paste0("comp", 1:ncol(Y))
    rownames(Y) <- rn

    bdCreate_hdf5_matrix(object = Y, filename = filename, group = "FINAL_RESULTS",
                         dataset = "Y", overwriteDataset = TRUE)

    if (scores) {
        Yd <- bdgetDatasetsList_hdf5(filename = filename, group = "FinalRes/Y")
        KX <- bdgetDatasetsList_hdf5(filename = filename, group = "KX")
        XK <- bdgetDatasetsList_hdf5(filename = filename, group = "XK")
        A <- productXKY_hdf5(filename=filename, Y=Yd, XK=XK, XKX=XKX, mc.cores)
        As <- getWeights_hdf5(filename, A=A, XX=X, K=K, KX=KX, initialGroup = group, mc.cores)
        scores <- getScores_hdf5(filename, XX = X, As = As, mc.cores)
        sapply(scores, function(sdataset) {
            # bdWriteDimnames_hdf5( filename, group = "FINAL_RESULTS/scores", dataset = sdataset,rownames = rn,colnames = colnames(Y))
            bdWrite_hdf5_dimnames( filename, group = "FINAL_RESULTS/scores",
                       dataset = sdataset,rownames = rn,colnames = colnames(Y))
        } )
    } else {
        scores <- NULL
    }

    lapply( X, function(dataX) {
        getCor_hdf5(filename = filename, x = dataX, Xgroup = "X", Ygroup = "FINAL_RESULTS",
                           y = "Y", byblocks = T, threads = mc.cores)
    })

    dst_corsy <- bdgetDatasetsList_hdf5( filename = filename,
                                         group = "/FINAL_RESULTS/corsY")

    get_AVE(filename = filename, group = "/FINAL_RESULTS/corsY",
            datasets = dst_corsy, nfac = nfac, p = p)

}
