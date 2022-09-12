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
#' @importFrom rfunctions geninv cgls
#' @importFrom MASS ginv
#' @importFrom RSpectra eigs eigs_sym

# mgcca_hdf5 <- function(x, datanames, nfac=2, scale=TRUE, pval=TRUE, scores=FALSE,
#                      method="solve", lambda, mc.cores=1, ...) {

mgcca_hdf5 <- function(x, filename, group, datasets, nfac=2, scale=TRUE, pval=TRUE, scores=FALSE,
                       method="solve", lambda, mc.cores=1, ...) {

  inv.type <- c("solve", "penalized")
  inv.method <- charmatch(method, inv.type, nomatch = 0)
  if (inv.method == 0)
    stop("method should be 'solve' or 'penalized' \n")

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

  ##.. HDF5 can't sotre <NA> --> We never foud NA in it. # nas <- sapply(x, function(x) any(is.na(x)))

  ##.. if (any(nas))
  ##..   stop("Missing values are not allowed. Either use 'impute' package or
  ##..       use tables with complete cases.")

  if (scale) {
      #..# x <- lapply(x, scale)
      x <- lapply(datasets, bdNormalize_hdf5, filename = filename, group = group, bcenter=TRUE, bscale=TRUE, force = TRUE, byrows = TRUE )
  }

  currentdatasets <- sapply(datasets, function(d, g) {
      if(scale) {
          daasetname <- paste("NORMALIZED", g, d, sep = "/")
      } else {
          daasetname <- paste( g, d, sep = "/")
      }
  }, g = group)

  distancia <- regexpr("\\/[^\\/]*$", currentdatasets[1])[[1]]
  ngroup <-  substr(currentdatasets[1],1,distancia-1)
  # print("Examinar currentdatasets per extreure posteriorment el grup")
  # print(distancia)
  # print(substr(currentdatasets[1],1,distancia-1))
  # browser()
  ## FINS AQUÍ OK !!! NORMALITZA LES DADES OK !!!

  #..  No te sentit ..# if (!is.list(x))
  #..  No te sentit ..#   stop("x must be a list containing the different matrices")

  if(length(currentdatasets)<=1)
    stop("we need more than one dataset to perform mgcca analysis")

  #..  No te sentit ..# if (any(unlist(lapply(x, function(x) !is.matrix(x)))))
  #..  No te sentit ..#   x <- lapply(x, as.matrix)


  #..# ns <- sapply(x, nrow)
  ns <- sapply( currentdatasets, function(dataset, file) { return( bdgetDim_hdf5(file, dataset)[2] ) }, file = filename )


  if(max(ns)==min(ns)) { # check whether there are missing individuals
      rn <- Reduce('union', sapply(datasets, getRowNames_hdf5, filename = filename, group = group))
  } else {
      rn <- sort(Reduce('union', sapply(datasets, getRowNames_hdf5, filename = filename, group = group)))
  }

  m <- length(rn)  # get the maximum number of individuals

  ##..## mclapply( datasets, getK_hdf5, ids=rn, m=m, mc.cores=mc.cores, filename = filename, group = group )
  mclapply( datasets, getK_hdf5, ids=rn, m=m, mc.cores=mc.cores, filename = filename, group = group, ngroup = ngroup )
  X <- bdgetDatasetsList_hdf5(filename = filename, group = "X")
  K <- bdgetDatasetsList_hdf5(filename = filename, group = "K")
  p <-  sapply( paste0( "X/",X), function(el, filename){
        res <- BigDataStatMeth::bdgetDim_hdf5(filename, el)
        return(res[2]) # number of variables per table
    }, filename = filename )

  numvars <- min(p) # minimum number of variables

  getXKX_hdf5(filename, X, K, inv.method, lambda=lambda, scores, mc.cores=mc.cores)

  # # Get the required XKX product and inverse that is computed multiple times
  # XKX <- getXKX_bd(X, K, inv.method, lambda=lambda, mc.cores=mc.cores)
  #..# Mi <- mclapply(1:n, solution, XX=X, K=K, XKX=XKX, mc.cores=mc.cores)

  XKX <- bdgetDatasetsList_hdf5(filename = filename, group = "XKX")
  Mi <- solution_hdf5( filename = filename, X = X, XKX = XKX, mc.cores)

  #..# M <- Reduce('+', Mi)
  BigDataStatMeth::bdReduce_matrix_hdf5(filename = filename, group = Mi, reducefunction = "+", outgroup = "FinalRes", outdataset = "M")

  #..# Ksum <- Reduce('+', K)
  BigDataStatMeth::bdReduce_matrix_hdf5(filename = filename, group = "K", reducefunction = "+", outgroup = "FinalRes", outdataset = "Ksum")

  ## OPTIMITZACIÓ DE CODI REPECTE ALTRES VERSIONS :
  ##    Eliminada la funció de càlcul de l'arrel quadrada de la diagonal
  ##    Modificada la función wXw
  ##            --> Enlloc dels passos anteriors : multipliquem la diagonal (directament) %*% matriuX
  ##            --> Estalviant el càlcul de l'arrel quadrada i el doble producte
  ##
  ##
  ##    Codi modificaat :
  #..# Ksum05 <- diag(Ksum)^(-0.5)
  #..# MKsum05 <- mult_wXw(M, Ksum05)

  Ksum <- BigDataStatMeth::bdgetDiagonal_hdf5(filename, "FinalRes", "Ksum")
  bdAdd_hdf5_matrix(as.matrix((Ksum^-0.5)^2), filename, "FinalRes", "Ksum05")
  bdWeightedProduct_hdf5(filename, group = "FinalRes", dataset = "M",
                         vectorgroup = "FinalRes", vectordataset = "Ksum05",
                         outdataset = "MKsum05", byrows = T,force = T)

  #..# eig <- BDSM::bdSVD_lapack(MKsum05, bcenter = FALSE, bscale = FALSE)
  BigDataStatMeth::bdSVD_hdf5( filename,group = "FinalRes", dataset = "MKsum05",bcenter = F,bscale = F)

  #..# Yast <- Re(eig$u[,1:nfac])
  Yast <- Re( rhdf5::h5read(filename,"SVD/MKsum05/u",))[,1:nfac]

  ##
  ## REPASSAR TOTA AQUESTA PART !!!! NO ACABEN DE QUADRAR ELS RESULTATS
  ## DE Y !!!!
  ##

  Y <- sqrt(n)* bdwproduct(Yast, Ksum^-0.5, "wX")
  colnames(Y) <- paste0("comp", 1:ncol(Y))
  rownames(Y) <- rn

  sapply(datasets, function(dataset) {
      bdAdd_hdf5_matrix(object = Y, filename = filename, group = "FinalRes/Y",
                        dataset = paste0(dataset,".Y"), force = T)
  } )

  ##
  ##
  ## REPASSAR TOT L'ANTERIOR !!!! PERQUÈ A PARTIR D'AQUÍ ELS RESULTATS NO SON CORRECTES !!!
  ##
  ##


  if (scores) {
      Yd <- bdgetDatasetsList_hdf5(filename = filename, group = "FinalRes/Y")
      XK <- bdgetDatasetsList_hdf5(filename = filename, group = "XK")

      #..# A <- mclapply(1:n, productXKY_bd, Y=Y, XKX=XKX, mc.cores=mc.cores)
      A <- productXKY_hdf5(filename=filename, Y=Yd, XK=XK, XKX=XKX, mc.cores)

      #..# As <- mclapply(1:n, getWeights_bd, A=A, XX=X, K=K, mc.cores=mc.cores)
      As <- getWeights_hdf5(filename=filename, A=A, XX=X, K=K, mc.cores)

      #..# scores <- mclapply(1:n, getScores_bd, dat=X, As=As, mc.cores=mc.cores)
      scores <- lapply(1:n, getScores_hdf5, dat = X, As = As, mc.cores)

    for (i in 1:n){
      rownames(A[[i]]) <- rownames(As[[i]]) <- colnames(x[[i]])
      colnames(A[[i]]) <- colnames(As[[i]]) <- paste0("comp", 1:ncol(A[[i]]))
    }
  }
  else {
    scores <- NULL
  }

  if(max(ns)==min(ns))
    corsY <- mclapply(x, function(x, y) cor(x, y), y=Y, mc.cores=mc.cores)
  else{
    ff <- function(x, y){
      o <- intersect(rownames(x), rownames(y))
      ans <- cor(x[o,], y[o,])
      ans
    }
    corsY <- mclapply(x, ff, y=Y, mc.cores=mc.cores)
  }

  if (pval)
    pval.cor <- mclapply(corsY, cor.test.p, n=m, mc.cores=mc.cores)
  else
    pval.cor <- NULL

  if(is.null(names(x)))
    names(x) <- paste0("df", 1:length(x))

  names(corsY) <- names(x)

  if(!is.null(scores))
    names(scores) <- names(x)
  if(pval)
    names(pval.cor) <- names(x)

  AVE_X <- lapply(corsY, function(x) apply(x^2, 2, mean))
  outer <- matrix(unlist(AVE_X), nrow = nfac)
  AVE_outer <- sapply(1:nfac, function(j, p) sum(p * outer[j,])/sum(p),
                      p=p)
  AVE_inner <- Re(eig$v[1:nfac])
  AVE <- list(AVE_X = AVE_X,
              AVE_outer_model = AVE_outer,
              AVE_inner_model = AVE_inner)

  ans <- list(Y=Y, corsY=corsY, scores=scores,
              pval.cor=pval.cor, AVE=AVE)
  class(ans) <- "mgcca"
  ans
}
