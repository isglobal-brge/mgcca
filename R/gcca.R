#' Generalized canonical correlation
#'
#' @param x list of matrices. Each matrix should have the ids in the
#' rownames. Missing is not allowed (see details)
#' @param nfac ...
#' @param scale ...
#' @param pval should p-values of correlation between variables and shared canonical variates be computed? Default is TRUE.
#' @param scores should canonical variables be computed for each table?
#'               Default is FALSE. See details
#' @param method ...
#' @param lambda ...
#' @param mc.cores ...
#' @details
#' eigen with Rfast hd.eigen select number of components speed up process (RSpectra?)
#' rownames matrices ... si hay distintos ordena - si no, no.
#' missings ...  dos opciones#'
#' scores .. si queremos correlation with each table!
#' inversa ...
#' "solve" para n>>p (pero lento) (solve function)
#' "penalized" penalized adding lambda*I (rfunctions cgls). Requires lambda

#' "geninv" n<<p generalized inverse (rfunctions geninv)
#' "ginv" generalized inverse (MASS ginv)
#'
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

gcca <- function(x, nfac=2, scale=TRUE, pval=TRUE, scores=FALSE,
                  method="solve", lambda, mc.cores=1, ...) {

  inv.type <- c("solve", "penalized", "geninv")
  inv.method <- charmatch(method, inv.type, nomatch = 0)
  if (inv.method == 0)
    stop("method should be 'solve' or 'penalized' \n")

  n <- length(x) # number of tables

  if (inv.method == 2){
    if (missing(lambda))
      stop("penalized method requires lambda parameter \n")
    else
      if(length(lambda)!=n)
        stop("lambda must be a vector of length equal to the number of tables \n")
  }

  nas <- sapply(x, function(x) any(is.na(x)))

  if (any(nas))
   stop("Missing values are not allowed. Either use 'impute' package or
        use tables with complete cases.")

  if (scale)
    x <- lapply(x, scale)


  if (!is.list(x))
    stop("x must be a list containing the different matrices")

  if (any(unlist(lapply(x, function(x) !is.matrix(x)))))
    x <- lapply(x, as.matrix)

  ns <- sapply(x, nrow)
  if(max(ns)!=min(ns)) # check whether there are missing individuals
    stop("Method requires complete cases or try 'mgcca' function \n
         to deal with missing individuals")

  p <- sapply(x, ncol) # number of variables per table
  numvars <- min(p) # minimum number of variables
  m <- nrow(x[[1]]) # number of individuals

  sol <- getSol(x, inv.method, lambda=lambda, mc.cores=mc.cores)
  Mi <- lapply(sol, "[[", 1)
  Minv <- lapply(sol, "[[", 2)

  M <- Reduce('+', Mi)

  if(isSymmetric(M))
    eig <- eigs_sym(M, k=nfac, ...)
  else
    eig <- eigs(M, k=nfac, ...)

  Y <- Re(eig$vectors)
  colnames(Y) <- paste0("comp", 1:ncol(Y))
  rownames(Y) <- rownames(x[[1]])

  if (scores) {
    scores <- mclapply(1:n, getScores0, Y=Y, XX=X, Minv=Minv, mc.cores=mc.cores)
    for (i in 1:n){
      rownames(scores[[i]]) <- rownames(x[[i]])
      colnames(scores[[i]]) <- paste0("comp", 1:nfac)
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
  AVE_inner <- Re(eig$values)
  AVE <- list(AVE_X = AVE_X,
              AVE_outer_model = AVE_outer,
              AVE_inner_model = AVE_inner)

  ans <- list(Y=Y, corsY=corsY, scores=scores,
              pval.cor=pval.cor, AVE=AVE)
  class(ans) <- "mgcca"
  ans
}
