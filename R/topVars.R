#' Variables most correlated with a shared component
#'
#' @description Lists, table by table, the variables at one end of a shared
#'   canonical component. By default it returns the \code{topN} variables with
#'   the largest correlation on that side; if \code{pval.cut} is supplied it
#'   returns instead every variable whose correlation is significant at that
#'   level and has the requested sign, which needs a fit carrying p-values.
#'
#' @param x an object of class 'mgcca'
#' @param axis axis to take in to account to get the top values
#' @param end which end of the component to look at: \code{"pos"} (default) for
#'   the most positively correlated variables, \code{"neg"} for the most
#'   negatively correlated ones.
#' @param topN number of values to return with topVars if not pval.cut is introduced
#' @param pval.cut significance level to cut data to return
#' @return A named list with one character vector per table of \code{x},
#'   holding the selected variable names. Without \code{pval.cut} each vector
#'   has the same length (\code{topN}, capped at the number of variables in the
#'   smallest table); with \code{pval.cut} the lengths vary and a vector may be
#'   empty.
#' @seealso \code{\link{getSignif}}, \code{\link{plotLoadings}}
#' @examples
#' data(cardiovascular)
#' u   <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))
#' sel <- u[1:150]
#' mk  <- function(d, cols = NULL) {
#'     m <- as.matrix(d[rownames(d) %in% sel, , drop = FALSE])
#'     if (!is.null(cols)) m <- m[, cols, drop = FALSE]
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = mk(X1, 1:20), clinical = mk(X2), other = mk(X3))
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3))
#'
#' ## three most positively correlated variables per table, first component
#' topVars(fit, axis = 1, topN = 3)
#'
#' ## and the negative end
#' topVars(fit, axis = 1, end = "neg", topN = 3)
#'
#' ## everything significant with a positive correlation
#' topVars(fit, axis = 1, pval.cut = 0.05)
#' @export
#' @importFrom utils head

topVars <- function(x, axis = 1, end = "pos", topN = 5, pval.cut){
  if (!inherits(x, "mgcca"))
    stop("x must be an object of class 'mgcca'")
  mm <- match(end, c("pos", "neg"), nomatch = NA)
  if (mm==2)
    side <- FALSE
  else
    side <- TRUE
  if (is.na(mm))
    stop("'end' must be 'pos' or 'neg'")

  if (missing(pval.cut)){
    a <- x$corsY
    nvars <- vapply(a, nrow, integer(1))
    topN <- min(topN, min(nvars))
  }
  else{
    a <- x$pval.cor
    a2 <- x$corsY
  }

  ntables <- length(a)
  tops <- list()
  for (i in seq_len(ntables)){
    if (missing(pval.cut)) {
      a.i <- a[[i]]
      top.i <- head(rownames(a.i)[order(a.i[,axis], decreasing = side)],
                  n=topN)
    }
    else {
      a.i <- a[[i]]
      a2.i <- a2[[i]]
      if (side)
        mask <- a.i[,axis] <= pval.cut & a2.i[,axis]>0
      else
        mask <- a.i[,axis] <= pval.cut & a2.i[,axis]<0
      top.i <- rownames(a.i)[mask]
    }
    tops[[i]] <- top.i
  }

  names(tops) <- names(a)
  tops
}
