#' Variables significantly correlated with the shared components
#'
#' @description Collects, from the per-table p-values stored in
#'   \code{x$pval.cor}, the variables whose correlation with \emph{at least one}
#'   shared canonical component is significant at level \code{pval.cut}. The fit
#'   must therefore carry p-values (see \code{\link{mgcca}}).
#'
#' @param x an object of class 'mgcca'
#' @param df tables to inspect: a vector of positions in \code{x$pval.cor} or of
#'   table names. The default \code{NA} means all tables.
#' @param pval.cut significance level. A variable is reported when any of its
#'   component p-values is strictly below this value. Default 0.05.
#' @param ... further arguments, currently ignored.
#' @return A \code{data.frame} with one row per selected variable and two
#'   character columns: \code{variable} (the variable name) and \code{table}
#'   (the table it comes from). A variable significant in more than one table
#'   appears once per table.
#' @seealso \code{\link{topVars}}, \code{\link{plotVariables}}
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
#' sig <- getSignif(fit, pval.cut = 0.01)
#' head(sig)
#' table(sig$table)
#'
#' ## only the clinical table
#' getSignif(fit, df = "clinical", pval.cut = 0.01)
#' @export


getSignif <- function(x, df=NA, pval.cut=0.05, ...){
  if (!inherits(x, "mgcca"))
    stop("mgcca object expected for 'x' argument")

  xx <- x$pval.cor

  if(is.null(xx))
      stop("no p-values are available. Run again 'mgcca' with 'pval=TRUE'")

  datasets <- names(xx)
  if(is.na(df))
    df <- seq_along(xx) else
      if (!all(df %in% seq_along(datasets)) & !all(df %in% datasets))
        stop("selected table in 'df' is not a valid name. Try a number or \n
             a correct name")

  ns.sig <- NULL
  for(i in df){
    idf <- xx[[i]]
    ns <- rownames(idf)
    if (is.character(i))
      datasets.i <- i
    else
      datasets.i <- datasets[i]
    ns.sig.i <- cbind(ns[apply(idf < pval.cut, 1, any)], datasets.i)
    ns.sig <- rbind(ns.sig, ns.sig.i)
  }
  ans <- data.frame(ns.sig)
  colnames(ans) <- c("variable", "table")
  ans
}
