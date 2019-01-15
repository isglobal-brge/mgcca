#' Generalized canonical correlation with missing individuals
#'
#' @param x an object of class 'mgcca'
#'
#' @export


getSignif <- function(x, df=NA, pval.cut=0.05, ...){
  if (!inherits(x, "mgcca"))
    stop("mgcca object expected for 'x' argument")

  xx <- x$pval.cor

  if(is.null(xx))
      stop("no p-values are available. Run again 'mgcca' with 'pval=TRUE'")

  datasets <- names(xx)
  if(is.na(df))
    df <- 1:length(xx) else
      if (!all(df %in% 1:length(datasets)) & !all(df %in% datasets))
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
