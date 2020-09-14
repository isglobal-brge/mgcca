#' Generalized canonical correlation with missing individuals
#'
#' @param x an object of class 'mgcca'
#' @param axis axis to take in to account to get the top values
#' @param topN number of values to return with topVars if not pval.cut is introduced
#' @param pval.cut significance level to cut data to return
#' @export

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
    nvars <- sapply(a, nrow)
    topN <- min(topN, min(nvars))
  }
  else{
    a <- x$pval.cor
    a2 <- x$corsY
  }

  ntables <- length(a)
  tops <- list()
  for (i in 1: ntables){
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
