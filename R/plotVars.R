#' Plot variables correlations with shared canonical components
#'
#' @param x an object of class 'mgcca'
#'
#'
#' @export
#' @importFrom made4 plotgenes



plotVars <- function(x, var=NA, axes=1:2,
                     var.col="red", # the length either 1 or length(var)
                     var.lab=FALSE, # T or F
                     bg.var.col="gray", # the length either 1 or length(df)
                     nlab=0,
                     df=NA, # either name of data.frame or numeric
                     layout=NA,
                     tit = "",
                     tit.pos = 0,
                     tit.cex = 1,...){


  if (!inherits(x, "mgcca"))
    stop("mgcca object expected for 'x' argument")
  if (length(axes) != 2)
    stop("you have to select (only) 2 axis")
  if (length(var.col) != 1 & length(var.col) != length(var))
    stop("the length of var.col could only be either 1 or length(var)")
  if (length(bg.var.col) != 1 & length(bg.var.col) != length(df))
    stop("the length of bg.var.col could only be either 1 or length(df)")

  datasets <- names(x$corsY)

  if (missing(df))
    df <- 1:length(datasets) else
      if (!all(df %in% 1:length(datasets)) & !all(df %in% datasets))
        stop("undefined data.frame selected")


  df.list <- lapply(x$corsY, function(x, axes) x[,axes], axes=axes)

  n <- length(df)
  n <- ceiling(n/2)*2

  #   ORIGINAL CODE
  # if (is.matrix(layout))
  #   layout(layout) else
  #     if (is.na(layout)) {
  #       if (length(df) == 1)
  #         layout(1) else
  #           if (length(df) == 2)
  #             layout(t(t(1:2))) else
  #               if (length(df) == 3)
  #                 layout(t(t(1:3))) else
  #                   if (length(df) > 3)
  #                     layout(matrix(1:n, n/2, byrow=T))
  # }

  if (is.matrix(layout)) {
    layout(layout)
  } else if (is.na(layout)) {
    if (length(df) == 1) {
      layout(1)
    } else if (length(df) == 2) {
      layout(t(t(1:2)))
    } else if (length(df) == 3) {
      layout(t(t(1:3)))
    } else if (length(df) > 3) {
      layout(matrix(1:n, n/2, byrow=T))
    }
  }

  vars <- var
  par(mar=c(0.1, .1, .1, .1))
  for (i in df){
    idf <- data.frame(df.list[[i]])
    ns <- rownames(idf)
    made4::plotgenes(idf, axis1=1, axis2=2, nlab=nlab, genelabels=ns,
              colpoints=bg.var.col, ...)
    ind <- ns %in% var
    if (any(ind)) {
      points(idf[ind, ], pch=20, col=var.col[na.omit(match(ns, var))])
      if (var.lab)
        text(idf[ind, 1], idf[ind, 2], ns[ind])
    }
    legend(x="bottomleft", bty="n", legend=datasets[i], x.intersp=-.5)
    title(main = tit, line = tit.pos, cex.main = tit.cex )
    vars <- cbind(vars, var %in% ns)
  }
  if (!is.na(vars[1])) {
    vars <- as.data.frame(vars)
    colnames(vars) <- c("Variables", "Dataset")
    if (!any(as.logical(vars[,2]))){
      cat("There are variables names not in your tables \n")
      print(vars)
    }
  }
}

