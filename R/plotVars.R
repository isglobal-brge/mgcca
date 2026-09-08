#' Base-R correlation plots of the variables, one panel per table
#'
#' @description Lays out one panel per selected table and, in each, plots the
#'   correlation of that table's variables with two shared canonical components
#'   (via \code{made4::plotgenes}). Variables named in \code{var} are
#'   highlighted, which makes it easy to follow a handful of variables of
#'   interest across tables. This is the original base-graphics plot;
#'   \code{\link{plotVariables}} is the modern per-table replacement and returns
#'   a \code{ggplot} object.
#'
#'   Requires the \pkg{made4} package (Suggests).
#'
#' @param x an object of class 'mgcca'
#' @param var character vector of variable names to highlight in every panel.
#'   Default \code{NA}, i.e. highlight nothing.
#' @param axes the two components to plot. Default \code{1:2}.
#' @param var.col colour(s) for the highlighted variables: either a single
#'   colour or one per element of \code{var}. Default \code{"red"}.
#' @param var.lab logical; write the names of the highlighted variables next to
#'   their points. Default \code{FALSE}.
#' @param bg.var.col colour(s) for the remaining ("background") variables:
#'   either a single colour or one per selected table. Default \code{"gray"}.
#' @param nlab number of the most extreme variables that
#'   \code{made4::plotgenes} labels in each panel. Default 0.
#' @param df tables to plot: a vector of positions in \code{x$corsY} or of table
#'   names. The default \code{NA} means all tables.
#' @param layout an explicit matrix passed to \code{\link[graphics]{layout}}. By
#'   default (\code{NA}) a layout is chosen from the number of tables.
#' @param tit panel title. Default \code{""}.
#' @param tit.pos title line position, as \code{line} in
#'   \code{\link[graphics]{title}}. Default 0.
#' @param tit.cex title character expansion. Default 1.
#' @param ... further arguments passed to \code{made4::plotgenes}.
#' @return No return value; called for the plots it draws on the current
#'   graphics device. When some names in \code{var} match no variable in any
#'   table, a note and a summary table are printed to the console.
#' @seealso \code{\link{plotVariables}}, \code{\link{plotInds}}
#' @examples
#' if (requireNamespace("made4", quietly = TRUE)) {
#'   data(cardiovascular)
#'   u   <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))
#'   sel <- u[1:150]
#'   mk  <- function(d, cols = NULL) {
#'       m <- as.matrix(d[rownames(d) %in% sel, , drop = FALSE])
#'       if (!is.null(cols)) m <- m[, cols, drop = FALSE]
#'       storage.mode(m) <- "double"
#'       m
#'   }
#'   X <- list(methylation = mk(X1, 1:20), clinical = mk(X2), other = mk(X3))
#'   fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3))
#'
#'   ## one table, highlighting two clinical variables
#'   plotVars(fit, df = 2, var = c("HDL", "Gluc"),
#'            var.col = c("red", "blue"), var.lab = TRUE)
#' }
#' @export
#' @importFrom graphics layout legend par points text title
#' @importFrom stats na.omit
# made4 is Suggests (plot-only); called as made4::plotgenes() below, guarded at call time.



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
    df <- seq_along(datasets) else
      if (!all(df %in% seq_along(datasets)) & !all(df %in% datasets))
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
      layout(t(t(seq_len(2))))
    } else if (length(df) == 3) {
      layout(t(t(seq_len(3))))
    } else if (length(df) > 3) {
      layout(matrix(seq_len(n), n/2, byrow=TRUE))
    }
  }

  vars <- var
  if (!requireNamespace("made4", quietly = TRUE))
    stop("plotVars() needs the 'made4' package (Suggests). Install it with ",
         "BiocManager::install('made4').")
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

