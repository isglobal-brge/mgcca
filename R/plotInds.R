#' Base-R scatter plot of the individuals on the shared components
#'
#' @description Draws the individuals in the plane of two shared canonical
#'   components, optionally coloured by a grouping factor and labelled with the
#'   individual identifiers. This is the original base-graphics plot;
#'   \code{\link{plotIndividuals}} is the modern replacement and returns a
#'   \code{ggplot} object instead of drawing on the current device.
#'
#' @param x an object of class 'mgcca'
#' @param group optional factor giving the group of each individual, in the row
#'   order of \code{x$Y}. When omitted every individual is put in a single
#'   group and no legend is drawn.
#' @param ax1,ax2 the components on the horizontal and vertical axes. Defaults 1
#'   and 2.
#' @param col.list vector of colours, one per level of \code{group}. When
#'   omitted a built-in palette is used; it must be at least as long as the
#'   number of levels.
#' @param print.labels logical; when \code{TRUE} the individual names are
#'   written instead of plotting points (using \pkg{wordcloud} to spread the
#'   labels out when that package is available). Default \code{FALSE}.
#' @param cex.label character expansion for those labels. Default 0.8.
#' @param pos.leg position of the legend, as in \code{\link[graphics]{legend}}.
#'   Default \code{"bottomright"}.
#' @param main plot title. Default \code{NULL}.
#' @param ... further graphical arguments passed to \code{\link[graphics]{plot}}.
#' @return No return value; called for the plot it draws on the current
#'   graphics device.
#' @seealso \code{\link{plotIndividuals}}, \code{\link{plotVars}}
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
#' plotInds(fit, main = "Individuals")
#'
#' ## coloured by an arbitrary two-level grouping
#' grp <- factor(rep(c("a", "b"), length.out = nrow(fit$Y)))
#' plotInds(fit, group = grp, col.list = c("steelblue", "tomato"))
#' @export
#' @importFrom ggthemes geom_rangeframe
#' @importFrom graphics abline grid legend points text
#' @importFrom grDevices colors


plotInds <- function(x, group, ax1=1, ax2=2, col.list, print.labels=FALSE,
                    cex.label=0.8, pos.leg = "bottomright", main = NULL, ...){
  if (!inherits(x, "mgcca"))
    stop("x must be an object of class 'rgcca'")
  if (missing(group))
    group <- as.factor(rep(1, nrow(x$Y)))
  else {
    if (!is.factor(group))
      stop("'group must be a factor variable")
  }

  ll <- levels(group)
  levs <- length(ll)
  comp1 <- x$Y[, ax1]
  comp2 <- x$Y[, ax2]
  if (missing(col.list)){
    mycols <-  c("red", "blue", "darkgreen", "orange", "violet", sample(colors()))
    col.list <- mycols[seq_len(levs)]
  }
  if (length(col.list) < levs)
    stop("'col.list' length should be equal to the levels of grouping variable")
  cols <- as.character(factor(group, labels=col.list))

  # p <- ggplot(as.data.frame(na.omit(x$Y)), aes( x = comp1,
  #                                      y = comp2,
  #                                      color = group)
  #             ) +
  #   geom_point(  size = 1 ) +
  #   labs(
  #     title = main,
  #     x = paste("Global component", ax1),
  #        y = paste("Global component", ax2)) +
  #   theme(
  #     text = element_text( size = 18),
  #     legend.background = element_blank(),
  #     legend.key = element_blank(),
  #     panel.background = element_blank(),
  #     panel.border = element_blank(),
  #     strip.background = element_blank(),
  #     plot.background = element_blank(),
  #     # axis.line = element_blank(),
  #     panel.grid = element_blank(),
  #     # plot.title = element_text(size = 12, hjust = 0.5),
  #     legend.title=element_blank(),
  #     legend.position = "bottom"
  #
  #     ) +
  #   geom_rangeframe() +
  #   geom_vline(xintercept = 0, linetype="dashed", color = "grey") +
  #   geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
  #   scale_color_manual(values=col.list)
  # print(p)


  # Equivalen codi superior amb ggplot2
  #
  plot(comp1, comp2, type="n", xlab=paste("Global component", ax1),
       ylab=paste("Global component", ax2), ...)
  if (!print.labels)
    points(comp1, comp2, pch=16, col=cols)
  abline(h=0)
  abline(v=0)
  grid(lty=3, col="gray80")
  if (length(ll) > 1)
    legend(pos.leg, legend=ll, pch=16, col=col.list)


  if (print.labels) {
    if (requireNamespace("wordcloud", quietly = TRUE)) {
      wordcloud::textplot(x = comp1, y = comp2,
                          words = names(comp1),
                          cex = cex.label,
                          new = FALSE, col=cols)
    } else {
      text(comp1, jitter(comp2), names(comp1),
           cex = cex.label, adj = 0, col=cols)
    }
  }
}
