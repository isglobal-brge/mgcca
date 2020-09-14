#' Plot variables correlations with shared canonical components
#'
#' @param x an object of class 'mgcca'
#'
#'
#' @export
#' @importFrom ggthemes geom_rangeframe


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
    col.list <- mycols[1:levs]
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
