plotInds <- function(x, group, ax1=1, ax2=2, col.list, print.labels=FALSE,
                    cex.label=0.8, pos.leg = "bottomright", ...){
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
  plot(comp1, comp2, type="n", ...)
  points(comp1, comp2, pch=16, col=cols, xlab=paste0("Global component", comp1),
         ylab=paste0("Global component", comp2))
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
