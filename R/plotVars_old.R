plotVars.old <- function(object, axis1=1, axis2=2, main, table.names,
                    tcolors, legend.show, cor.thres, lab.cex=2, point.cex=1) {
  x <- object$corsY
  all0 <- do.call(rbind, x)
  ii <- which(duplicated(rownames(all0)))
  rownames(all0)[ii] <- paste0(rownames(all0)[ii], "rep")

  n <- sapply(x, nrow)
  if(missing(table.names))
    table.names <- paste0("table", 1:length(x))
  all <- data.frame(x=all0[,axis1], y = all0[,axis2],
                    feature=rep(table.names, n))
  rownames(all) <- rownames(all0)


  ## /

  ## Filters for 0s
  f1 <- all[ , 1] != 0 & all[ , 2] != 0
  f2 <- all[ , 1] == 0 | all[ , 2] == 0
  f3 <- abs(all[ , 1]) > cor.thres | abs(all[ , 2]) > cor.thres
  ## /

  ## Draw Circle
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- data.frame(xcircle = cos(theta), ycircle = sin(theta))
  cplot <- ggplot2::ggplot(data = circle,
                           ggplot2::aes_string(paste0("Comp", axis1),
                                               paste0("Comp", axis2))) +
    ggplot2::geom_path(ggplot2::aes_string("xcircle", "ycircle"), color="darkgray") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color="darkcyan") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color="darkcyan")
  ## /

  ## Draw important points
  if(sum(f1) != 0) {
    cplot <- cplot + ggplot2::geom_point(data = all[f1, ],
                                         ggplot2::aes_string("x", "y",
                                                             color = "feature"),
                                         shape = 19, size = point.cex)
  }
  ## /

  ## Draw points zeros?
  if(sum(f2) != 0) {
    cplot <- cplot + ggplot2::geom_point(data = all[f2, ],
                                         ggplot2::aes_string("x", "y",
                                                             color = "feature"),
                                         shape = 19, size = point.cex) +
      ggplot2::scale_color_manual(values=tcolors) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(paste("Axis", axis2)) +
      ggplot2::scale_x_continuous(paste("Axis", axis1))
  }
  ## /

  ## Add labels for features
  all$label <- rownames(all)
  cplot <- cplot + ggrepel::geom_text_repel(
    data = subset(all, abs(all$x) >= cor.thres | abs(all$y) >= cor.thres),
    ggplot2::aes_string("x", "y", label="label"),
    size = lab.cex,
    box.padding = ggplot2::unit(0.35, "lines"),
    point.padding = ggplot2::unit(0.3, "lines"),
    color="#222222",
    segment.color="#BBBBBB"
  ) + ggplot2::theme_bw() + ggplot2::theme(legend.position="none")
  ## /

  ## Extra bar charts
  uplot <- ggplot2::ggplot(all, ggplot2::aes_string(x="x", y="x", color="feature",
                                                    fill="feature")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::scale_color_manual(values=tcolors) +
    ggplot2::scale_fill_manual(values=tcolors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position="none",
      axis.title.x=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_blank(),
      axis.ticks.x=ggplot2::element_blank()
    ) + ggplot2::ylab(paste("Axis", axis1)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color="darkcyan")
  rplot <- ggplot2::ggplot(all, ggplot2::aes_string(x="y", y="y", color="feature", fill="feature")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::scale_color_manual(values=tcolors) +
    ggplot2::scale_fill_manual(values=tcolors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title=ggplot2::element_blank(),
      axis.title.y=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank()
    ) + ggplot2::ylab(paste("Axis", axis2)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color="darkcyan") +
    ggplot2::coord_flip()
  ##

  ## Add titles
  if(missing(main)) {
    uplot <- uplot + ggplot2::ggtitle(paste0("GCCA\n", paste(names(x), collapse = " - "))) +
      ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 0.8, face = "bold"))
  } else {
    uplot <- uplot + ggplot2::ggtitle(main) +
      ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 0.8, face = "bold"))
  }
  ## /

  ## Add legend
  if(legend.show == FALSE) {
    rplot <- rplot +  ggplot2::theme(legend.position="none")
  }
  ##

  ## Empty
  empty <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(1,1), colour="white")+
    ggplot2::theme(
      axis.ticks=ggplot2::element_blank(),
      panel.background=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.title.x=ggplot2::element_blank(),
      axis.title.y=ggplot2::element_blank()
    )
  ##
  return(gridExtra::grid.arrange(uplot, empty, cplot, rplot, ncol=2, nrow=2,
                                 widths=c(4, 1), heights=c(1, 4)))
}

