plot.mgcca <- function(x, table=1, nvars=5) {
  if(!inherits(x, "mgcca"))
   stop(" 'x' must be an object of class mgcca")

  As <- x$A
  kk <- table
  ww1 <- order(abs(As[[kk]][,1]))[1:nvars]
  ww2 <- order(abs(As[[kk]][,2]))[1:nvars]
  ww <- union(ww1,ww2)
  ww <- c(1:nrow(As[[kk]]))%in%ww
  plot(As[[kk]][,1], As[[kk]][,2], col="white", pch=19,
         xlab=expression(As[1]), ylab=expression(As[2]))
  abline(h=0,v=0,lty=2,col="grey")
  text(As[[kk]][,1], As[[kk]][,2], rownames(As[[kk]]), col=ifelse(ww,"red","black"),cex=0.8)
  invisible()
}




