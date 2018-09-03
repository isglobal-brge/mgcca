cor.test.p <- function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(pt(abs(t),(n-2), lower.tail=FALSE))
  p
}

