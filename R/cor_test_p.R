#' S3 method for cor.test on class 'p'
#' @description Compute two-sided p-value from correlation and sample size.
#' @param r Correlation coefficient.
#' @param n Sample size (>= 3).
#' @param ... Further arguments (ignored).
#' @return Numeric p-value (two-sided).
#' @export
#' @method cor.test p
cor.test.p <- function(r, n)
{

  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(pt(abs(t),(n-2), lower.tail=FALSE))
  p

}

