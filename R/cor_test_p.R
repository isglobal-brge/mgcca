#' Two-sided p-value of a correlation coefficient
#'
#' @description Compute the two-sided p-value of a Pearson correlation
#'   coefficient from the coefficient and the sample size, without needing the
#'   underlying data. Used internally by \code{\link{mgcca_results}} to turn the
#'   between-table correlations into p-values, and exported because it is useful
#'   on its own.
#'
#' @details Despite its name (kept for backward compatibility) this is a plain
#'   function, not an S3 method: \code{stats::cor.test} is not a generic.
#'
#' @param r Correlation coefficient.
#' @param n Sample size (>= 3).
#'
#' @return A numeric vector of the same length as \code{r} with the two-sided
#'   p-values from a t distribution on \code{n - 2} degrees of freedom.
#'
#' @examples
#' cor.test.p(0.35, 100)
#' cor.test.p(c(0.1, 0.5, 0.9), 50)
#'
#' @importFrom stats pt
#' @export cor.test.p
cor.test.p <- function(r, n)
{

  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(pt(abs(t),(n-2), lower.tail=FALSE))
  p

}
