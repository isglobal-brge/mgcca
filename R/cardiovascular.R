#' Cardiovascular multi-table example data
#'
#' @description A small, self-contained example for \code{\link{mgcca}}: three
#'   tables measured on partly overlapping sets of individuals, illustrating the
#'   missing-individuals setting mgcca is designed for. Loading the data set
#'   creates three matrices in the environment.
#'
#' @details Each table has individuals in rows (identified by their row names,
#'   sample barcodes) and variables in columns. The three tables do not share
#'   exactly the same individuals (union of 646; 606 common to all three):
#'   \describe{
#'     \item{\code{X1}}{618 x 69 methylation beta-values (CpG columns).}
#'     \item{\code{X2}}{633 x 6 clinical measurements (e.g. HDL, LDL, BMI, glucose).}
#'     \item{\code{X3}}{645 x 4 cell-type composition estimates.}
#'   }
#'
#' @format Three numeric matrices, \code{X1}, \code{X2} and \code{X3}.
#' @name cardiovascular
#' @docType data
#' @keywords datasets
#' @usage data(cardiovascular)
#' @examples
#' data(cardiovascular)
#' X <- list(methylation = as.matrix(X1),
#'           clinical    = as.matrix(X2),
#'           other       = as.matrix(X3))
#' sapply(X, dim)
NULL
