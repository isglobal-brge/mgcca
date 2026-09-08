#' Cardiovascular multi-table example data
#'
#' @description A small, self-contained example for \code{\link{mgcca}}: three
#'   tables measured on partly overlapping sets of individuals, illustrating the
#'   missing-individuals setting mgcca is designed for. Loading the data set
#'   creates three data frames, \code{X1}, \code{X2} and \code{X3}, in the
#'   calling environment.
#'
#' @details Each table has individuals in rows (identified by their row names,
#'   Illumina sample barcodes) and variables in columns; every column is
#'   numeric and there are no missing values within a table. The three tables
#'   do not share exactly the same individuals -- the union is 646 individuals
#'   and only 606 are common to all three, which is the incomplete-blocks
#'   situation \code{\link{mgcca}} handles:
#'   \describe{
#'     \item{\code{X1}}{618 x 69 data frame of methylation beta values (columns
#'       \code{cg1}--\code{cg69}).}
#'     \item{\code{X2}}{633 x 6 data frame of clinical measurements: \code{TC},
#'       \code{LDL}, \code{HDL}, \code{BMI}, \code{cintura} (waist
#'       circumference) and \code{Gluc} (glucose).}
#'     \item{\code{X3}}{645 x 4 data frame of blood cell-type composition
#'       estimates: \code{CD4T}, \code{NK}, \code{Mono} and \code{Gran}.}
#'   }
#'
#'   \code{\link{mgcca}} expects numeric matrices, so coerce the tables with
#'   \code{as.matrix()} before fitting (see the examples).
#'
#' @format Three numeric data frames with individuals in rows and variables in
#'   columns:
#'   \describe{
#'     \item{\code{X1}}{\code{data.frame} with 618 rows and 69 columns
#'       (methylation beta values).}
#'     \item{\code{X2}}{\code{data.frame} with 633 rows and 6 columns (clinical
#'       measurements).}
#'     \item{\code{X3}}{\code{data.frame} with 645 rows and 4 columns
#'       (cell-type composition).}
#'   }
#' @name cardiovascular
#' @aliases cardiovascular X1 X2 X3
#' @docType data
#' @keywords datasets
#' @usage data(cardiovascular)
#' @examples
#' data(cardiovascular)
#' X <- list(methylation = as.matrix(X1),
#'           clinical    = as.matrix(X2),
#'           other       = as.matrix(X3))
#' vapply(X, dim, integer(2))
NULL
