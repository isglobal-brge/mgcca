#' mgcca: Generalized Canonical Correlation Analysis with missing individuals
#'
#' Generalized Canonical Correlation Analysis (van de Velden and Takane, 2012)
#' for multi-block omics integration, running its numerical core in C++ on HDF5
#' files via the BigDataStatMeth API. See \code{\link{mgcca}} to get started.
#'
#' @name mgcca-package
#' @useDynLib mgcca
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
