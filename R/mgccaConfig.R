#' mgcca: Generalized Canonical Correlation Analysis with missing individuals
#'
#' Generalized Canonical Correlation Analysis for tables with missing
#' individuals (whole missing rows; van de Velden and Bijmolt, 2006) for
#' multi-block omics integration, running its numerical core in C++ on HDF5
#' files via the BigDataStatMeth API. Missing cells within a measured
#' individual are outside the estimator's scope and must be handled before
#' import. See \code{\link{mgcca}} to get started.
#'
#' @name mgcca-package
#' @useDynLib mgcca, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom BigDataStatMeth hdf5_matrix hdf5_create_matrix hdf5_close_all
#'   bdgetDatasetsList_hdf5
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
