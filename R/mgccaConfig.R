#' mgcca: Generalized Canonical Correlation Analysis with missing individuals
#'
#' Generalized Canonical Correlation Analysis for tables with missing
#' individuals (whole missing rows; van de Velden and Bijmolt, 2006) for
#' multi-block omics integration, running its numerical core in C++ on HDF5
#' files via the BigDataStatMeth API. Missing cells within a measured
#' individual are outside the estimator's scope and must be handled before
#' import. See \code{\link{mgcca}} to get started.
#'
#' See \code{vignette("mgcca_example", package = "mgcca")} for a worked
#' analysis and \code{vignette("mgcca_reliability", package = "mgcca")} for
#' the sensitivity and stability tools.
#'
#' Exported names follow one convention: \code{mgcca_*} for the core workflow
#' (the fit and the analyses built on it), and camelCase for the helpers
#' (data preparation, extraction and plotting).
#'
#' @name mgcca-package
#' @useDynLib mgcca, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom BigDataStatMeth hdf5_matrix hdf5_create_matrix hdf5_close_all
#'   bdgetDatasetsList_hdf5
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
