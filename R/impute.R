#' Imputes each assay matrix data from a MultiAssayExperiment object
#'
#' @description Fills in the missing values of every assay of a
#'   \code{MultiAssayExperiment} by k-nearest-neighbour imputation
#'   (\code{\link[impute]{impute.knn}}), replacing each assay in place. This is
#'   imputation \emph{within} a table (scattered missing entries); it is not the
#'   whole-row missingness across tables that \code{\link{mgcca}} handles
#'   natively.
#'
#' @param multiassayexperiment a \code{MultiAssayExperiment} object whose assays
#'   are numeric matrices.
#' @param method imputation method, given as a character string; partially
#'   matched against \code{"knn"} and \code{"hmisc"}. Only \code{"knn"} is
#'   implemented -- \code{"hmisc"} is accepted but currently leaves the object
#'   untouched. An unmatched or ambiguous value is an error.
#' @param remove.col logical; when \code{TRUE}, columns with a proportion of
#'   missing values of \code{colmax} or more are dropped from an assay before it
#'   is imputed. Default \code{FALSE}.
#' @param impute.zero logical; when \code{TRUE}, exact zeros are treated as
#'   missing values (recoded to \code{NA}) before imputation. Default
#'   \code{FALSE}.
#' @param rowmax,colmax numeric tolerances handed on to
#'   \code{\link[impute]{impute.knn}}. Note that they are passed positionally,
#'   so \code{rowmax} reaches \code{impute.knn} as its \code{k} (number of
#'   neighbours) and \code{colmax} as its \code{rowmax}. Defaults 0.5 and 0.8.
#' @param ... further arguments passed on to \code{\link[impute]{impute.knn}}.
#' @return The input \code{MultiAssayExperiment} with each assay replaced by its
#'   imputed version. When \code{method} matches \code{"hmisc"} the object is
#'   returned unchanged, that branch not being implemented.
#' @examples
#' if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
#'   data(cardiovascular)
#'   sel <- rownames(X2)[1:30]
#'   a1 <- t(as.matrix(X1[rownames(X1) %in% sel, 1:8, drop = FALSE]))
#'   a1[1, 1:3] <- NA                     # a few missing entries to fill in
#'   a2 <- t(as.matrix(X2[rownames(X2) %in% sel, , drop = FALSE]))
#'   mae <- MultiAssayExperiment::MultiAssayExperiment(
#'       MultiAssayExperiment::ExperimentList(
#'           list(methylation = a1, clinical = a2)))
#'
#'   anyNA(mae[[1]])
#'   mae.imp <- impute(mae, method = "knn")
#'   anyNA(mae.imp[[1]])
#' }
#' @export
#' @importFrom impute impute.knn

impute <- function(multiassayexperiment, method, remove.col = FALSE,
                   impute.zero = FALSE, rowmax = 0.5, colmax = 0.8, ...){

  # Check that the input is a MultiAssayExperiment
  if (!inherits(multiassayexperiment, "MultiAssayExperiment"))
    stop("Input must be a 'MultiAssayExperiment' object \n")

  # Check that method is provided
  inv.type <- c("knn", "hmisc")
  inv.method <- charmatch(method, inv.type, nomatch = 0)
  if (inv.method == 0)
    stop("method should be 'knn' or 'hmisc' \n")

  # KNN method
  if (inv.method == 1)
    for (assay in seq_along(multiassayexperiment)) {
      matrix_to_impute <- as.matrix(assays(multiassayexperiment)[[assay]])
      if (impute.zero)
        # Transform 0 to NA (doesn't work for 0.000?)
        matrix_to_impute <- replace(matrix_to_impute, which(matrix_to_impute == 0), NA)
      if (remove.col)
        # Remove columns (samples) with more than colmax of NA
        matrix_to_impute <- matrix_to_impute[, which(colMeans(is.na(matrix_to_impute)) < colmax)]
      imputed_matrix <- impute.knn(matrix_to_impute, rowmax, colmax, ...)
      multiassayexperiment[[assay]] <- imputed_matrix$data
    }

  # Hmisc method
  ### Do the hmisc method...

  multiassayexperiment

}

### no me funciona el @ImportFrom
### Remove rows with a lot of NA's/0?













