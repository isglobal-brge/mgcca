#' Imputes each assay matrix data from a MultiAssayExperiment object
#' @param method ...
#' @param remove.col remove cols
#' @param remove.row remove rows
#' @export
#' @importFrom impute impute.knn

impute <- function(multiassayexperiment, method, remove.col = FALSE,
                   impute.zero = FALSE, rowmax = 0.5, colmax = 0.8, ...){

  # Check that the input is a MultiAssayExperiment
  if (!class(multiassayexperiment) == "MultiAssayExperiment")
    stop("Input must be a 'MultiAssayExperiment' object \n")

  # Check that method is provided
  inv.type <- c("knn", "hmisc")
  inv.method <- charmatch(method, inv.type, nomatch = 0)
  if (inv.method == 0)
    stop("method should be 'knn' or 'hmisc' \n")

  # KNN method
  if (inv.method == 1)
    for (assay in 1:length(multiassayexperiment)) {
      matrix_to_impute <- as.matrix(assays(multiassayexperiment)[[assay]])
      if (impute.zero)
        # Transform 0 to NA (doesn't work for 0.000?)
        matrix_to_impute = replace(matrix_to_impute, which(matrix_to_impute == 0), NA)
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













