#' Convert a character table in a list of tables to numeric
#'
#' @description Some assays (e.g. methylation from TCGA) are stored as character
#'   matrices; \code{mgcca} needs numeric tables. This converts the table at
#'   position \code{index} in a list of tables to numeric, preserving dimnames.
#' @param listMAE a list of tables (e.g. from \code{\link{getTables}}).
#' @param index integer position of the table to convert.
#' @return The input list of tables (still of class \code{"ListMAE"}) with the
#'   table at position \code{index} replaced by a numeric matrix of the same
#'   dimensions, row names and column names. All other tables are returned
#'   unchanged.
#' @seealso \code{\link{getTables}}
#' @examples
#' if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
#'   data(cardiovascular)
#'   sel <- rownames(X2)[1:20]
#'   ## two assays, features in rows, individuals in columns
#'   a1 <- t(as.matrix(X1[rownames(X1) %in% sel, 1:8, drop = FALSE]))
#'   a2 <- t(as.matrix(X2[rownames(X2) %in% sel, , drop = FALSE]))
#'   ## a character assay, as methylation downloaded from TCGA often is
#'   storage.mode(a1) <- "character"
#'   mae <- MultiAssayExperiment::MultiAssayExperiment(
#'       MultiAssayExperiment::ExperimentList(
#'           list(methylation = a1, clinical = a2)))
#'
#'   tabs <- getTables(mae)
#'   class(tabs[[1]][1, 1])          # "character"
#'   tabs <- matrix.chr2num(tabs, 1)
#'   class(tabs[[1]][1, 1])          # "numeric"
#' }
#' @export
matrix.chr2num <- function(listMAE, index) {

  # Check that the input is of ListMAE class
  if (!inherits(listMAE, "ListMAE"))
    stop("Input must be a 'ListMAE' object \n")

  # Check index parameter has been given
  if (missing(index))
    stop("Index parameter indicating the position of the matrix in the list must be given \n")

  # First, keep rownames and colnames from the original matrix
  chr.colnames <- colnames(listMAE[[index]])
  chr.rownames <- rownames(listMAE[[index]])
  chr.numeric <- as.numeric(listMAE[[index]])
  # Create new matrix with the same values and dimensions of the original matrix. By default, values will be numeric.
  num.matrix <- matrix(data=chr.numeric, ncol=dim(listMAE[[index]])[2],
              nrow=dim(listMAE[[index]])[1])
  # Reassign colnames and rownames from the original matrix
  colnames(num.matrix) <- chr.colnames
  rownames(num.matrix) <- chr.rownames
  # Update the matrix in the list
  listMAE[[index]] <- num.matrix
  listMAE
}



