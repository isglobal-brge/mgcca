#' Convert a character table in a list of tables to numeric
#'
#' @description Some assays (e.g. methylation from TCGA) are stored as character
#'   matrices; \code{mgcca} needs numeric tables. This converts the table at
#'   position \code{index} in a list of tables to numeric, preserving dimnames.
#' @param listMAE a list of tables (e.g. from \code{\link{getTables}}).
#' @param index integer position of the table to convert.
#' @return The same list with table \code{index} coerced to a numeric matrix.
#' @export
matrix.chr2num <- function(listMAE, index) {

  # Check that the input is of ListMAE class
  if (!class(listMAE) == "ListMAE")
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



