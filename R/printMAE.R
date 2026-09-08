#' Print method for 'listMAE' objects
#'
#' @description Compact print for a MultiAssayExperiment-like list, as returned
#'   by \code{\link{getTables}}: reports the class, the number of assays and
#'   their names.
#'
#' @param x Object of class \code{"listMAE"}.
#' @param ... Further arguments (ignored).
#'
#' @return Invisibly the assay names of \code{x}; called for the side effect of
#'   printing the summary.
#'
#' @examples
#' obj <- structure(list(a = matrix(1:4, 2), b = matrix(1:6, 2)),
#'                  class = "listMAE")
#' print(obj)
#'
#' @export
#' @method print listMAE
print.listMAE <- function(x, ...) {
    print(paste("Object of class ", class(x), sep = ""))
    print(paste(length(x), " assays in the MultiAssayExperiment list:", sep = ""))
    print(names(x))
}
