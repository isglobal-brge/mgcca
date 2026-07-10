#' Print method for 'listMAE' objects
#' @description Compact print for a MultiAssayExperiment-like list.
#' @param listMAE Object of class 'listMAE'.
#' @param ... Further arguments (ignored).
#' @export
#' @method print listMAE
print.listMAE <- function(listMAE) {
    print(paste("Object of class ", class(listMAE), sep = ""))
    print(paste(length(listMAE), " assays in the MultiAssayExperiment list:", sep = ""))
    print(names(listMAE))
}
