#' Print method for 'ListMAE' objects
#'
#' @description Compact print for a MultiAssayExperiment-like list, as returned
#'   by \code{\link{getTables}}: reports the class, the number of assays and
#'   their names.
#'
#' @param x Object of class \code{"ListMAE"}, as returned by
#'   \code{\link{getTables}}.
#' @param ... Further arguments (ignored).
#'
#' @return Invisibly the assay names of \code{x}; called for the side effect of
#'   printing the summary.
#'
#' @seealso \code{\link{getTables}}, which produces the object printed here.
#'
#' @examples
#' if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
#'   data(cardiovascular)
#'   sel <- rownames(X2)[1:20]
#'   a1 <- t(as.matrix(X1[rownames(X1) %in% sel, 1:8, drop = FALSE]))
#'   a2 <- t(as.matrix(X2[rownames(X2) %in% sel, , drop = FALSE]))
#'   mae <- MultiAssayExperiment::MultiAssayExperiment(
#'       MultiAssayExperiment::ExperimentList(
#'           list(methylation = a1, clinical = a2)))
#'
#'   tabs <- getTables(mae)
#'   class(tabs)          # "ListMAE"
#'   tabs                 # dispatches to print.ListMAE()
#' }
#'
#' @export
#' @method print ListMAE
print.ListMAE <- function(x, ...) {
    print(paste("Object of class ", class(x), sep = ""))
    print(paste(length(x), " assays in the MultiAssayExperiment list:", sep = ""))
    print(names(x))
}
