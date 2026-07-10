#' Splits a MultiAssayExperiment in a list of tables
#'
#' @param multiassayexperiment MultiAssayExperiment
#' @return List of tables
#' @importFrom SummarizedExperiment assays
#' @export
getTables <- function(multiassayexperiment){

  # Check that the input is a MultiAssayExperiment
  if (!class(multiassayexperiment) == "MultiAssayExperiment")
    stop("Input must be a 'MultiAssayExperiment' object \n")

  tables.list = list()

  for (assay in 1:length(multiassayexperiment)) {
    matrix.add <- as.matrix(assays(multiassayexperiment)[[assay]])
    tables.list[[names(assays(multiassayexperiment)[assay])]] <- t(matrix.add)
  }

  class(tables.list) <- "ListMAE"
  tables.list

}
