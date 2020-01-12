#' Splits a MultiAssayExperiment in a list of tables

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

print.listMAE <- function(listMAE) {
  print(paste("Object of class ", class(listMAE), sep = ""))
  print(paste(length(listMAE), " assays in the MultiAssayExperiment list:", sep = ""))
  print(names(listMAE))
  }
  
  
  