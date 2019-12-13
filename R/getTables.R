#' Splits a MultiAssayExperiment in a list of tables

getTables <- function(multiassayexperiment){

  tables_list = list()
  
  for (assay in 1:length(multiassayexperiment)) {
    
    matrix_to_add <- as.matrix(assays(multiassayexperiment)[[assay]]) 
    print(names(assays(multiassayexperiment)[assay]))
    print(class(assays(multiassayexperiment)[[assay]]))
    tables_list[[names(assays(multiassayexperiment)[assay])]] <- t(matrix_to_add) 
    
  }
  
  tables_list
  
}