#' Splits a MultiAssayExperiment in a list of tables

getTables <- function(multiassayexperiment){

  tables_list = list()
  
  for (assay in 1:length(multiassayexperiment)) {
    
    matrix_to_add <- as.matrix(assays(multiassayexperiment)[[assay]]) 
    print(names(assays(multiassayexperiment)[assay]))
    print(class(assays(multiassayexperiment)[[assay]]))
    # This is for TGCA data, to obtain patient ID's, because patient full names contains
    # the type of data aswell (RNA-seq, miRNA, etc), the first 12 characters is the ID.
    colnames(matrix_to_add) <- substr(colnames(matrix_to_add),1,12)
    tables_list[[names(assays(multiassayexperiment)[assay])]] <- t(matrix_to_add) 
    
  }
  
  ans <- tables_list
  
}

### --> TO DO
# Samples IDs should be the same for all tables, and in TGCA it is not