#' Imputes each assay matrix data from a MultiAssayExperiment object
#' @param method ...
#' @importFrom impute impute.knn

### --> TO DO 

impute <- function(multiassayexperiment, method){
  
  ### Check that the object is a multiassay experiment...
  
  # Check that method is provided
  inv.type <- c("knn", "hmisc")
  inv.method <- charmatch(method, inv.type, nomatch = 0)
  if (inv.method == 0)
    stop("method should be 'knn' or 'hmisc' \n")
  
  ### Remove columns with more than 80% of NA / 0 data...??
  #m_list[[1]][, -which(colMeans(is.na(m_list[[1]])) > 0.5)]
  #little[,-which(colMeans(is.na(little))>=0.8)]
  
  if (inv.method == 1)
    for (assay in 1:length(multiassayexperiment)) {
      matrix_to_impute <- as.matrix(assays(multiassayexperiment)[[assay]]) 
      imputed_matrix <- impute.knn(matrix_to_impute, k = 10, rowmax = 0.5, colmax = 0.8)
      multiassayexperiment[[assay]] <- imputed_matrix$data
    }
  
  ### Do the hmisc method...
    
  ans <- multiassayexperiment
}



















