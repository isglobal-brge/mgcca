#' Get hdf5 columns and rownames
#'
#' @param filename character array indicating the name of the file to create
#' @param group path to dataset, character array indicating the complete route to the dataset
#' @param daasetname character array with dataset name
#' @return List of tables
#' @export
getDimNames_hdf5 <- function(filename, group, datasetname, sort = FALSE) {

  if(!file.exists(filename)){
    stop("file does not exists")
  }

  full_datasetname <- paste0( group, "/.", datasetname, "_dimnames")

  if( sort == TRUE ){
    return(list( rownames = sort(rhdf5::h5read( filename, paste0(full_datasetname, "/2"))),
                 colnames = sort(rhdf5::h5read( filename, paste0(full_datasetname, "/1"))) )
           )
  } else {
    return(list( rownames = rhdf5::h5read( filename, paste0(full_datasetname, "/2")),
                 colnames = rhdf5::h5read( filename, paste0(full_datasetname, "/1")) )
    )
  }
}

#' Get hdf5 row names
#'
#' @param filename character array indicating the name of the file to create
#' @param group path to dataset, character array indicating the complete route to the dataset
#' @param daasetname character array with dataset name
#' @return List of tables
#' @export
getRowNames_hdf5 <- function(filename, group, datasetname, sort = FALSE) {

  if(!file.exists(filename)){
    stop("file does not exists")
  }

  full_datasetname <- paste0( group, "/.", datasetname, "_dimnames")

  if( sort == TRUE ){
    return( sort(rhdf5::h5read( filename, paste0(full_datasetname, "/2"))) )
  } else {
    return( rhdf5::h5read( filename, paste0(full_datasetname, "/2")) )
  }

}



#' Get hdf5 column names
#'
#' @param filename character array indicating the name of the file to create
#' @param group path to dataset, character array indicating the complete route to the dataset
#' @param daasetname character array with dataset name
#' @return List of tables
#' @export
getColNames_hdf5 <- function(filename, group, datasetname, sort = FALSE) {

  if(!file.exists(filename)){
    stop("file does not exists")
  }

  full_datasetname <- paste0( group, "/.", datasetname, "_dimnames")

  if( sort == TRUE ) {
    return( sort(rhdf5::h5read( filename, paste0(full_datasetname, "/1"))) )
  } else {
    return( rhdf5::h5read( filename, paste0(full_datasetname, "/1")) )
  }

}

