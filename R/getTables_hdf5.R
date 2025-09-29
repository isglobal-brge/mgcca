#' Split a MultiAssayExperiment into HDF5 tables
#'
#' @description Write each assay from a \code{MultiAssayExperiment} to an
#'   HDF5 file (group \code{"MGCCA_IN"}) and return a list of tables.
#'
#' @param multiassayexperiment A \code{MultiAssayExperiment} (or compatible).
#' @param filename Character string. Path to the target HDF5 file.
#' @param overwriteFile Logical. If \code{TRUE}, overwrite/create the file.
#' @param overwriteDataset Logical. If \code{TRUE}, overwrite existing datasets.
#' @param debugnrows \code{NULL} or single integer. If \code{NULL}, use all
#'   rows; otherwise keep the first \code{debugnrows} rows of each assay.
#' @param debugncols \code{NULL} or single integer. If \code{NULL}, use all
#'   columns; otherwise keep the first \code{debugncols} columns of each assay.
#'
#' @return A named list with one entry per assay (tables written to HDF5).
#'
#' @export
getTables_hdf5 <- function(multiassayexperiment, filename,
           overwriteFile = FALSE, overwriteDataset = FALSE,
           debugncols = NULL, debugnrows = NULL)
{

    # Check if route exists
    if( !dir.exists(basename(dirname(filename))) ) {
        stop("Destination path does not exists")
    }
    # Check that the input is a MultiAssayExperiment
    if (!class(multiassayexperiment) == "MultiAssayExperiment")
        stop("Input must be a 'MultiAssayExperiment' object \n")

    tables.list <- character()
    for (assay in 1:length(multiassayexperiment)) {

        if(!is.null(debugncols) & is.numeric(debugncols)) {  nc_subset = debugncols  }
        else { nc_subset = ncol(assays(multiassayexperiment)[[assay]])  }

        if(!is.null(debugnrows) & is.numeric(debugnrows)) {  nr_subset = debugnrows }
        else { nr_subset = nrow(assays(multiassayexperiment)[[assay]]) }

        tables.list <- c(tables.list, names(assays(multiassayexperiment)[assay]))

        if(assay == 1) {
            bdCreate_hdf5_matrix( object = as.matrix(assays(multiassayexperiment)[[assay]][ 1:nr_subset, 1:nc_subset]),
                  filename = filename, group = "MGCCA_IN",
                  dataset = names(assays(multiassayexperiment)[assay]),
                  overwriteFile = overwriteFile, overwriteDataset = overwriteDataset)
            # tables.list[[assay]] <- names(assays(multiassayexperiment)[assay])
          ##..## rownames(assays(multiassayexperiment)[[1]])
          ##..## colnames(assays(multiassayexperiment)[[1]])
        } else {
            bdCreate_hdf5_matrix( object = as.matrix(assays(multiassayexperiment)[[assay]][ 1:nr_subset, 1:nc_subset]),
                                  filename = filename, group = "MGCCA_IN",
                                  dataset = names(assays(multiassayexperiment)[assay]),
                                  overwriteFile = FALSE, overwriteDataset = overwriteDataset)
            tables.list[[assay]] <- names(assays(multiassayexperiment)[assay])}
        }

    # Aquí potser estaria bé retornar els noms de les taules creades...
    # es a dir els noms dels datasets que s'han creat ??
    tables.list

}


