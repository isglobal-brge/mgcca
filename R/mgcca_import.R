#' Import data of any supported type into an HDF5 file for mgcca
#'
#' @description Versatile entry point that normalizes the input tables to the
#'   mgcca HDF5 convention and writes them under a single group. Accepts data
#'   already stored in HDF5, a \code{MultiAssayExperiment}, an
#'   \code{ExpressionSet}/\code{SummarizedExperiment}, or a plain \code{list} of
#'   R matrices. Each table is stored as **individuals x variables** (R-view)
#'   with \code{rownames} = individual IDs (required to align tables and detect
#'   missing individuals) and \code{colnames} = variable names.
#'
#' @details Dimnames are written in the BigDataStatMeth layout
#'   (\code{<group>/.<dataset>_dimnames/1}=colnames, \code{/2}=rownames) by
#'   \code{\link[BigDataStatMeth]{bdCreate_hdf5_matrix}}, which preserves the
#'   \code{dimnames} attribute of the object automatically. For assay-based
#'   inputs (features x samples), assays are transposed to individuals x
#'   variables. A plain \code{list} is assumed to already be individuals x
#'   variables (as produced by \code{\link{getTables}}).
#'
#' @param x Input data: a length-1 character path to an existing HDF5 file, a
#'   \code{MultiAssayExperiment}, an \code{ExpressionSet}/
#'   \code{SummarizedExperiment}, or a named \code{list} of matrices.
#' @param filename Target HDF5 file. Ignored (and may be missing) when \code{x}
#'   is itself an HDF5 path.
#' @param group HDF5 group to write the tables under. Default \code{"MGCCA_IN"}.
#' @param datasets Optional character vector of dataset names. When \code{x} is
#'   an existing HDF5 file, restricts/orders the tables to use; otherwise names
#'   the written datasets (defaults to the names of \code{x}).
#' @param overwriteDataset Logical, overwrite datasets if present. Default FALSE.
#' @param overwriteFile Logical, overwrite/create the file. Default FALSE.
#'
#' @return A descriptor \code{list} with elements \code{filename}, \code{group}
#'   and \code{datasets} (the dataset names written/used), for downstream mgcca.
#' @seealso \code{\link{getTables}}
#' @export
mgcca_import_hdf5 <- function(x, filename, group = "MGCCA_IN", datasets = NULL,
                              overwriteDataset = FALSE, overwriteFile = FALSE) {

    # ---- (a) already an HDF5 file: use as-is -------------------------------
    if (is.character(x) && length(x) == 1L) {
        if (!file.exists(x))
            stop("'x' looks like an HDF5 path but the file does not exist: ", x)
        ds <- datasets
        if (is.null(ds))
            ds <- BigDataStatMeth::bdgetDatasetsList_hdf5(filename = x, group = group)
        return(list(filename = x, group = group, datasets = as.character(ds)))
    }

    # ---- other input types: coerce and write ONE table at a time -----------
    # A lazy source lets us materialise each table (individuals x variables) only
    # when it is about to be written, and free it before the next -- so the peak
    # extra RAM is one table (+ its transpose), not the whole list at once.
    src <- .mgcca_table_source(x)
    nms <- src$names
    if (is.null(nms) || any(nms == "")) nms <- paste0("table", seq_along(nms))
    if (!is.null(datasets)) {
        if (!all(datasets %in% nms))
            stop("'datasets' must be a subset of the table names in 'x'")
        sel <- match(datasets, nms); nms <- datasets
    } else {
        sel <- seq_along(nms)
    }
    if (length(sel) < 2L)
        stop("mgcca needs at least two tables")

    if (missing(filename) || is.null(filename))
        stop("'filename' is required unless 'x' is an existing HDF5 file")
    if (!dir.exists(dirname(filename)))
        stop("destination directory does not exist: ", dirname(filename))
    if (isTRUE(overwriteFile) && file.exists(filename))
        file.remove(filename)

    # hdf5_create_matrix() writes the dimnames from the object's dimnames
    # attribute; orientation is handled internally, so we pass individuals x
    # variables (R-view) as-is.
    for (k in seq_along(sel)) {
        m <- src$get(sel[k])
        if (is.null(rownames(m)))
            stop("table '", nms[k], "' has no rownames (individual IDs required)")
        hm <- BigDataStatMeth::hdf5_create_matrix(
            filename = filename,
            dataset  = paste0(group, "/", nms[k]),
            data     = m,
            dtype    = "double",
            overwrite = overwriteDataset)
        # close the read handle the R6 object leaves open, then free the table
        if (!is.null(hm) && is.function(hm$close)) hm$close()
        rm(m); invisible(gc(FALSE))
    }

    list(filename = filename, group = group, datasets = nms)
}


# Lazy table source: returns the table names plus a getter that materialises one
# table (individuals x variables, storage double) on demand. Lets the importer
# process tables one at a time instead of holding them all in RAM at once.
.mgcca_table_source <- function(x) {
    dbl <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }

    if (is.list(x) && !methods::is(x, "MultiAssayExperiment")) {
        # plain list of matrices, assumed individuals x variables already
        return(list(names = names(x), get = function(j) dbl(x[[j]])))
    }
    if (methods::is(x, "MultiAssayExperiment")) {
        if (!requireNamespace("MultiAssayExperiment", quietly = TRUE))
            stop("package 'MultiAssayExperiment' is required for this input")
        ex <- MultiAssayExperiment::experiments(x)
        return(list(names = names(ex),
                    get = function(j) dbl(t(SummarizedExperiment::assay(ex[[j]])))))
    }
    if (methods::is(x, "SummarizedExperiment"))
        return(list(names = "assay",
                    get = function(j) dbl(t(SummarizedExperiment::assay(x)))))
    if (methods::is(x, "ExpressionSet"))
        return(list(names = "exprs",
                    get = function(j) dbl(t(Biobase::exprs(x)))))

    stop("unsupported input type for mgcca_import_hdf5: ", class(x)[1])
}

# Eager variant: coerce a supported input into a named list of all tables.
# Used where every table is needed at once (e.g. mgcca_permtest()).
.mgcca_as_tables <- function(x) {
    src <- .mgcca_table_source(x)
    nms <- src$names
    if (is.null(nms) || any(nms == "")) nms <- paste0("table", seq_along(nms))
    out <- lapply(seq_along(nms), function(j) src$get(j))
    names(out) <- nms
    out
}
