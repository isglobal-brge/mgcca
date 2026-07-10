## Provenance manifest persistence and reload for mgcca results.
##
## The numeric results already live in the HDF5 file under FINAL_RESULTS. What
## was NOT persisted is the run's provenance (method, lambda, route, ...). We
## store it as native HDF5 attributes so a results file is self-describing and
## reloadable without the original R session:
##   * run-level facts        -> attributes of the group  FINAL_RESULTS
##   * per-table facts (hybrid)-> attributes of each dataset FINAL_RESULTS/corsY/<ds>
##     (lambda, route_dual), because they are genuinely one-per-table.
## Logical values are stored as integer 0/1 (HDF5 has no native boolean) and
## rebuilt with as.logical() here.

## Write the provenance manifest for a finished run. `res` is the descriptor
## returned by mgcca_rcpp(); `method`/`lambda` are the mgcca() call arguments.
.mgcca_write_manifest <- function(res, method, lambda) {
    fg <- if (is.null(res$final_group)) "FINAL_RESULTS" else res$final_group
    datasets <- as.character(res$datasets)

    run <- list(
        mgcca_version = as.character(utils::packageVersion("mgcca")),
        mgcca_date    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        method        = as.character(method),
        route         = as.character(res$route),
        nfac          = as.integer(res$nfac),
        scores        = as.integer(isTRUE(res$scores)),  # logical -> int
        m             = as.integer(res$m),
        datasets      = datasets,
        eig_values    = as.numeric(res$eig_values))
    mgcca_write_attrs_rcpp(res$filename, fg, "", run)
    BigDataStatMeth::hdf5_close_all()

    rd <- res$route_dual                      # named logical vector or NULL
    for (i in seq_along(datasets)) {
        ds  <- datasets[i]
        val <- if (!is.null(rd) && ds %in% names(rd)) rd[[ds]] else rd[i]
        per <- list(route_dual = as.integer(isTRUE(val)))
        if (identical(as.character(method), "penalized") && !is.null(lambda))
            per$lambda <- as.numeric(lambda[[min(i, length(lambda))]])
        mgcca_write_attrs_rcpp(res$filename, paste0(fg, "/corsY"), ds, per)
        BigDataStatMeth::hdf5_close_all()
    }
    invisible(TRUE)
}

## Auto-discover table names from FINAL_RESULTS/corsY subgroups (fallback for
## files written before the manifest existed). Drops the .<ds>_dimnames helpers.
.mgcca_discover_datasets <- function(file, group) {
    kids <- mgcca_list_group_rcpp(file, paste0(group, "/corsY"))
    BigDataStatMeth::hdf5_close_all()
    kids <- kids[!grepl("^\\.", kids)]
    as.character(kids)
}

#' Reload an mgcca result straight from its HDF5 file
#'
#' @description Rebuilds a full \code{"mgcca"} object from a results file written
#'   by \code{\link{mgcca}}, using \emph{only} the file: it reads the provenance
#'   manifest stored as native HDF5 attributes (method, lambda, route,
#'   route_dual, nfac, datasets, eigenvalues, ...) and combines it with
#'   \code{\link{mgcca_results}} (which reads the numeric blocks). No original R
#'   session, descriptor list or \code{.rds} is needed -- the file is
#'   self-describing.
#'
#'   If the file predates the manifest (no attributes), \code{mgcca_load} falls
#'   back to auto-discovering the table names from the \code{corsY} subgroups and
#'   returns what it can, with the provenance fields set to \code{NA}.
#'
#' @param file Path to the HDF5 results file.
#' @param group HDF5 group holding the results. Default \code{"FINAL_RESULTS"}.
#'
#' @return An object of class \code{"mgcca"} (see \code{\link{mgcca_results}}),
#'   with the reconstructed provenance descriptor on \code{attr(., "desc")}.
#' @seealso \code{\link{mgcca}}, \code{\link{mgcca_results}}
#' @examples
#' \dontrun{
#' fit <- mgcca(X, filename = "res.h5", method = "penalized", lambda = c(0.75, 1))
#' rm(fit)
#' fit2 <- mgcca_load("res.h5")   # same object, straight from disk
#' }
#' @export
mgcca_load <- function(file, group = "FINAL_RESULTS") {
    if (!is.character(file) || length(file) != 1L)
        stop("'file' must be a single HDF5 file path")
    if (!file.exists(file)) stop("HDF5 file not found: ", file)

    man <- tryCatch(mgcca_read_attrs_rcpp(file, group, ""),
                    error = function(e) list())
    BigDataStatMeth::hdf5_close_all()

    if (length(man) && !is.null(man$datasets)) {
        datasets <- as.character(man$datasets)

        route_dual <- stats::setNames(rep(NA, length(datasets)), datasets)
        lambda     <- stats::setNames(rep(NA_real_, length(datasets)), datasets)
        for (ds in datasets) {
            pa <- tryCatch(
                mgcca_read_attrs_rcpp(file, paste0(group, "/corsY"), ds),
                error = function(e) list())
            BigDataStatMeth::hdf5_close_all()
            if (!is.null(pa$route_dual)) route_dual[ds] <- as.logical(pa$route_dual)
            if (!is.null(pa$lambda))     lambda[ds]     <- as.numeric(pa$lambda)
        }
        if (all(is.na(lambda))) lambda <- NULL

        num  <- function(v) if (is.null(v)) NULL else as.numeric(v)
        int  <- function(v) if (is.null(v)) NULL else as.integer(v)
        chr  <- function(v) if (is.null(v)) NA_character_ else as.character(v)

        desc <- list(
            filename      = file,
            datasets      = datasets,
            nfac          = int(man$nfac),
            m             = int(man$m),
            eig_values    = num(man$eig_values),
            scores        = if (is.null(man$scores)) TRUE else as.logical(man$scores),
            route         = chr(man$route),
            route_dual    = route_dual,
            final_group   = group,
            method        = chr(man$method),
            lambda        = lambda,
            mgcca_version = chr(man$mgcca_version),
            mgcca_date    = chr(man$mgcca_date))
    } else {
        ## no manifest -> best effort from the file layout
        datasets <- .mgcca_discover_datasets(file, group)
        if (!length(datasets))
            stop("no manifest and no corsY tables found under '", group,
                 "' in ", file)
        top <- tryCatch(mgcca_list_group_rcpp(file, group),
                        error = function(e) character(0))
        BigDataStatMeth::hdf5_close_all()
        desc <- list(
            filename    = file,
            datasets    = datasets,
            final_group = group,
            scores      = "scores" %in% top,
            route       = NA_character_,
            method      = NA_character_,
            lambda      = NULL)
    }

    obj <- mgcca_results(desc)
    attr(obj, "desc") <- desc
    obj
}
