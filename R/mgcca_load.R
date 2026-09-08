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
.mgcca_write_manifest <- function(res, method, lambda,
                                  input_group = NULL, scale = NULL) {
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
    # The results manifest records where the RESULTS live. The reliability API
    # (mgcca_sensitivity / mgcca_block_profile / mgcca_stability) also needs the
    # SOURCE blocks -- the sensitivity kernel reads them from HDF5 and stability
    # refits -- so the INPUT group must be recorded too, or a reloaded fit cannot
    # be used for anything that touches the data. Written here rather than assumed
    # to be "MGCCA_IN", which is only the default and is user-configurable.
    if (!is.null(input_group)) run$input_group <- as.character(input_group)
    if (!is.null(scale))       run$scale       <- as.integer(isTRUE(scale))
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
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:100]
#' num <- function(d, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% ids, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = num(X1, 1:10), clinical = num(X2), other = num(X3))
#'
#' ## A file-backed fit. collect = FALSE leaves the results in the file only.
#' h5 <- tempfile(fileext = ".h5")
#' desc <- mgcca(X, filename = h5, nfac = 2, method = "penalized",
#'               lambda = rep(0.1, 3), collect = FALSE)
#' rm(desc)
#'
#' ## Rebuilt from the file alone -- provenance included.
#' fit <- mgcca_load(h5)
#' fit
#' attr(fit, "desc")$method
#' attr(fit, "desc")$lambda
#' unlink(h5)
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
            mgcca_date    = chr(man$mgcca_date),
            input_group   = if (is.null(man$input_group)) NA_character_
                            else as.character(man$input_group),
            scale         = if (is.null(man$scale)) NA
                            else as.logical(as.integer(man$scale)))
        # NB `input_group` stays NA when the fit predates it. mgcca_load() does NOT
        # warn here: it is the general loader, and provenance of the INPUT blocks
        # only matters to the reliability functions. Warning every caller about a
        # field they will never use is noise, and it changed the behaviour of a
        # documented, tested path. The check belongs where it is actionable --
        # the reliability entry points fail closed on a missing input group.
    } else {
        ## no manifest -> best effort from the file layout
        datasets <- .mgcca_discover_datasets(file, group)
        if (!length(datasets))
            stop("no manifest and no corsY tables found under '", group,
                 "' in ", file)
        top <- tryCatch(mgcca_list_group_rcpp(file, group),
                        error = function(e) character(0))
        BigDataStatMeth::hdf5_close_all()
        # No manifest at all: provenance is unknown. `input_group` stays NA on
        # purpose so the reliability entry points can fail closed on it; this
        # function stays silent, exactly as before.
        desc <- list(
            filename    = file,
            datasets    = datasets,
            final_group = group,
            scores      = "scores" %in% top,
            route       = NA_character_,
            method      = NA_character_,
            lambda      = NULL,
            input_group = NA_character_,
            scale       = NA)
    }

    obj <- mgcca_results(desc)
    attr(obj, "desc") <- desc
    obj
}
