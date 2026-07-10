#' Collect mgcca results from HDF5 into an in-memory 'mgcca' object
#'
#' @description The C++/HDF5 core (\code{\link{mgcca_rcpp}} / \code{\link{mgcca}})
#'   writes its results to the HDF5 file under \code{FINAL_RESULTS} rather than
#'   returning them in memory. This reader materialises those results into a
#'   compact list of class \code{"mgcca"} -- the same shape the legacy in-memory
#'   \code{mgcca()} returned -- so the post-analysis helpers
#'   \code{\link{plotInds}}, \code{\link{plotVars}}, \code{\link{getSignif}} and
#'   \code{\link{topVars}} work unchanged on the port's output.
#'
#'   The shared components \code{Y} and the per-table \code{corsY} / \code{pval}
#'   are small (individuals or variables \eqn{\times} \code{nfac}), so collecting
#'   them into RAM does not defeat the out-of-core design; the large tables and
#'   intermediates stay on disk.
#'
#' @param x Either the descriptor list returned by \code{mgcca()} /
#'   \code{mgcca_rcpp()} (with elements \code{filename}, \code{datasets},
#'   \code{final_group}), or a length-1 character path to the HDF5 file.
#' @param datasets Optional table names to read. Defaults to \code{x$datasets}
#'   when \code{x} is a descriptor, otherwise required.
#' @param final_group HDF5 group holding the results. Default
#'   \code{x$final_group} or \code{"FINAL_RESULTS"}.
#' @param scores,pval Whether to read the per-table scores / p-values (they exist
#'   only if \code{mgcca()} was called with \code{scores = TRUE} / the
#'   defaults). Default: infer from the descriptor, otherwise \code{TRUE}.
#'
#' @return A list of class \code{"mgcca"} with elements \code{Y} (matrix,
#'   individuals \eqn{\times} components), \code{corsY} and \code{pval.cor}
#'   (named lists of variable \eqn{\times} component matrices), \code{scores}
#'   (named list or \code{NULL}) and \code{AVE} (\code{AVE_X},
#'   \code{AVE_outer_model}, \code{AVE_inner_model}).
#' @seealso \code{\link{mgcca}}, \code{\link{plotInds}}, \code{\link{plotVars}}
#' @export
mgcca_results <- function(x, datasets = NULL, final_group = NULL,
                          scores = NULL, pval = TRUE) {

    if (is.list(x)) {
        filename    <- x$filename
        if (is.null(datasets))    datasets    <- x$datasets
        if (is.null(final_group)) final_group <- x$final_group
        if (is.null(scores))      scores      <- isTRUE(x$scores)
    } else if (is.character(x) && length(x) == 1L) {
        filename <- x
        if (is.null(scores)) scores <- TRUE
    } else {
        stop("'x' must be an mgcca descriptor list or an HDF5 file path")
    }
    if (is.null(final_group)) final_group <- "FINAL_RESULTS"
    if (is.null(datasets) || !length(datasets))
        stop("'datasets' is required when 'x' does not carry them")
    if (!file.exists(filename)) stop("HDF5 file not found: ", filename)

    open_handles <- list()
    read_mat <- function(path) {
        hm <- BigDataStatMeth::hdf5_matrix(filename, path)
        open_handles[[length(open_handles) + 1L]] <<- hm
        m  <- as.matrix(hm)
        dn <- dimnames(hm)
        if (!is.null(dn))
            dimnames(m) <- list(if (length(dn[[1]])) as.character(dn[[1]]) else NULL,
                                if (length(dn[[2]])) as.character(dn[[2]]) else NULL)
        m
    }
    on.exit({
        for (h in open_handles)
            if (!is.null(h) && is.function(h$close)) try(h$close(), silent = TRUE)
    }, add = TRUE)

    per_table <- function(sub) {
        out <- lapply(datasets, function(ds) read_mat(paste0(final_group, "/", sub, "/", ds)))
        names(out) <- datasets
        out
    }
    # optional per-table blocks that may be absent (scores=FALSE / scale=FALSE)
    try_per_table <- function(sub) {
        out <- lapply(datasets, function(ds)
            tryCatch(read_mat(paste0(final_group, "/", sub, "/", ds)),
                     error = function(e) NULL))
        if (all(vapply(out, is.null, logical(1)))) return(NULL)
        names(out) <- datasets
        out
    }

    Y       <- read_mat(paste0(final_group, "/Y"))
    corsY   <- per_table("corsY")
    pvalc   <- if (isTRUE(pval))   per_table("pval")   else NULL
    scr     <- if (isTRUE(scores)) per_table("scores") else NULL
    weights <- try_per_table("weights")     # p x nfac (present iff scores = TRUE)
    scaling <- try_per_table("scaling")     # p x 2 [center, scale] (iff scale = TRUE)

    AVE <- list(
        AVE_X           = as.matrix(read_mat(paste0(final_group, "/AVE/AVE_X"))),
        AVE_outer_model = as.numeric(read_mat(paste0(final_group, "/AVE/AVE_outer"))),
        AVE_inner_model = as.numeric(read_mat(paste0(final_group, "/AVE/AVE_inner"))))

    ans <- list(Y = Y, corsY = corsY, scores = scr, pval.cor = pvalc,
                weights = weights, scaling = scaling, AVE = AVE)
    class(ans) <- "mgcca"
    ans
}
