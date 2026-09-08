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
#' @param tmp_group HDF5 group holding the intermediates, where the per-table
#'   presence masks \code{K} live. Only used for the overlap report; if the group
#'   is absent the report is \code{NULL} and nothing else changes. Default
#'   \code{"MGCCA_TMP"}.
#' @param outputs Which components to collect. \code{NULL} (default) collects
#'   everything, which is exactly what this function has always done. Otherwise
#'   a character vector naming a subset of \code{"Y"}, \code{"corsY"},
#'   \code{"scores"}, \code{"pval"}, \code{"weights"}, \code{"scaling"},
#'   \code{"AVE"}, \code{"eigen"}, \code{"overlap"} (\code{"pval.cor"} and
#'   \code{"eig"} are accepted as synonyms of the last two). The datasets behind
#'   the components you do not ask for are never read: on a large fit collected
#'   over a slow filesystem that is where the time goes. The choice is purely
#'   subtractive -- every component you do ask for is bit-for-bit the one a full
#'   collection would have produced, because selecting outputs changes what is
#'   read and nothing about what was computed. Components not asked for come back
#'   \code{NULL}, and the object keeps its shape, so \code{$corsY} on a partial
#'   object is \code{NULL} rather than an error. What was requested is recorded
#'   on \code{attr(., "outputs")}.
#'
#' @return A list of class \code{"mgcca"} with elements \code{Y} (matrix,
#'   individuals \eqn{\times} components), \code{corsY} and \code{pval.cor}
#'   (named lists of variable \eqn{\times} component matrices), \code{scores}
#'   (named list or \code{NULL}) and \code{AVE} (\code{AVE_X},
#'   \code{AVE_outer_model}, \code{AVE_inner_model}).
#'
#'   Two read-only diagnostics are attached as well, and neither takes part in
#'   any computation: \code{eigen} (the external spectral gap
#'   \eqn{\mu_L - \mu_{L+1}} of the mgcca operator, absolute and relative, plus
#'   the eigen-residual) and \code{overlap} (the number of individuals shared by
#'   each pair of tables, \eqn{n_{jk}}, its normalised version
#'   \eqn{\alpha_{jk} = n_{jk}/n} over the union, and the dispersion of
#'   \eqn{\alpha} across pairs). \code{\link{summary}} prints both.
#'
#'   Each \code{scores} matrix spans the union of individuals, but an individual
#'   that is not part of a given table has \code{NA} in that table's scores (the
#'   HDF5 layer stores it as \code{NaN}). Only \code{Y} is defined for everyone.
#' @seealso \code{\link{mgcca}}, \code{\link{plotInds}}, \code{\link{plotVars}}
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
#' ## collect = FALSE returns the descriptor; the results stay in the file.
#' h5 <- tempfile(fileext = ".h5")
#' desc <- mgcca(X, filename = h5, nfac = 2, method = "penalized",
#'               lambda = rep(0.1, 3), collect = FALSE)
#'
#' res <- mgcca_results(desc)
#' dim(res$Y)
#' names(res$corsY)
#' res$overlap$pairs
#'
#' ## Subtractive collection: only the shared components are read back, and
#' ## they are bit-for-bit the ones a full collection would have returned.
#' resY <- mgcca_results(desc, outputs = "Y")
#' identical(resY$Y, res$Y)
#' names(Filter(Negate(is.null), resY))
#' unlink(h5)
#' @export
mgcca_results <- function(x, datasets = NULL, final_group = NULL,
                          scores = NULL, pval = TRUE,
                          tmp_group = "MGCCA_TMP", outputs = NULL) {

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

    # What to read. With outputs = NULL this is every component, i.e. exactly the
    # behaviour of every previous version; a subset only ever REMOVES reads.
    want <- .mgcca_output_plan(outputs, scores = scores, pval = pval)

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

    Y       <- if (want[["Y"]])     read_mat(paste0(final_group, "/Y")) else NULL
    corsY   <- if (want[["corsY"]]) per_table("corsY") else NULL
    pvalc   <- if (want[["pval"]])  per_table("pval")  else NULL
    # Individuals absent from a table have no score in it; the C++ layer marks
    # them NaN (HDF5 has no NA), which becomes NA_real_ here.
    na_absent <- function(S) { S[is.nan(S)] <- NA_real_; S }
    scr     <- if (want[["scores"]]) lapply(per_table("scores"), na_absent) else NULL
    # p x nfac (present iff scores = TRUE)
    weights <- if (want[["weights"]]) try_per_table("weights") else NULL
    # p x 2 [center, scale] (present iff scale = TRUE)
    scaling <- if (want[["scaling"]]) try_per_table("scaling") else NULL

    AVE <- if (want[["AVE"]]) list(
        AVE_X           = as.matrix(read_mat(paste0(final_group, "/AVE/AVE_X"))),
        AVE_outer_model = as.numeric(read_mat(paste0(final_group, "/AVE/AVE_outer"))),
        AVE_inner_model = as.numeric(read_mat(paste0(final_group, "/AVE/AVE_inner")))) else NULL

    # Read-only diagnostics. They are collected AFTER every result above and are
    # not an input to any of them: the estimator has finished by the time these
    # run, and both come from datasets the pipeline had already written. A file
    # that does not carry them (results-only, or written by an older version)
    # simply gets NULL.
    eig <- if (want[["eigen"]])   .mgcca_eigen_report(filename) else NULL
    ovl <- if (want[["overlap"]])
        .mgcca_overlap_report(filename, datasets, tmp_group) else NULL

    ans <- list(Y = Y, corsY = corsY, scores = scr, pval.cor = pvalc,
                weights = weights, scaling = scaling, AVE = AVE,
                eigen = eig, overlap = ovl)
    class(ans) <- "mgcca"
    # Only a partial collection is labelled: a full one must stay byte-identical
    # to what every previous version returned, attributes included.
    if (!is.null(outputs))
        attr(ans, "outputs") <- .mgcca_output_canonical(outputs)
    ans
}
