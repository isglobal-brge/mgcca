#' Generalized Canonical Correlation Analysis with missing individuals
#'
#' @description High-level, one-call entry point for the C++/HDF5 MGCCA pipeline:
#'   import any supported input into HDF5, then run the single-call orchestrator
#'   \code{\link{mgcca_rcpp}}. Results are written to the HDF5 file under
#'   \code{FINAL_RESULTS} (\code{Y}, \code{corsY}, \code{pval}, \code{AVE}, and
#'   \code{scores} if requested); with \code{collect = TRUE} they are also
#'   returned in memory as an object of class \code{"mgcca"} ready for
#'   \code{\link{plotIndividuals}} / \code{\link{plotVariables}}.
#'
#' @param x Input data: an HDF5 file path (data already under \code{group}), a
#'   \code{MultiAssayExperiment}, an \code{ExpressionSet}/
#'   \code{SummarizedExperiment}, or a named \code{list} of matrices
#'   (individuals x variables, rownames = individual IDs).
#' @param filename Target HDF5 file (ignored when \code{x} is an HDF5 path).
#' @param group Input HDF5 group. Default \code{"MGCCA_IN"}.
#' @param datasets Optional dataset names (subset/order).
#' @param nfac Number of shared components. Default 2.
#' @param scale If TRUE (default), column-center+scale each table before
#'   analysis. Scaling is performed \emph{inside} the HDF5 file (out-of-core,
#'   base-R \code{scale()} semantics) after a raw, untransformed import, so it
#'   applies uniformly to every input type -- including inputs that are already
#'   an HDF5 file -- and never loads a full table into RAM.
#' @param method Inversion method: \code{"solve"} (SPD Cholesky),
#'   \code{"penalized"} (requires \code{lambda}), or \code{"geninv"}/\code{"ginv"}
#'   (Moore-Penrose pseudoinverse).
#' @param lambda Numeric vector (length = number of tables) for
#'   \code{method = "penalized"}.
#' @param scores If TRUE, also compute per-table weights and scores.
#' @param route Per-table algebra route: \code{"auto"} (default) picks per table
#'   by \code{min(n, p)} (covariance \code{X'X} when \code{p < n}, Gram
#'   \code{XX'} when \code{p >= n}); \code{"cov"} or \code{"dual"} force one route
#'   for all tables. The Gram/dual route is the one that scales to \code{p >> n}
#'   (e.g. full methylation) without ever forming a \code{p x p} matrix.
#' @param threads Optional thread count.
#' @param overwrite If TRUE (default), overwrite the file/datasets on import.
#' @param collect If TRUE (default), read the results back from HDF5 and return an
#'   in-memory object of class \code{"mgcca"} (via \code{\link{mgcca_results}})
#'   ready for \code{\link{plotIndividuals}} / \code{\link{plotVariables}} /
#'   \code{\link{getSignif}}; the HDF5 descriptor is kept on its \code{"desc"}
#'   attribute. With \code{collect = FALSE} the lightweight descriptor list is
#'   returned instead (results still live in the HDF5 file).
#'
#' @return With \code{collect = TRUE} (default), an object of class
#'   \code{"mgcca"} (see \code{\link{mgcca_results}}); its HDF5 descriptor is on
#'   \code{attr(., "desc")}. With \code{collect = FALSE}, the descriptor list
#'   (\code{filename}, \code{datasets}, \code{nfac}, \code{m}, \code{eig_values},
#'   \code{route}, ...). Results always live in the file under
#'   \code{FINAL_RESULTS}.
#' @seealso \code{\link{mgcca_results}}, \code{\link{plotIndividuals}},
#'   \code{\link{mgcca_import_hdf5}}, \code{\link{mgcca_rcpp}}
#' @examples
#' \dontrun{
#' data(cardiovascular)
#' X <- list(methylation = as.matrix(X1),
#'           clinical    = as.matrix(X2),
#'           other       = as.matrix(X3))
#' fit <- mgcca(X, filename = tempfile(fileext = ".h5"),
#'              method = "solve", scores = TRUE, collect = TRUE)
#' plotIndividuals(fit)
#' }
#' @export
mgcca <- function(x, filename, group = "MGCCA_IN", datasets = NULL, nfac = 2,
                  scale = TRUE, method = "penalized", lambda = NULL,
                  scores = FALSE, route = c("auto", "cov", "dual"),
                  threads = NULL, overwrite = TRUE, collect = TRUE) {

    route <- match.arg(route)
    inv <- switch(method, solve = 1L, penalized = 2L, geninv = 3L, ginv = 3L,
                  stop("method must be 'solve', 'penalized', 'geninv' or 'ginv'"))
    if (inv == 2L && (is.null(lambda)))
        stop("method 'penalized' requires 'lambda'")

    # Import RAW, untransformed: any transformation belongs in HDF5, not in R
    # memory. Scaling (if requested) is done out-of-core by mgcca_rcpp before
    # getK, so it applies uniformly whatever the input type -- including inputs
    # that are already an HDF5 file.
    desc <- mgcca_import_hdf5(x, filename = filename, group = group,
                              datasets = datasets, overwriteFile = overwrite,
                              overwriteDataset = overwrite)

    res <- mgcca_rcpp(desc$filename, desc$group, desc$datasets,
                      nfac = as.integer(nfac), inv_method = inv, lambda = lambda,
                      scores = scores, scale = isTRUE(scale), route = route,
                      threads = threads)

    # Persist the provenance manifest as native HDF5 attributes so the results
    # file is self-describing and reloadable via mgcca_load() without this
    # session. Done whether or not we collect, so collect = FALSE files carry it.
    .mgcca_write_manifest(res, method = method, lambda = lambda)

    if (!isTRUE(collect))
        return(res)

    obj <- mgcca_results(res)
    attr(obj, "desc") <- res
    obj
}
