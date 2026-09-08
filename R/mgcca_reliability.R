#' Present-only, feature-streamed participant block Gram (reliability layer, K1)
#'
#' @description Computes the participant Gram \eqn{G = Z^\top Z} (\code{n_pr x
#'   n_pr}) of one HDF5-backed block, where \eqn{Z} standardizes each feature
#'   (row) to mean 0 and unit variance using the \strong{present participants
#'   only} and the \eqn{(n_{pr}-1)} sample-SD denominator. Features are streamed
#'   in chunks so the full \eqn{p \times n_{pr}} standardized matrix is never
#'   materialized (only the active chunk and the \code{n_pr x n_pr} Gram are
#'   retained). This is the first kernel (K1) of the mgcca reliability-layer
#'   C++/HDF5 port; it reproduces the sealed R oracle
#'   (\code{34_methyl_gram_dense.R}) within scale-aware tolerance.
#'
#' @details Block-absent participants must NEVER enter centring/scaling/cross
#'   products as observed zero rows (the frozen zero-padding regression); pass
#'   the true present-only participant set in \code{present_ids}. Degenerate
#'   (constant) features trigger a hard stop reporting their index, matching the
#'   oracle convention (no silent variance floor, drop, or rescale).
#'
#' @param file,group,dataset HDF5 location of the block, stored R-view
#'   \code{p x N} (features x all-participants) with participant IDs as colnames.
#' @param present_ids Character vector of present participant IDs in the frozen
#'   order; selects the columns that enter the Gram.
#' @param chunk Integer number of features per streamed chunk (default 4096).
#' @param var_eps Degenerate-feature threshold on the per-feature SD (default
#'   1e-8, the oracle value); \code{sd <= var_eps} stops with the feature index.
#' @param out Optional \code{list(file, group, dataset, compression)} to persist
#'   \eqn{G} to HDF5 entirely in C++ (the file is created if missing, participant
#'   dimnames carried, compression default 0).
#'
#' @return A list with \code{G} (n_pr x n_pr, participant dimnames), \code{ids},
#'   \code{p_eff}, \code{n_pr}, \code{N_all}, \code{chunk}, \code{n_chunks},
#'   \code{last_chunk}, \code{eigen_threads}, an \code{invariants} list (e_trace,
#'   e_center, sym_err, neg_mass, min_eig_ratio, num_rank), and a per-chunk
#'   \code{ledger}.
#' @keywords internal
mgcca_block_gram <- function(file, group, dataset, present_ids,
                             chunk = 4096L, var_eps = 1e-8, out = NULL) {
    stopifnot(is.character(present_ids), length(present_ids) >= 2L,
              length(chunk) == 1L, chunk >= 1L)
    of <- og <- od <- ""; oc <- 0L
    if (!is.null(out)) {
        stopifnot(all(c("file", "dataset") %in% names(out)))
        of <- out$file; og <- if (is.null(out$group)) "" else out$group
        od <- out$dataset; oc <- if (is.null(out$compression)) 0L else as.integer(out$compression)
    }
    reliability_gram_hdf5(file, group, dataset, as.character(present_ids),
                          as.integer(chunk), as.numeric(var_eps),
                          out_file = of, out_group = og, out_dataset = od,
                          out_compression = oc)
}
