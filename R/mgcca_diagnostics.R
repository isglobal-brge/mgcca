## Fit diagnostics: they REPORT on a fit, they never change it.
##
## Two quantities, both already implied by what the pipeline computes, neither of
## them previously visible to the user:
##
##   (1) the EXTERNAL SPECTRAL GAP mu_L - mu_{L+1} of the mgcca operator, which
##       run_eigen already writes to EIGEN/MKsum05 (absolute, relative, plus the
##       eigen-residual). When that gap sits at machine zero the leading
##       L-dimensional subspace is not identified: the solver converged, but the
##       pair of eigenvectors it returned is arbitrary. Convergence is not
##       identifiability, and only the gap tells them apart.
##
##   (2) the PAIRWISE OVERLAP n_jk (individuals shared by two tables) and its
##       normalised version alpha_jk = n_jk / n over the union. The dispersion of
##       alpha across pairs is the structural quantity that governs how much the
##       treatment of missing individuals can matter: with alpha flat and close
##       to 1 every sensible treatment agrees, and the choice only starts to bite
##       when the pairs disagree about who they have in common.
##
## Both read small datasets (1 x 1, and m x 1 per table) and neither touches the
## estimator, so no computed quantity of a fit can move because of them.


## ---- thresholds ------------------------------------------------------------
## Degeneracy screen for the external gap. The band between the two regimes is
## empty and it is wide: every broken configuration measured so far sits at
## machine zero -- 2.2e-16 to 8.9e-16 on the fixtures shipped with the package,
## and up to 5.8e-14 as the worst of 200 simulated replicates -- while the
## smallest gap of a configuration that is NOT broken is 2.3e-4 (the penalized
## twin of the broken p > n fixture; the healthy fixtures run 1.7e-3 to 9.9e-2).
## That is more than nine orders of magnitude of empty space, and 1e-8 sits
## inside it with about five orders of clearance above the broken side and four
## below the healthy side, so the screen does not depend on where exactly it is
## placed. It is also the cut the package already applies to the same phenomenon
## when it decides whether the top eigenvalues of a fixture are degenerate
## (tests/testthat/test-methods.R).
##
## The relative gap is the primary screen because it is scale-free. The absolute
## one is a safety net for the pathological case where mu_L itself collapses, in
## which the ratio stops being informative; the operator's spectrum is bounded by
## the number of tables, so an absolute cut is meaningful for it.
.mgcca_gap_tol_rel <- 1e-8
.mgcca_gap_tol_abs <- 1e-10

## Reporting heuristics for the overlap notes. These carry NO inferential status:
## they decide whether a sentence is printed, nothing else, and the quantity each
## one is based on is printed beside it so the reader can judge for themselves.
##
##   * Heterogeneity is screened on the COEFFICIENT OF VARIATION of alpha, not on
##     its standard deviation, and the reason is mechanical rather than cosmetic.
##     What a missing-data treatment does to the pairwise associations is
##     reweight each of them by its own overlap; a reweighting that is the same
##     for every pair is a global scalar and cancels out of the eigenproblem up
##     to scale. So a set of alphas all equal to 0.4 is as benign as a set all
##     equal to 0.95, and what bites is the spread of alpha RELATIVE to its
##     level. cv > 0.25 -- the pairs differ by more than a quarter of their mean
##     overlap -- is the cut; it is a display threshold, not a test.
##   * min n_jk < 30 -- the conventional floor below which a pairwise association
##     is barely estimable at all.
.mgcca_alpha_cv_flag <- 0.25
.mgcca_n_jk_flag     <- 30L


## ---- small HDF5 readers ----------------------------------------------------
## Read a 1 x 1 (or first) numeric value; NA_real_ if it is not there.
.mgcca_read_scalar <- function(filename, path) {
    hm <- BigDataStatMeth::hdf5_matrix(filename, path)
    on.exit(if (is.function(hm$close)) try(hm$close(), silent = TRUE), add = TRUE)
    v <- as.numeric(as.matrix(hm))
    if (length(v)) v[1] else NA_real_
}

## Presence matrix (m x J logical): the non-zero entries of the diagonal of K_j,
## i.e. exactly the individuals table j holds. This is the same source the C++
## layer uses for every per-block statistic (mgcca::present_rows), so the report
## and the estimator cannot disagree about who is present.
.mgcca_presence <- function(filename, datasets, tmp_group = "MGCCA_TMP") {
    handles <- list()
    on.exit({
        for (h in handles)
            if (!is.null(h) && is.function(h$close)) try(h$close(), silent = TRUE)
    }, add = TRUE)
    cols <- lapply(datasets, function(ds) {
        hm <- BigDataStatMeth::hdf5_matrix(filename, paste0(tmp_group, "/K/", ds))
        handles[[length(handles) + 1L]] <<- hm
        as.numeric(as.matrix(hm))
    })
    ms <- vapply(cols, length, integer(1))
    if (!length(ms) || length(unique(ms)) != 1L)
        stop("the K masks of the tables do not share a single union of individuals")
    matrix(unlist(cols) != 0, nrow = ms[1], ncol = length(datasets),
           dimnames = list(NULL, datasets))
}


## ---- (1) external spectral gap ---------------------------------------------
## The diagnostics run_eigen wrote for this fit. NULL when the file does not
## carry them (a results-only file, or one written before they existed).
.mgcca_eigen_report <- function(filename, group = "EIGEN/MKsum05") {
    out <- tryCatch(
        list(gap      = .mgcca_read_scalar(filename, paste0(group, "/eiggap")),
             gap_rel  = .mgcca_read_scalar(filename, paste0(group, "/eiggap_rel")),
             residual = .mgcca_read_scalar(filename, paste0(group, "/eig_residual"))),
        error = function(e) NULL)
    if (is.null(out) || all(vapply(out, function(v) !is.finite(v), logical(1))))
        return(out)
    out
}

## TRUE when the external gap says the leading subspace is not identified.
.mgcca_subspace_degenerate <- function(eig) {
    if (is.null(eig)) return(FALSE)
    gr <- eig$gap_rel; ga <- eig$gap
    (is.finite(gr) && gr < .mgcca_gap_tol_rel) ||
        (is.finite(ga) && ga < .mgcca_gap_tol_abs)
}

## The "chivato de rotura": a warning, and only a warning. The fit is returned
## exactly as computed -- nothing is altered, nothing is aborted.
.mgcca_warn_degenerate <- function(eig, nfac) {
    msg <- paste0(
        "the leading ", nfac, "-dimensional subspace is not numerically ",
        "identified at this size.\n",
        "  External spectral gap mu_", nfac, " - mu_", nfac + 1L, " = ",
        format(eig$gap, digits = 3), " (", format(eig$gap_rel, digits = 3),
        " relative), i.e. at machine zero: the eigenvalues that define the top\n",
        "  subspace coincide, so which directions the solver returned is ",
        "arbitrary even though it converged.\n",
        "  This is the regime in which a table carries about as many variables ",
        "as it has observed rows (p_j vs observed rows).\n",
        "  The results of this configuration are unreliable. ",
        "method = \"penalized\" (with a lambda on the scale of the tables)\n",
        "  restores the separation and is the recommended configuration here. ",
        "The fit itself is returned unchanged.")
    warning(warningCondition(msg, class = "mgcca_degenerate_subspace",
                             call = NULL))
    invisible(TRUE)
}


## ---- (2) pairwise overlap --------------------------------------------------
## n_jk, alpha_jk = n_jk / n and the dispersion of alpha across pairs. Computed
## in R from the K masks after the fit: they are m x 1 vectors, so this costs a
## J-column read and a J x J crossproduct whatever the size of the tables.
.mgcca_overlap_report <- function(filename, datasets, tmp_group = "MGCCA_TMP") {
    P <- tryCatch(.mgcca_presence(filename, datasets, tmp_group),
                  error = function(e) NULL)
    .mgcca_overlap_core(P, datasets)
}

## The arithmetic of the overlap report, factored out so the in-memory path can
## feed the same presence matrix (m x J logical, columns = datasets) it holds in
## RAM and get a byte-identical report. The file-backed reader above builds P from
## the HDF5 K masks; this is otherwise unchanged from what it always computed.
.mgcca_overlap_core <- function(P, datasets) {
    if (is.null(P) || ncol(P) < 2L) return(NULL)

    # base:: on purpose: BigDataStatMeth is attached (Depends) and masks
    # crossprod/diag for the user; the package namespace resolves to base today,
    # but an import(BigDataStatMeth) added later would silently change what these
    # two lines mean. They are plain J x J integer arithmetic and must stay so.
    n <- nrow(P)
    N <- base::crossprod(P + 0)                 # J x J counts; diag = n_present
    n_present <- stats::setNames(as.integer(base::diag(N)), datasets)

    ij     <- utils::combn(ncol(P), 2L)
    pairs  <- data.frame(table1 = datasets[ij[1, ]],
                         table2 = datasets[ij[2, ]],
                         n_jk   = as.integer(N[t(ij)]),
                         stringsAsFactors = FALSE)
    pairs$alpha_jk <- pairs$n_jk / n

    # With a single pair there is nothing to disperse; sd() would be NA and an NA
    # would then have to be special-cased by every consumer. 0 is the honest
    # value: the one pair cannot disagree with another.
    a  <- pairs$alpha_jk
    sd <- if (length(a) > 1L) stats::sd(a) else 0
    cv <- if (mean(a) > 0) sd / mean(a) else 0

    list(n             = as.integer(n),
         n_present     = n_present,
         pairs         = pairs,
         alpha_sd      = sd,
         alpha_cv      = cv,
         alpha_min     = min(a),
         alpha_max     = max(a),
         n_jk_min      = min(pairs$n_jk),
         heterogeneous = isTRUE(cv > .mgcca_alpha_cv_flag),
         sparse_pair   = isTRUE(min(pairs$n_jk) < .mgcca_n_jk_flag))
}
