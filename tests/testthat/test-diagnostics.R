# --- Bioconductor build time -------------------------------------------------
# This file is one of the four heaviest in the suite. Bioconductor's builders
# cap `R CMD check` at 10 minutes and this suite is the long pole, so the file
# is skipped THERE only (`IS_BIOC_BUILD_MACHINE`); it runs in full everywhere
# else, including on CRAN-style checks and in development.
testthat::skip_on_bioc()

# Fit diagnostics that report on the fit rather than change it.
#
# (1) External spectral gap.  run_eigen already computes mu_L - mu_{L+1} (and its
#     relative version) and stores it in EIGEN/MKsum05.  When that gap is at
#     machine zero the leading L-dimensional subspace is NOT identified: the
#     solver converged, but which pair of eigenvectors it returned is arbitrary.
#     That is the regime p_j >~ (rows observed in block j).  The fit must SAY so
#     and must not change because of it.
#
# (2) Pairwise overlap.  n_jk = individuals shared by blocks j and k, and
#     alpha_jk = n_jk / n over the union.  This is the structural quantity that
#     governs how much the missing-data treatment can matter; it is reported,
#     never acted upon.
#
# Both are NEW OUTPUT.  No computed quantity of the fit may move.

# ---- fixtures --------------------------------------------------------------
# Paired by construction: same union (24), same missingness (4 absent per
# block), only p_j changes.  p = 10 puts p_j above half the observed rows and
# the two rank-p_j projectors then share a subspace of dimension > nfac, so the
# top eigenvalues coincide and the external gap collapses.  p = 6 does not.
dg_pair <- function(p, seed = 20260902) {
    set.seed(seed)
    ids <- sprintf("ind%02d", 1:24)
    mk <- function(drop) {
        Z <- matrix(rnorm(24 * p), 24, p,
                    dimnames = list(ids, paste0("v", seq_len(p))))
        Z[-drop, , drop = FALSE]
    }
    list(B1 = mk(1:4), B2 = mk(21:24))
}

dg_fit <- function(X, method = "solve", lambda = NULL, ...) {
    f <- mgcca(X, filename = tempfile(fileext = ".h5"), nfac = 2,
               method = method, lambda = lambda, scores = TRUE,
               collect = TRUE, ...)
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
    f
}

# Three blocks over a union of 30 with KNOWN pairwise overlaps:
#   B1 all 30; B2 drops {1,4,9,16,25,30}; B3 drops {2,5,11,22}
#   n_12 = 24, n_13 = 26, n_23 = 30 - 10 = 20
dg_overlap_blocks <- function() {
    set.seed(20260902)
    ids <- sprintf("ind%02d", 1:30)
    B1 <- matrix(rnorm(30 * 5), 30, 5, dimnames = list(ids, paste0("v", 1:5)))
    B2 <- matrix(rnorm(30 * 4), 30, 4, dimnames = list(ids, paste0("v", 1:4)))
    B3 <- matrix(rnorm(30 * 3), 30, 3, dimnames = list(ids, paste0("v", 1:3)))
    list(B1 = B1,
         B2 = B2[-c(1, 4, 9, 16, 25, 30), , drop = FALSE],
         B3 = B3[-c(2, 5, 11, 22), , drop = FALSE])
}


# ---- (1) the gap is collected and reported ---------------------------------

test_that("the collected fit carries the external spectral gap", {
    f <- dg_fit(dg_pair(6))
    expect_false(is.null(f$eigen))
    expect_true(all(c("gap", "gap_rel", "residual") %in% names(f$eigen)))
    expect_true(is.finite(f$eigen$gap_rel))
    expect_gt(f$eigen$gap_rel, 1e-8)

    # it must be the number run_eigen wrote, not a recomputation
    h5 <- attr(f, "desc")$filename
    on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
    expect_equal(f$eigen$gap,
                 as.numeric(as.matrix(BigDataStatMeth::hdf5_matrix(
                     h5, "EIGEN/MKsum05/eiggap"))))
    expect_equal(f$eigen$gap_rel,
                 as.numeric(as.matrix(BigDataStatMeth::hdf5_matrix(
                     h5, "EIGEN/MKsum05/eiggap_rel"))))
})

test_that("print() shows the gap", {
    f <- dg_fit(dg_pair(6))
    expect_output(print(f), "spectral gap")
})


# ---- (2) the breakdown warning ---------------------------------------------

test_that("method='solve' warns when the leading subspace is not identified", {
    X <- dg_pair(10)
    expect_warning(f <- dg_fit(X, "solve"),
                   "not numerically identified",
                   class = "mgcca_degenerate_subspace")
    # the warning describes a real degeneracy, it is not a false alarm
    expect_lt(f$eigen$gap_rel, 1e-8)
    # and it is a warning ONLY: the fit is complete and unchanged in shape
    expect_s3_class(f, "mgcca")
    expect_equal(dim(f$Y), c(24L, 2L))
    expect_true(all(is.finite(f$Y)))
})

test_that("the warning names the mechanism and the remedy", {
    X <- dg_pair(10)
    w <- tryCatch(dg_fit(X, "solve"), warning = function(w) w)
    expect_s3_class(w, "mgcca_degenerate_subspace")
    msg <- conditionMessage(w)
    expect_match(msg, "method = \"penalized\"", fixed = TRUE)
    expect_match(msg, "unreliable")
    expect_match(msg, "observed rows")
})

test_that("a healthy solve fit does not warn", {
    expect_no_warning(dg_fit(dg_pair(6), "solve"))
})

test_that("penalized restores the gap on the same fixture and does not warn", {
    X <- dg_pair(10)
    expect_no_warning(f <- dg_fit(X, "penalized", lambda = c(1, 1)))
    expect_gt(f$eigen$gap_rel, 1e-8)
})


# ---- (3) the overlap report -------------------------------------------------

test_that("the overlap report has the right n_j, n_jk and alpha_jk", {
    f <- dg_fit(dg_overlap_blocks(), "solve")
    ov <- f$overlap
    expect_false(is.null(ov))
    expect_equal(ov$n, 30L)
    expect_equal(ov$n_present, c(B1 = 30L, B2 = 24L, B3 = 26L))

    p <- ov$pairs
    expect_equal(nrow(p), 3L)
    key <- paste(p$table1, p$table2, sep = "-")
    expect_equal(p$n_jk[match("B1-B2", key)], 24L)
    expect_equal(p$n_jk[match("B1-B3", key)], 26L)
    expect_equal(p$n_jk[match("B2-B3", key)], 20L)
    expect_equal(p$alpha_jk, p$n_jk / 30)

    expect_equal(ov$alpha_sd, stats::sd(p$alpha_jk))
    expect_equal(ov$alpha_cv, stats::sd(p$alpha_jk) / mean(p$alpha_jk))
    expect_equal(ov$alpha_min, min(p$alpha_jk))
    expect_equal(ov$n_jk_min, 20L)
})

test_that("blocks with identical individuals give alpha = 1 and zero dispersion", {
    set.seed(20260902)
    ids <- sprintf("ind%02d", 1:24)
    X <- list(B1 = matrix(rnorm(24 * 5), 24, 5,
                          dimnames = list(ids, paste0("a", 1:5))),
              B2 = matrix(rnorm(24 * 4), 24, 4,
                          dimnames = list(ids, paste0("b", 1:4))))
    ov <- dg_fit(X, "solve")$overlap
    expect_equal(ov$pairs$n_jk, 24L)
    expect_equal(ov$pairs$alpha_jk, 1)
    expect_equal(ov$alpha_sd, 0)          # a single pair -> defined as 0
    expect_false(ov$heterogeneous)
})

test_that("heterogeneous overlap is flagged, homogeneous overlap is not", {
    # B3 shares only 6 individuals with B1/B2 -> alpha spread far apart
    set.seed(20260902)
    ids <- sprintf("ind%02d", 1:30)
    mk <- function(keep, p, tag) {
        Z <- matrix(rnorm(30 * p), 30, p,
                    dimnames = list(ids, paste0(tag, seq_len(p))))
        Z[keep, , drop = FALSE]
    }
    X <- list(B1 = mk(1:30, 5, "a"), B2 = mk(1:28, 4, "b"), B3 = mk(25:30, 3, "c"))
    ov <- dg_fit(X, "solve")$overlap
    expect_true(ov$heterogeneous)
    expect_true(ov$sparse_pair)                     # min n_jk = 4 or 6, < 30
    expect_false(dg_fit(dg_overlap_blocks(), "solve")$overlap$heterogeneous)
})

test_that("summary() prints the overlap table and its note", {
    f <- dg_fit(dg_overlap_blocks(), "solve")
    expect_output(summary(f, top = 2), "Pairwise overlap")
    expect_output(summary(f, top = 2), "n_jk")
    expect_output(summary(f, top = 2), "alpha")
    expect_output(summary(f, top = 2), "B2-B3")
})


# ---- (4) the diagnostics must not have changed anything --------------------

test_that("the reported quantities are untouched by the diagnostics", {
    X <- dg_overlap_blocks()
    f <- dg_fit(X, "solve")
    # the diagnostics are additive: every pre-existing element is still there
    expect_true(all(c("Y", "corsY", "scores", "pval.cor", "weights",
                      "scaling", "AVE") %in% names(f)))
    expect_true(all(is.finite(f$Y)))
    expect_equal(nrow(f$Y), 30L)
    # and the union size the descriptor reports is the union size the overlap
    # report uses -- one and the same n
    expect_equal(attr(f, "desc")$m, f$overlap$n)
})
