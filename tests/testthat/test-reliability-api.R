# End-to-end tests for the public reliability API.
#
# Self-contained and deterministic: small synthetic blocks with a fixed seed, no
# external dataset and no download. The blocks deliberately have partly
# overlapping individuals, which is the setting the estimator exists for and the
# one where the availability weighting matters.

rel_api_fixture <- function(n = 40L, L = 2L, seed = 3L) {
    set.seed(seed)
    ids <- sprintf("s%02d", seq_len(n))
    Z <- matrix(rnorm(n * L), n, L)
    mk <- function(p, keep) {
        X <- Z %*% matrix(rnorm(L * p), L, p) + matrix(rnorm(n * p), n, p)
        rownames(X) <- ids
        X[keep, , drop = FALSE]
    }
    h5 <- tempfile(fileext = ".h5")
    fit <- mgcca(list(a = mk(8L, seq_len(n)), b = mk(6L, 3:n)),
                 filename = h5, nfac = L, method = "penalized", lambda = c(1, 1))
    grp <- factor(rep(c("x", "y"), length.out = n)); names(grp) <- ids
    list(fit = fit, ids = ids, h5 = h5, Z = Z, L = L, n = n, group = grp,
         query = setNames(as.numeric(Z[, 1]) + rnorm(n, 0, 0.4), ids))
}

test_that("a fit that cannot locate its input blocks is refused, with the reason", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    expect_error(mgcca_sensitivity(structure(list(), class = "mgcca"), f$query),
                 "provenance descriptor")
    bad <- f$fit; attr(bad, "desc")$input_group <- NA_character_
    expect_error(mgcca_sensitivity(bad, f$query), "INPUT blocks")
})

test_that("participants are aligned by identifier and never by position", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    expect_error(mgcca_sensitivity(f$fit, unname(f$query)), "identifiers")
    q <- f$query; names(q)[2] <- names(q)[1]
    expect_error(mgcca_sensitivity(f$fit, q), "duplicated")
    Q <- cbind(one = f$query, two = f$query)
    colnames(Q) <- c("a", "a")
    expect_error(mgcca_sensitivity(f$fit, Q), "duplicated column names")
})

test_that("the memory and HDF5 backends agree, including with forced blocking", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    m <- mgcca_sensitivity(f$fit, f$query, group = f$group, backend = "memory")
    h <- mgcca_sensitivity(f$fit, f$query, group = f$group, backend = "hdf5")
    b <- mgcca_sensitivity(f$fit, f$query, group = f$group, backend = "hdf5",
                           block_size = 3L)
    expect_equal(h$overall$T, m$overall$T, tolerance = 1e-8)
    expect_equal(b$overall$T, m$overall$T, tolerance = 1e-8)
    expect_equal(b$overall$S_total, m$overall$S_total, tolerance = 1e-6)
    expect_identical(m$settings$backend, "memory")
    expect_identical(h$settings$backend, "hdf5")
})

test_that("the grouped and ungrouped results differ only in the group-dependent fields", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    u <- mgcca_sensitivity(f$fit, f$query)
    g <- mgcca_sensitivity(f$fit, f$query, group = f$group)
    expect_equal(u$overall$T, g$overall$T, tolerance = 1e-10)
    expect_equal(u$overall$S_total, g$overall$S_total, tolerance = 1e-10)
    # ABSENT MEANS ABSENT: without a grouping the columns are not present at all.
    # An earlier version filled them with NA and this test asserted that; a column
    # of NAs is a column that should not be there, so both the code and the
    # expectation changed together.
    expect_false("S_between" %in% names(u$overall))
    expect_false("R_between" %in% names(u$overall))
    expect_false("S_between" %in% names(u$by_block))
    expect_true("S_between" %in% names(g$overall))
    expect_false(is.na(g$overall$S_between))
    expect_gte(g$overall$R_between, 0); expect_lte(g$overall$R_between, 1)
})

test_that("the alignment is bounded and discriminates signal from noise", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    Q <- cbind(signal = f$query, noise = setNames(rnorm(f$n), f$ids))
    rownames(Q) <- f$ids
    s <- mgcca_sensitivity(f$fit, Q, group = f$group)
    expect_true(all(s$overall$T >= 0 & s$overall$T <= 1))
    expect_gt(s$overall$T[s$overall$query == "signal"],
              s$overall$T[s$overall$query == "noise"])
    # the per-block decomposition adds up to the total it decomposes
    for (q in s$overall$query) {
        bb <- s$by_block[s$by_block$query == q, ]
        expect_equal(sum(bb$S_total), s$overall$S_total[s$overall$query == q],
                     tolerance = 1e-8)
        expect_equal(sum(bb$share_total), 1, tolerance = 1e-8)
    }
})

test_that("a query spanned by its covariates is an error, not a number from rounding", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    cv <- matrix(f$query, ncol = 1, dimnames = list(f$ids, "same"))
    expect_error(mgcca_sensitivity(f$fit, f$query, covariates = cv),
                 "no variation left")
})

test_that("the block profile extracts and never recomputes", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s <- mgcca_sensitivity(f$fit, f$query, group = f$group)
    p <- mgcca_block_profile(s)
    expect_s3_class(p, "mgcca_block_profile")
    expect_identical(p$profile, s$by_block)          # the same object, reformatted
    expect_error(mgcca_block_profile(s, query = "nope"), "unknown query")
    expect_error(mgcca_block_profile(f$fit), "mgcca_sensitivity")
})

test_that("stability runs both arms, bounds its metrics and archives its plan", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    expect_error(mgcca_stability(f$fit, method = "loco"), "`group` is required")
    s <- mgcca_stability(f$fit, method = "both", group = f$group, B = 6L, seed = 1L)
    ok <- s$resamples$ok
    expect_true(all(s$resamples$overlap[ok] <= 1 + 1e-8))
    expect_true(all(s$resamples$rank_margin[ok] <= 1 + 1e-8))
    expect_setequal(unique(s$resamples$method), c("subsample", "loco"))
    expect_setequal(s$resamples$left_out[s$resamples$method == "loco"], levels(f$group))
    expect_equal(length(s$plan), nrow(s$resamples))
    # the same seed reproduces the run
    a <- mgcca_stability(f$fit, B = 5L, seed = 42L)$resamples$overlap
    b <- mgcca_stability(f$fit, B = 5L, seed = 42L)$resamples$overlap
    expect_equal(a, b)
})

test_that("persistence is opt-in and refuses to overwrite", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s0 <- mgcca_sensitivity(f$fit, f$query, group = f$group)
    expect_null(s0$storage)
    s1 <- mgcca_sensitivity(f$fit, f$query, group = f$group, store = TRUE)
    expect_identical(s1$storage$group, "RELIABILITY/sensitivity")
    expect_error(mgcca_sensitivity(f$fit, f$query, group = f$group, store = TRUE),
                 "already exists")
})

test_that("the S3 methods print and the plots are ggplots", {
    skip_if_not_installed("BigDataStatMeth")
    f <- rel_api_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s <- mgcca_sensitivity(f$fit, f$query, group = f$group)
    expect_output(print(s), "mgcca sensitivity")
    expect_output(summary(s), "Per-block decomposition")
    expect_s3_class(plot(s), "ggplot")
    expect_s3_class(plot(s, type = "overall"), "ggplot")
    expect_s3_class(plot(mgcca_block_profile(s)), "ggplot")
    st <- mgcca_stability(f$fit, B = 4L, seed = 2L)
    expect_output(print(st), "subspace stability")
    expect_s3_class(plot(st), "ggplot")
})
