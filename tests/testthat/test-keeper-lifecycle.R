# The handle `mgcca_sensitivity()` and `mgcca_stability()` hold open on their
# HDF5 file for the length of the analysis: who may open the file while it is
# held, and who may not.
#
# A file open through one statically linked HDF5 instance cannot be opened
# through the other, so the four R-level BigDataStatMeth opens of this package
# refuse rather than let "(checkHDF5File) Cannot open file: H5Fopen failed"
# through. BigDataStatMeth's own close-all route is the other side of the same
# coin: it sweeps ITS instance, so a handle held in mgcca's must survive it.

h5_cleanup_on_exit()

# A held window, armed in the caller's frame in the order the exported analyses
# use: the release first, the handle second, so no point of the window owns an
# open file whose release is not yet scheduled.
hold_file <- function(file, frame = parent.frame()) {
    ctx <- new.env(parent = emptyenv())
    ctx$file <- normalizePath(file, mustWork = FALSE)
    withr::defer(mgcca:::.mgcca_rel_release(ctx), envir = frame)
    ctx$handle <- mgcca:::.mgcca_rel_hold(ctx$file)
    ctx
}

.keeper_cache <- new.env(parent = emptyenv())

# One small file-backed fit, shared by every test here.
keeper_fixture <- function() {
    if (!is.null(.keeper_cache$fx)) return(.keeper_cache$fx)
    set.seed(20261002)
    n <- 24L; L <- 2L
    ids <- sprintf("id%03d", seq_len(n))
    p <- c(blkA = 5L, blkB = 4L, blkC = 6L)
    Z <- matrix(stats::rnorm(n * L), n, L)
    tabs <- lapply(seq_along(p), function(j) {
        X <- Z %*% matrix(stats::rnorm(L * p[j]), L, p[j]) +
            matrix(stats::rnorm(n * p[j]), n, p[j])
        dimnames(X) <- list(ids, sprintf("%s_v%02d", names(p)[j],
                                         seq_len(p[j])))
        X
    })
    names(tabs) <- names(p)
    tabs$blkB <- tabs$blkB[-seq.int(n - 3L, n), , drop = FALSE]
    h5 <- h5_tmp()
    fit <- suppressMessages(mgcca(tabs, filename = h5, nfac = L, scale = TRUE,
                                  method = "solve", scores = TRUE))
    BigDataStatMeth::hdf5_close_all()
    .keeper_cache$fx <- list(file = h5, fit = fit, ids = ids,
                             datasets = names(p), L = L,
                             query = stats::setNames(as.numeric(Z[, 1]), ids))
    .keeper_cache$fx
}

test_that("the presence reader refuses a file an mgcca handle holds open", {
    fx <- keeper_fixture()
    hold_file(fx$file)
    expect_error(mgcca:::.mgcca_presence(fx$file, fx$datasets),
                 class = "mgcca_file_held")
    expect_error(mgcca:::.mgcca_presence(fx$file, fx$datasets),
                 "an mgcca HDF5 handle is still open on this file")
})

test_that("the eigenvalue report refuses a file an mgcca handle holds open", {
    fx <- keeper_fixture()
    hold_file(fx$file)
    # Not turned into a missing diagnostic: a held file is the caller's
    # mistake, and the NULL an absent dataset gets would hide it.
    expect_error(mgcca:::.mgcca_eigen_report(fx$file),
                 class = "mgcca_file_held")
})

test_that("mgcca_results() refuses a file an mgcca handle holds open", {
    fx <- keeper_fixture()
    hold_file(fx$file)
    expect_error(mgcca_results(attr(fx$fit, "desc")), class = "mgcca_file_held")
})

test_that("mgcca_import_hdf5() refuses a file an mgcca handle holds open", {
    fx <- keeper_fixture()
    hold_file(fx$file)
    # Listing the datasets of an existing file is one open; writing a table
    # into it is the other.
    expect_error(mgcca_import_hdf5(fx$file, filename = fx$file,
                                   group = "MGCCA_IN"),
                 class = "mgcca_file_held")
    one <- matrix(1, 3L, 2L, dimnames = list(fx$ids[1:3], c("a", "b")))
    expect_error(mgcca_import_hdf5(list(t1 = one, t2 = one),
                                   filename = fx$file, group = "KEEPER_IN"),
                 class = "mgcca_file_held")
})

test_that("a file the handle does not hold is read as usual", {
    fx <- keeper_fixture()
    other <- h5_tmp()
    one <- matrix(c(1, 2, 3, 4, 5, 6), 3L, 2L,
                  dimnames = list(c("a", "b", "c"), c("v1", "v2")))
    desc <- mgcca_import_hdf5(list(t1 = one, t2 = one), filename = other,
                              group = "MGCCA_IN", overwriteFile = TRUE)
    hold_file(fx$file)
    expect_identical(sort(desc$datasets), c("t1", "t2"))
    again <- mgcca_import_hdf5(other, filename = other, group = "MGCCA_IN")
    expect_identical(sort(again$datasets), c("t1", "t2"))
})

test_that("only one handle is held at a time", {
    fx <- keeper_fixture()
    hold_file(fx$file)
    expect_error(mgcca:::.mgcca_rel_hold(fx$file),
                 class = "mgcca_nested_file_handle")
    expect_error(suppressMessages(mgcca_sensitivity(fx$fit, query = fx$query)),
                 class = "mgcca_nested_file_handle")
})

test_that("a held handle survives BigDataStatMeth's close-all route", {
    fx <- keeper_fixture()
    ctx <- hold_file(fx$file)
    expect_true(mgcca:::mgcca_file_is_open_rcpp(ctx$handle))
    # close_all sweeps BigDataStatMeth's own HDF5 instance; the handle lives in
    # mgcca's, which is the whole reason the reverse collision exists.
    suppressMessages(BigDataStatMeth::hdf5_close_all())
    expect_true(mgcca:::mgcca_file_is_open_rcpp(ctx$handle))
    expect_identical(mgcca:::.mgcca_held_file$path, ctx$file)
    d <- mgcca:::mgcca_read_dimensions_rcpp(ctx$file, "MGCCA_IN", fx$datasets)
    expect_identical(as.integer(d$nrow[["blkA"]]), 24L)
    mgcca:::.mgcca_rel_release(ctx)
    expect_false(mgcca:::mgcca_file_is_open_rcpp(ctx$handle))
    expect_null(mgcca:::.mgcca_held_file$path)
})

test_that("releasing twice, and releasing nothing, are both no-ops", {
    fx <- keeper_fixture()
    ctx <- hold_file(fx$file)
    mgcca:::.mgcca_rel_release(ctx)
    expect_silent(mgcca:::.mgcca_rel_release(ctx))
    expect_silent(mgcca:::.mgcca_rel_release(NULL))
    expect_null(mgcca:::.mgcca_held_file$path)
    empty <- new.env(parent = emptyenv())
    expect_silent(mgcca:::.mgcca_rel_release(empty))
})

test_that("the analyses leave nothing held, on return and after an error", {
    fx <- keeper_fixture()
    on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
    sens <- suppressMessages(mgcca_sensitivity(fx$fit, query = fx$query))
    expect_identical(nrow(sens$by_block), 3L)
    expect_null(mgcca:::.mgcca_held_file$path)
    stab <- suppressMessages(mgcca_stability(fx$fit, method = "subsample",
                                             B = 4L))
    expect_identical(nrow(stab$resamples), 4L)
    expect_null(mgcca:::.mgcca_held_file$path)
    # A query a covariate spans has no alignment: the analysis stops inside the
    # held window.
    cv <- matrix(fx$query, ncol = 1L, dimnames = list(fx$ids, "same"))
    expect_error(suppressMessages(
        mgcca_sensitivity(fx$fit, query = fx$query, covariates = cv)))
    expect_null(mgcca:::.mgcca_held_file$path)
    expect_true(is.matrix(as.matrix(
        BigDataStatMeth::hdf5_matrix(fx$file, "FINAL_RESULTS/Y"))))
})
