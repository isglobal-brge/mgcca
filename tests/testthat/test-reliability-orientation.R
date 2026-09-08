# HDF5-vs-R ORIENTATION CONTRACT for the reliability layer.
#
# WHY THIS TEST EXISTS. HDF5 stores transposed with respect to R's view, and the
# package's own standards call this "the number one source of bugs -- watch it in
# every product". The reliability layer builds participant Grams out-of-core, so
# every one of its results depends on this contract holding. Reasoning about it
# from code comments is not enough; this measures it.
#
# The fixture is deliberately NON-SQUARE (7 individuals x 4 variables) so that a
# transposition error cannot hide behind a matching dimension, and the negative
# control asserts that the test can tell the two apart at all.

test_that("the R view survives the HDF5 round trip", {
    skip_if_not_installed("BigDataStatMeth")
    set.seed(11)
    n <- 7L; p <- 4L
    X <- matrix(rnorm(n * p), n, p,
                dimnames = list(sprintf("i%d", seq_len(n)), sprintf("v%d", seq_len(p))))
    h5 <- tempfile(fileext = ".h5")
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    hm <- BigDataStatMeth::hdf5_create_matrix(h5, "IN/blk", data = X, overwrite = TRUE)

    expect_identical(dim(hm), c(n, p))
    expect_equal(as.matrix(hm), X, tolerance = 1e-12)
})

test_that("the out-of-core participant Gram equals the in-memory one", {
    skip_if_not_installed("BigDataStatMeth")
    set.seed(11)
    n <- 7L; p <- 4L
    X <- matrix(rnorm(n * p), n, p,
                dimnames = list(sprintf("i%d", seq_len(n)), sprintf("v%d", seq_len(p))))
    h5 <- tempfile(fileext = ".h5")
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    hm <- BigDataStatMeth::hdf5_create_matrix(h5, "IN/blk", data = X, overwrite = TRUE)

    G <- as.matrix(BigDataStatMeth::tcrossprod(hm))
    expect_identical(dim(G), c(n, n))                 # individuals x individuals
    expect_equal(G, tcrossprod(X), tolerance = 1e-10, ignore_attr = TRUE)

    # Negative control: the test must be able to tell the orientations apart.
    # crossprod(X) is p x p, so a silently transposed pipeline could not pass above.
    expect_false(identical(dim(crossprod(X)), dim(tcrossprod(X))))
})

test_that("standardise-then-Gram agrees between HDF5 and R", {
    skip_if_not_installed("BigDataStatMeth")
    set.seed(11)
    n <- 7L; p <- 4L
    X <- matrix(rnorm(n * p), n, p,
                dimnames = list(sprintf("i%d", seq_len(n)), sprintf("v%d", seq_len(p))))
    h5 <- tempfile(fileext = ".h5")
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    hm <- BigDataStatMeth::hdf5_create_matrix(h5, "IN/blk", data = X, overwrite = TRUE)

    # This is the exact chain the reliability layer will use to build a block Gram:
    # standardise the block (present-only, because a stored input block contains
    # only its own individuals), then tcrossprod out-of-core.
    G_h5 <- as.matrix(BigDataStatMeth::tcrossprod(
        BigDataStatMeth::scale(hm, center = TRUE, scale = TRUE)))
    G_r  <- tcrossprod(scale(X, center = TRUE, scale = TRUE))
    expect_equal(G_h5, G_r, tolerance = 1e-8, ignore_attr = TRUE)
})
