# The C++ entry points back to back on one file, in the shape Windows refuses:
# read dimensions, read blocks, save a group, read row names, over and over.
# Inside the window of a held handle every one of them must succeed and return
# exactly what the same call returns with nothing held -- the handle is a
# lifetime, not a different way of reading.
#
# Then the window is left through an error rather than a return, and the release
# is checked three ways: the registry is empty, the handle no longer holds the
# file, and BigDataStatMeth can open the same file again.

h5_cleanup_on_exit()

churn_file <- function() {
    set.seed(20261002)
    n <- 20L; L <- 2L
    ids <- sprintf("id%03d", seq_len(n))
    p <- c(blkA = 4L, blkB = 3L)
    Z <- matrix(stats::rnorm(n * L), n, L)
    tabs <- lapply(seq_along(p), function(j) {
        X <- Z %*% matrix(stats::rnorm(L * p[j]), L, p[j]) +
            matrix(stats::rnorm(n * p[j]), n, p[j])
        dimnames(X) <- list(ids, sprintf("%s_v%02d", names(p)[j],
                                         seq_len(p[j])))
        X
    })
    names(tabs) <- names(p)
    h5 <- h5_tmp()
    suppressMessages(mgcca(tabs, filename = h5, nfac = L, scale = TRUE,
                           method = "solve", scores = TRUE))
    BigDataStatMeth::hdf5_close_all()
    list(file = normalizePath(h5), datasets = names(p), ids = ids)
}

test_that("the C++ entry points churn unchanged inside a held window", {
    fx <- churn_file()
    ds <- fx$datasets
    # Reference values, taken with nothing held: one open per call, the way
    # every version before this one worked.
    ref_dim  <- mgcca:::mgcca_read_dimensions_rcpp(fx$file, "MGCCA_IN", ds)
    ref_blk  <- mgcca:::mgcca_read_blocks_rcpp(fx$file, "MGCCA_IN", ds)
    ref_rown <- mgcca:::mgcca_read_rownames_rcpp(fx$file, "MGCCA_IN", ds)

    ctx <- new.env(parent = emptyenv())
    ctx$file <- fx$file
    on.exit(mgcca:::.mgcca_rel_release(ctx), add = TRUE)
    ctx$handle <- mgcca:::.mgcca_rel_hold(ctx$file)

    for (round in 1:3) {
        grp <- paste0("KEEPER_CHURN/round", round)
        expect_identical(
            mgcca:::mgcca_read_dimensions_rcpp(ctx$file, "MGCCA_IN", ds),
            ref_dim)
        expect_identical(
            mgcca:::mgcca_read_blocks_rcpp(ctx$file, "MGCCA_IN", ds),
            ref_blk)
        # A write between the reads, into a group of its own each round: the
        # writer refuses to overwrite, so repeating it needs a new group.
        mgcca:::mgcca_save_audit_rcpp(ctx$file, grp,
                                      stats::setNames(ref_blk$values, ds),
                                      stats::setNames(list(), character(0)))
        expect_identical(
            mgcca:::mgcca_read_rownames_rcpp(ctx$file, "MGCCA_IN", ds),
            ref_rown)
        # A read-only entry point in the same window: with a handle live the
        # file is opened read-write whatever the entry point asked for, which
        # must change nothing about what it reports.
        expect_identical(sort(mgcca:::mgcca_list_group_rcpp(ctx$file, grp)),
                         sort(ds))
        # And the written group reads back as what went in.
        expect_identical(
            mgcca:::mgcca_read_blocks_rcpp(ctx$file, grp, ds)$values,
            ref_blk$values)
    }
    expect_true(mgcca:::mgcca_file_is_open_rcpp(ctx$handle))
    expect_identical(mgcca:::.mgcca_held_file$path, ctx$file)
})

test_that("an error inside the window still gives the handle up", {
    fx <- churn_file()
    kept <- NULL
    # The window of a function that stops halfway, armed exactly as the
    # exported analyses arm theirs.
    broken <- function(file) {
        ctx <- new.env(parent = emptyenv())
        ctx$file <- file
        on.exit(mgcca:::.mgcca_rel_release(ctx), add = TRUE)
        ctx$handle <- mgcca:::.mgcca_rel_hold(ctx$file)
        kept <<- ctx$handle
        mgcca:::mgcca_read_dimensions_rcpp(ctx$file, "MGCCA_IN", fx$datasets)
        stop("stopped halfway, on purpose")
    }
    expect_error(broken(fx$file), "stopped halfway, on purpose")

    # (1) the registry is empty
    expect_null(mgcca:::.mgcca_held_file$path)
    expect_null(mgcca:::.mgcca_held_file$handle)
    # (2) the handle no longer holds the file
    expect_false(mgcca:::mgcca_file_is_open_rcpp(kept))
    # (3) BigDataStatMeth opens the same file again
    hm <- BigDataStatMeth::hdf5_matrix(fx$file, "FINAL_RESULTS/Y")
    expect_identical(nrow(as.matrix(hm)), 20L)
    hm$close()
    suppressMessages(BigDataStatMeth::hdf5_close_all())
})

test_that("a second analysis on the same file runs from zero", {
    fx <- churn_file()
    ds <- fx$datasets
    ref <- mgcca:::mgcca_read_blocks_rcpp(fx$file, "MGCCA_IN", ds)
    # One whole window, returning the handle it held so the caller can see it
    # was given up on the way out.
    one_window <- function(file) {
        ctx <- new.env(parent = emptyenv())
        ctx$file <- file
        on.exit(mgcca:::.mgcca_rel_release(ctx), add = TRUE)
        ctx$handle <- mgcca:::.mgcca_rel_hold(ctx$file)
        expect_identical(
            mgcca:::mgcca_read_blocks_rcpp(ctx$file, "MGCCA_IN", ds), ref)
        ctx$handle
    }
    for (round in 1:2) {
        handle <- one_window(fx$file)
        expect_null(mgcca:::.mgcca_held_file$path)
        expect_false(mgcca:::mgcca_file_is_open_rcpp(handle))
    }
})
