# Semantics of individuals that are ABSENT from a block.
#
# mgcca aligns the tables on the union of individual identifiers and pads the
# missing rows of X_j with zeros, so that K_j X_j = X_j.  Those padded rows are
# an algebraic device, not data.  They must not leak into
#   (a) the standard deviation that normalises the block scores      [D-A],
#   (b) the degrees of freedom of the variable/component p-values    [D-B],
#   (c) the block score matrix returned to the user                  [D-C].
# The consensus Y is a different matter: it IS defined on the union (built by
# K-weighting the observed blocks) and must not change.

# ---- fixture: 3 blocks over a union of 30 individuals; B2 misses 6, B3 misses 4
mr_blocks <- function() {
    set.seed(20260901)
    ids <- sprintf("ind%02d", 1:30)
    B1 <- matrix(rnorm(30 * 5), 30, 5, dimnames = list(ids, paste0("v", 1:5)))
    B2 <- matrix(rnorm(30 * 4), 30, 4, dimnames = list(ids, paste0("v", 1:4)))
    B3 <- matrix(rnorm(30 * 3), 30, 3, dimnames = list(ids, paste0("v", 1:3)))
    list(B1 = B1,
         B2 = B2[-c(1, 4, 9, 16, 25, 30), , drop = FALSE],
         B3 = B3[-c(2, 5, 11, 22), , drop = FALSE])
}

# One fit reused by every block below (fitting is the expensive part).
.mr_cache <- new.env(parent = emptyenv())
mr_fit <- function() {
    if (is.null(.mr_cache$fit)) {
        .mr_cache$fit <- mgcca(mr_blocks(), filename = tempfile(fileext = ".h5"),
                               nfac = 2, method = "solve", scores = TRUE,
                               collect = TRUE)
        try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
    }
    .mr_cache$fit
}

# Two-sided p-value of a Pearson correlation on n observations (cor.test's).
mr_pval <- function(r, n)
    2 * stats::pt(abs(r * sqrt(n - 2) / sqrt(1 - r^2)), n - 2, lower.tail = FALSE)


test_that("block scores are standardised over the individuals PRESENT in the block", {
    # D-A: As = A / sd(X_j A). The sd must be taken over the rows of block j
    # that hold data; taken over the whole union it is deflated by the zeros,
    # by a factor sqrt((n_present - 1) / (n_union - 1)) when the block scores
    # are centred, so the returned scores are systematically too large.
    fit <- mr_fit()
    X   <- mr_blocks()
    for (j in names(X)) {
        pres <- rownames(X[[j]])
        S    <- as.matrix(fit$scores[[j]])[pres, , drop = FALSE]
        expect_equal(unname(apply(S, 2, stats::sd)), rep(1, ncol(S)),
                     tolerance = 1e-10, info = j)
    }
})


test_that("p-values use each block's own number of present individuals", {
    # D-B: corsY_j is already computed on the present individuals only, so its
    # p-value must carry n_present - 2 degrees of freedom, not n_union - 2.
    # Using the union inflates the df and makes the test anticonservative.
    fit <- mr_fit()
    X   <- mr_blocks()
    for (j in names(X)) {
        npr <- nrow(X[[j]])
        r   <- as.matrix(fit$corsY[[j]])
        expect_equal(unname(as.matrix(fit$pval.cor[[j]])), unname(mr_pval(r, npr)),
                     tolerance = 1e-10, info = j)
    }
    # and the direction of the old defect: union df gives smaller p-values
    m <- nrow(fit$Y)
    r <- as.matrix(fit$corsY[["B2"]])
    expect_true(all(mr_pval(r, m) <= mr_pval(r, nrow(X$B2)) + 1e-15))
})


test_that("individuals absent from a block get NA in that block's scores", {
    # D-C: an exact 0 is indistinguishable from a genuinely average measured
    # individual. Absent individuals must be marked, not silently placed at the
    # origin. Rows, row names and column names still span the union.
    fit  <- mr_fit()
    X    <- mr_blocks()
    ids  <- rownames(fit$Y)
    for (j in names(X)) {
        S    <- as.matrix(fit$scores[[j]])
        pres <- rownames(X[[j]])
        abs_ <- setdiff(ids, pres)
        expect_identical(rownames(S), ids, info = j)
        expect_identical(colnames(S), colnames(fit$Y), info = j)
        expect_false(anyNA(S[pres, , drop = FALSE]), info = j)
        if (length(abs_))
            expect_true(all(is.na(S[abs_, , drop = FALSE])), info = j)
    }
    expect_equal(length(setdiff(ids, rownames(X$B2))), 6L)
})


test_that("the present block scores are the block's own projection", {
    # The values kept for present individuals are exactly X_j (centred/scaled
    # with the training parameters) times the block weights -- computed here
    # from the fixture and the fit's own weights/scaling, not from the scores.
    fit <- mr_fit()
    X   <- mr_blocks()
    for (j in names(X)) {
        W    <- as.matrix(fit$weights[[j]])
        cs   <- as.matrix(fit$scaling[[j]])[rownames(W), , drop = FALSE]
        sc   <- cs[, "scale"]; sc[sc == 0] <- 1
        Xs   <- sweep(sweep(X[[j]][, rownames(W), drop = FALSE], 2,
                            cs[, "center"], "-"), 2, sc, "/")
        pres <- rownames(X[[j]])
        expect_equal(unname(Xs %*% W),
                     unname(as.matrix(fit$scores[[j]])[pres, , drop = FALSE]),
                     tolerance = 1e-8, info = j)
    }
})


test_that("the consensus Y still scores every individual of the union", {
    # Y is built by K-weighting the observed blocks: it is defined for everyone
    # and must not acquire NAs, nor depend on whether scores were requested.
    fit <- mr_fit()
    X   <- mr_blocks()
    ids <- Reduce(union, lapply(X, rownames))
    expect_equal(nrow(fit$Y), length(ids))
    expect_false(anyNA(fit$Y))

    noscores <- mgcca(X, filename = tempfile(fileext = ".h5"), nfac = 2,
                      method = "solve", scores = FALSE, collect = TRUE)
    on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
    expect_equal(unname(as.matrix(noscores$Y)), unname(as.matrix(fit$Y)),
                 tolerance = 1e-12)
})


test_that("plotScores leaves the absent individuals out of the plot", {
    skip_if_not_installed("ggplot2")
    fit <- mr_fit()
    expect_message(p <- plotScores(fit, table = "B2"), "absent")
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 24L)                 # 30 - 6 absent
    expect_silent(plotScores(fit, table = "B1"))    # complete block: no message
})
