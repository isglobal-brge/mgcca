# Unit tests for the internal reliability algebra (R/reliability_algebra.R).
#
# Deterministic and self-contained: a fixed seed, small synthetic blocks, no HDF5
# and no download. Two of the checks are analytical oracles rather than
# regression values -- a query lying inside the retained subspace must give
# T = 1 and one orthogonal to it must give T = 0 -- so the test detects a wrong
# answer, not merely a changed one.

rel_fixture <- function(n = 40L, J = 3L, L = 2L, p = 12L, seed = 1L) {
    set.seed(seed)
    ids <- sprintf("id%02d", seq_len(n))
    Z <- matrix(rnorm(n * L), n, L)
    Glist <- vector("list", J); present <- vector("list", J)
    for (j in seq_len(J)) {
        X <- Z %*% matrix(rnorm(L * p), L, p) + matrix(rnorm(n * p), n, p)
        m <- rep(TRUE, n); m[sample.int(n, 5L)] <- FALSE; X[!m, ] <- 0
        G <- tcrossprod(X); dimnames(G) <- list(ids, ids)
        Glist[[j]] <- G; present[[j]] <- m
    }
    list(fit = mgcca:::.mgcca_rel_fit(Glist, present, rep(1e-3, J), L, ids),
         ids = ids, n = n, L = L)
}

test_that("the participant-space reference fit is well formed", {
    f <- rel_fixture()
    expect_s3_class(f$fit, "mgcca_rel_fit")
    expect_equal(f$fit$n, f$n)
    expect_lt(max(abs(f$fit$S - t(f$fit$S))), 1e-12)
    expect_false(is.unsorted(rev(f$fit$mu)))
    expect_true(all(is.finite(f$fit$D)))
})

test_that("a participant present in no block is a hard error, not an Inf weight", {
    f <- rel_fixture()
    present <- f$fit$present
    for (j in seq_along(present)) present[[j]][1L] <- FALSE
    expect_error(
        mgcca:::.mgcca_rel_fit(f$fit$Glist, present, f$fit$lambda, f$L, f$ids),
        "present in no block")
})

test_that("the alignment functional matches its analytical values", {
    f <- rel_fixture(); fit <- f$fit
    z <- setNames(rnorm(f$n), f$ids)
    Tv <- mgcca:::.mgcca_rel_T(fit, z)
    expect_gte(Tv, 0); expect_lte(Tv, 1)
    # inside the retained subspace -> 1 ; orthogonal to it -> 0
    z_in  <- setNames(as.numeric(fit$D * fit$V[, 1L]),  f$ids)
    z_out <- setNames(as.numeric(fit$D * fit$V[, f$n]), f$ids)
    expect_equal(mgcca:::.mgcca_rel_T(fit, z_in),  1, tolerance = 1e-10)
    expect_equal(mgcca:::.mgcca_rel_T(fit, z_out), 0, tolerance = 1e-10)
})

test_that("participants are aligned by identifier and never by position", {
    f <- rel_fixture(); fit <- f$fit
    z <- setNames(rnorm(f$n), f$ids)
    expect_error(mgcca:::.mgcca_rel_T(fit, unname(z)), "identifiers")
    z_dup <- z; names(z_dup)[2L] <- names(z_dup)[1L]
    expect_error(mgcca:::.mgcca_rel_T(fit, z_dup), "duplicated")
    z_unk <- z; names(z_unk)[3L] <- "not_in_fit"
    expect_error(mgcca:::.mgcca_rel_T(fit, z_unk), "absent from the fit")
    # partial coverage is legitimate, and is counted rather than assumed
    expect_equal(mgcca:::.mgcca_rel_zstar(fit, z[1:30])$n_used, 30L)
    expect_error(mgcca:::.mgcca_rel_zstar(fit, setNames(rep(NA_real_, f$n), f$ids)),
                 "no finite value")
})

test_that("a grouping keeps its labels and rejects unusable levels", {
    f <- rel_fixture()
    g <- setNames(rep(c("a", "b", "c"), length.out = f$n), f$ids)
    G <- mgcca:::.mgcca_rel_group(g, f$ids)
    expect_identical(G$levels, c("a", "b", "c"))
    expect_equal(range(G$code), c(0L, 2L))
    expect_true(min(G$index) >= 0L)
    g_singleton <- g; g_singleton[1L] <- "solo"
    expect_error(mgcca:::.mgcca_rel_group(g_singleton, f$ids), "fewer than 2")
    g_na <- g; g_na[1:5] <- NA
    expect_equal(mgcca:::.mgcca_rel_group(g_na, f$ids)$n_dropped, 5L)
})

test_that("the perturbation coefficients and per-block inputs have the kernel's shapes", {
    f <- rel_fixture(); fit <- f$fit
    z <- setNames(rnorm(f$n), f$ids)
    zs <- mgcca:::.mgcca_rel_zstar(fit, z)$zstar
    Cm <- mgcca:::.mgcca_rel_C(fit, tcrossprod(zs) / sum(zs^2))
    expect_equal(dim(Cm), c(f$L, f$n - f$L))
    expect_true(all(is.finite(Cm)))
    bi <- mgcca:::.mgcca_rel_block_inputs(fit, 1L)
    expect_equal(dim(bi$Rr), c(f$n, f$L))
    expect_equal(dim(bi$Sm), c(f$n, f$n - f$L))
    expect_equal(bi$alpha, fit$lambda[1L])
})

test_that("validity reports the eigenvalue separation the expansion depends on", {
    f <- rel_fixture()
    v <- mgcca:::.mgcca_rel_validity(f$fit)
    expect_true(is.finite(v$rel_gap))
    expect_true(is.logical(v$valid))
    if (isTRUE(v$valid)) expect_true(is.na(v$reason)) else expect_true(nzchar(v$reason))
})
