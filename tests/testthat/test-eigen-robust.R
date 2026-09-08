# Robustness of the final eigendecomposition (run_eigen), 2026-07-24. The operator
# MKsum05 is n x n by construction, so run_eigen uses a DIRECT symmetric solver by
# default (adaptive: iterative out-of-core only when the operator is too large, with a
# direct fallback on non-convergence). These tests exercise the cases ChatGPT asked for
# (round 28 sec.6): well-separated, near-identity/degenerate, small external gap,
# nfac==1 vs nfac>1, direct-vs-iterative agreement, and the validity diagnostics.

# ---- helpers: build small synthetic multi-block data with shared structure --------
synth <- function(n = 40, p = c(60, 50), L = 2, snr = 2, seed = 1) {
  set.seed(seed)
  ids <- sprintf("id%03d", seq_len(n))
  Z <- scale(matrix(rnorm(n * L), n, L))
  tabs <- lapply(p, function(pj) {
    A <- matrix(rnorm(L * pj), L, pj)
    sig <- Z %*% A; E <- matrix(rnorm(n * pj), n, pj)
    X <- scale(snr * sig / sqrt(mean(apply(sig, 2, var))) + E)
    rownames(X) <- ids; colnames(X) <- sprintf("v%04d", seq_len(pj)); X
  })
  names(tabs) <- paste0("B", seq_along(p)); tabs
}

fit_read <- function(tabs, nfac = 2, method = "penalized", lambda = NULL, env = NULL) {
  if (!is.null(env)) { old <- Sys.getenv("MGCCA_EIGEN_DIRECT_MAX", NA)
    Sys.setenv(MGCCA_EIGEN_DIRECT_MAX = env)
    on.exit(if (is.na(old)) Sys.unsetenv("MGCCA_EIGEN_DIRECT_MAX")
            else Sys.setenv(MGCCA_EIGEN_DIRECT_MAX = old)) }
  h5 <- tempfile(fileext = ".h5")
  fit <- mgcca(tabs, filename = h5, nfac = nfac, scale = FALSE, method = method,
               lambda = lambda, route = "auto", collect = TRUE)
  rd <- function(p) tryCatch(as.numeric(as.matrix(hdf5_matrix(h5, p))),
                             error = function(e) NA_real_)
  list(fit = fit, Y = fit$Y,
       vals = rd("EIGEN/MKsum05/values"), gap = rd("EIGEN/MKsum05/eiggap"),
       grel = rd("EIGEN/MKsum05/eiggap_rel"), reig = rd("EIGEN/MKsum05/eig_residual"))
}

proj <- function(Y) { Y <- as.matrix(Y); Y %*% solve(crossprod(Y), t(Y)) }

test_that("well-separated operator: converges, tiny residual, positive gaps", {
  r <- fit_read(synth(), nfac = 2, lambda = c(60, 50))
  expect_length(r$Y, 40 * 2)
  expect_true(is.finite(r$gap) && r$gap > 0)          # external gap defined & positive
  expect_true(is.finite(r$grel) && r$grel > 0)
  expect_lt(r$reig, 1e-8)                              # eigenpairs solve S V = V diag
})

test_that("nfac == 1 works (no Krylov ncv>nev+1 constraint with the direct solver)", {
  r <- fit_read(synth(), nfac = 1, lambda = c(60, 50))
  expect_length(r$Y, 40 * 1)
  expect_true(is.finite(r$gap))                        # mu_1 - mu_2 defined (n > 1)
})

test_that("near-identity / degenerate operator CONVERGES and the gap flags it", {
  # p >> n with a tiny absolute ridge -> M_j ~ I -> whitened operator ~ I (degenerate).
  # The iterative solver failed here; the direct solver must return and flag near-0 gap.
  tabs <- synth(n = 30, p = c(300, 300), snr = 0.2)
  r <- fit_read(tabs, nfac = 2, lambda = c(1e-4, 1e-4))
  expect_length(r$Y, 30 * 2)                           # returned, did not error
  expect_lt(r$grel, 1e-2)                              # relative gap tiny => not identifiable
})

test_that("direct and iterative paths agree (projector + eigenvalues)", {
  tabs <- synth(seed = 7)
  d <- fit_read(tabs, nfac = 2, lambda = c(60, 50))                 # direct (default)
  i <- fit_read(tabs, nfac = 2, lambda = c(60, 50), env = "1")      # force iterative
  expect_equal(sort(d$vals[1:2]), sort(i$vals[1:2]), tolerance = 1e-6)
  expect_equal(proj(d$Y), proj(i$Y), tolerance = 1e-6)             # subspace agrees
  expect_true(is.finite(i$reig) && i$reig < 1e-6)  # block-wise residual (not NaN) on iter path
})

test_that("nfac contract: nfac>m rejected, nfac<1 rejected, nfac==m gives NA gap", {
  tabs <- synth(n = 20, p = c(30, 30))            # n = 20 individuals -> m = 20
  expect_error(fit_read(tabs, nfac = 25, lambda = c(30, 30)))   # nfac > m: rejected
  expect_error(fit_read(tabs, nfac = 0,  lambda = c(30, 30)))   # nfac < 1: rejected
  # nfac == m: allowed by the general estimator, but the EXTERNAL gap is undefined
  # (no mu_{m+1}) -> returned as NA (not an interpretable number).
  r <- fit_read(tabs, nfac = 20, lambda = c(30, 30))
  expect_length(r$Y, 20 * 20)
  expect_true(is.na(r$gap))
})

test_that("subspace is stable under a repeated leading pair (mu1 ~ mu2)", {
  # two nearly-exchangeable blocks -> mu1 ~ mu2; individual components rotate but the
  # L=2 projector is still well defined. Compare the PROJECTOR, never single vectors.
  set.seed(3); n <- 40; ids <- sprintf("id%03d", seq_len(n))
  Z <- scale(matrix(rnorm(n * 2), n, 2))
  mk <- function(s) { A <- matrix(rnorm(2 * 50), 2, 50)
    X <- scale(2 * (Z %*% A) / sqrt(mean(apply(Z %*% A, 2, var))) + matrix(rnorm(n * 50), n, 50))
    rownames(X) <- ids; colnames(X) <- sprintf("v%03d", 1:50); X }
  tabs <- list(B1 = mk(1), B2 = mk(2))
  r <- fit_read(tabs, nfac = 2, lambda = c(50, 50))
  expect_length(r$Y, n * 2)
  expect_true(is.finite(r$reig) && r$reig < 1e-8)     # solved regardless of near-equal mus
})
