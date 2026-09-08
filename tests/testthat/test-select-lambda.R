# =============================================================================
# mgcca_select_lambda() -- the CV cross-block agreement guardrail (wave 2).
#
# The criterion under test is the one the thesis calibration establishes in
# reliability/05_lambda_selection.md and implements in
# reliability/calibration/estimator.R:445-476 (cv_crossblock) and :501-671
# (cv_crossblock_gamma, the frozen-contract version over the dimensionless
# ridge). Its content: lambda is chosen by held-out cross-block score agreement,
# NEVER by minimising predicted instability -- in the specific-dominant regime
# the stability-optimal lambda sits at the over-regularised end of the grid and
# costs 0.666 of recovery, while the CV-optimal lambda costs 0.002.
#
# These tests reproduce that guardrail qualitatively on synthetic fixtures whose
# truth is known, and pin the API shape. Everything here is deterministic.
# =============================================================================

## ---- self-contained fixture generator --------------------------------------
## Mirrors MGCCA_Simulations/R/simulate.R: X_j = Z A_j + specstr * F_j B_j + E_j
## with Z the SHARED latent (the truth) and F_j a block-SPECIFIC axis that a
## sufficiently over-regularised operator latches onto instead.
sl_sim <- function(n, p, L, snr, specstr, seed) {
  set.seed(seed)
  J   <- length(p)
  ids <- paste0("ind", formatC(seq_len(n), width = nchar(n), flag = "0"))
  on. <- function(M) qr.Q(qr(M))[, seq_len(ncol(M)), drop = FALSE]
  Z <- base::scale(on.(matrix(rnorm(n * L), n, L)))
  rownames(Z) <- ids
  tabs <- vector("list", J); names(tabs) <- paste0("table", seq_len(J))
  for (j in seq_len(J)) {
    pj  <- p[j]
    sig <- Z %*% matrix(rnorm(L * pj), L, pj)
    Fj  <- base::scale(on.(matrix(rnorm(n * 1L), n, 1L)))
    sig <- sig + specstr * (Fj %*% matrix(rnorm(pj), 1L, pj))
    vs  <- mean(apply(sig, 2, stats::var))
    Xj  <- sig + matrix(rnorm(n * pj, sd = sqrt(vs / snr)), n, pj)
    dimnames(Xj) <- list(ids, paste0("t", j, "_v", seq_len(pj)))
    tabs[[j]] <- Xj
  }
  list(tables = tabs, Z = Z, ids = ids)
}

## Drop individuals from each table; the union still spans everybody.
sl_missing <- function(tables, frac, seed) {
  set.seed(seed)
  ids <- rownames(tables[[1]]); n <- length(ids)
  keep <- lapply(tables, function(m) sort(sample(n, round(n * (1 - frac)))))
  orph <- setdiff(seq_len(n), Reduce(union, keep))
  for (i in orph) { j <- sample(length(keep), 1L)
                    keep[[j]] <- sort(union(keep[[j]], i)) }
  Map(function(m, k) m[k, , drop = FALSE], tables, keep)
}

## RV coefficient (rotation/sign/scale invariant), as in
## MGCCA_Simulations/R/metrics.R:23-30 -- used here only to score recovery.
sl_rv <- function(A, B) {
  A <- base::scale(as.matrix(A), TRUE, FALSE)
  B <- base::scale(as.matrix(B), TRUE, FALSE)
  SA <- base::tcrossprod(A); SB <- base::tcrossprod(B)
  sum(SA * SB) / sqrt(sum(SA * SA) * sum(SB * SB))
}

## The decisive regime of 05_lambda_selection.md C2: specific-dominant.
sl_c2 <- function() {
  s <- sl_sim(n = 60, p = c(120, 90), L = 2, snr = 1.5, specstr = 3, seed = 4242)
  list(tables = sl_missing(s$tables, 0.3, 707), Z = s$Z)
}

sl_grid <- c(1, 10, 100, 1000, 1e4, 1e5, 1e6)

## In this regime the CV curve falls monotonically in lambda over the whole grid,
## so the best candidate sits on the LOWER edge and mgcca_select_lambda warns --
## correctly -- about a boundary optimum. It is the same behaviour
## 05_lambda_selection.md records as its honest caveat for regime C2 (argmax-CV
## at the smallest grid point). Expected here, so it is suppressed at the call
## and asserted explicitly in the guardrail test below.
sl_select <- function(...) suppressWarnings(mgcca_select_lambda(...))


## ---- API shape --------------------------------------------------------------

test_that("mgcca_select_lambda is exported and returns the documented object", {
  expect_true(is.function(mgcca_select_lambda))
  expect_true("mgcca_select_lambda" %in% getNamespaceExports("mgcca"))

  f  <- sl_c2()
  sel <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                   scale_matched = FALSE, K = 4, seed = 11)

  expect_s3_class(sel, "mgcca_lambda")
  # the selected lambda: one per table, named, positive, on the grid
  expect_type(sel$lambda, "double")
  expect_length(sel$lambda, length(f$tables))
  expect_equal(names(sel$lambda), names(f$tables))
  expect_true(all(sel$lambda > 0))
  expect_true(all(sel$lambda %in% sl_grid))
  # the agreement diagnostic
  expect_true(is.numeric(sel$agreement) && length(sel$agreement) == 1L)
  expect_true(is.finite(sel$agreement) && sel$agreement >= -1 && sel$agreement <= 1)
  # the candidate grid trace
  expect_s3_class(sel$grid, "data.frame")
  expect_true(all(c("lambda", "cv_mean", "cv_se", "n_cells", "in_1se", "selected")
                  %in% names(sel$grid)))
  expect_equal(nrow(sel$grid), length(sl_grid))
  expect_equal(sum(sel$grid$selected), 1L)
  expect_equal(sel$agreement, sel$grid$cv_mean[sel$grid$selected])
  # the record of how the choice was made
  expect_identical(sel$criterion, "cv_crossblock_agreement")
  expect_identical(sel$rule, "argmax")
  expect_identical(sel$nfac, 2L)
  expect_identical(sel$K, 4L)
  expect_identical(sel$seed, 11)
  expect_true(is.character(sel$status) && length(sel$status) == 1L)
  expect_true(!is.null(sel$call))
})

test_that("the selection is deterministic in the seed", {
  f <- sl_c2()
  a <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                 scale_matched = FALSE, K = 4, seed = 11)
  b <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                 scale_matched = FALSE, K = 4, seed = 11)
  a$call <- b$call <- NULL
  expect_identical(a, b)
  c2 <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                  scale_matched = FALSE, K = 4, seed = 12)
  expect_false(identical(a$fold, c2$fold))
})

test_that("print() shows the choice and the criterion", {
  f <- sl_c2()
  sel <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                   scale_matched = FALSE, K = 4, seed = 11)
  out <- paste(capture.output(print(sel)), collapse = "\n")
  expect_match(out, "cross-block")
  expect_match(out, "agreement")
  expect_match(out, "lambda")
})


## ---- THE GUARDRAIL ----------------------------------------------------------
## 05_lambda_selection.md section 2: in the specific-dominant regime the
## stability-optimal lambda is the most regularised point of the grid, and it is
## catastrophic for recovery. The CV criterion must not go there.

test_that("CV agreement rejects the over-regularised (stability-optimal) end", {
  f <- sl_c2()
  sel <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                   scale_matched = FALSE, K = 4, seed = 11)
  g <- sel$grid
  lam_max <- max(sl_grid)

  # (i) the trap is NOT chosen
  expect_false(isTRUE(all(sel$lambda == lam_max)))
  expect_false(g$selected[which.max(g$lambda)])

  # (ii) CV collapses at over-regularisation -- this is WHY it is not chosen
  cv_sel <- sel$agreement
  cv_top <- g$cv_mean[which.max(g$lambda)]
  expect_true(cv_top < cv_sel - 0.05)

  # (iii) and the trap is outside the 1-SE near-optimal region as well
  expect_false(g$in_1se[which.max(g$lambda)])

  # (iv) the criterion falls all the way to the lower edge here, so the honest
  # report is a boundary optimum -- and that is what comes back.
  expect_identical(sel$status, "boundary_optimum")
})

test_that("the CV-chosen lambda keeps recovery; the trap loses it", {
  f <- sl_c2()
  sel <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                   scale_matched = FALSE, K = 4, seed = 11)
  J <- length(f$tables)
  rec <- vapply(sl_grid, function(l) {
    fit <- mgcca:::.mgcca_fit_penalized_dense(f$tables, rep(l, J), 2L)
    sl_rv(fit$Y, f$Z[fit$ids, , drop = FALSE])
  }, numeric(1))

  best     <- max(rec)
  loss_cv  <- best - rec[match(sel$lambda[1], sl_grid)]
  loss_trap<- best - rec[which.max(sl_grid)]

  # the shape of the thesis result (0.002 vs 0.666): CV nearly free, trap costly
  expect_lt(loss_cv, 0.05)
  expect_gt(loss_trap, 0.20)
  expect_gt(loss_trap, 4 * max(loss_cv, 1e-3))
})


## ---- the internal estimator the criterion is built on ------------------------

test_that("the dense penalized fit is the in-memory twin of method='penalized'", {
  skip_if_not(requireNamespace("BigDataStatMeth", quietly = TRUE))
  tabs <- readRDS(fixture_path("mgcca_subset.rds"))
  lam  <- rep(0.1, length(tabs))

  h5  <- tempfile(fileext = ".h5")
  ref <- mgcca(tabs, filename = h5, nfac = 2, method = "penalized",
               lambda = lam, scores = FALSE, collect = TRUE)
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)

  fit <- mgcca:::.mgcca_fit_penalized_dense(tabs, lam, 2L)
  # Same individuals. NOT necessarily in the same order: the C++ pipeline sorts
  # the union whenever the tables differ in size (src/mgcca_phases.hpp:166),
  # while the dense twin keeps the first-seen union order of the calibration
  # code (estimator.R:49). Row order is a labelling of the same fit -- everything
  # downstream aligns by ID, as here.
  expect_setequal(fit$ids, rownames(ref$Y))
  # same shared subspace (RV is rotation/sign invariant, which is all a subspace
  # is defined up to)
  expect_gt(sl_rv(fit$Y, ref$Y[fit$ids, , drop = FALSE]), 0.999)
})

test_that("block weights use the Gram (dual) form and equal the primal form", {
  # B_j = (X'X + lambda I)^-1 X'Y  ==  X'(XX' + lambda I)^-1 Y  (push-through).
  # estimator.R:461 forms the p x p primal version; the package must use the dual
  # one, which is the same number and the only one that survives p >> n.
  f   <- sl_sim(n = 40, p = c(15, 12), L = 2, snr = 3, specstr = 1, seed = 99)
  fit <- mgcca:::.mgcca_fit_penalized_dense(f$tables, c(2, 2), 2L)
  for (j in 1:2) {
    X  <- fit$Xlist[[j]]
    primal <- base::solve(base::crossprod(X) + 2 * base::diag(ncol(X)),
                          base::crossprod(X, fit$Y))
    dual   <- mgcca:::.mgcca_block_weights(fit, j)
    expect_equal(dual, primal, tolerance = 1e-8)
  }
})

test_that("the scale-matched block scale is tr(G_j)/rank(G_j)", {
  # estimator.R:96-124
  f <- sl_sim(n = 30, p = c(20, 10), L = 2, snr = 3, specstr = 1, seed = 7)
  s <- mgcca:::.mgcca_lambda_scale(f$tables)
  expect_length(s, 2L)
  for (j in 1:2) {
    ev <- base::eigen(base::tcrossprod(base::scale(f$tables[[j]])),
                      symmetric = TRUE, only.values = TRUE)$values
    ev[ev < 0] <- 0
    r <- sum(ev > 1e-8 * max(ev))
    expect_equal(unname(s[j]), sum(ev) / r, tolerance = 1e-8)
  }
})


## ---- the scale-matched (dimensionless gamma) parameterisation ---------------

test_that("scale_matched = TRUE selects gamma and reports lambda = gamma * s_j", {
  f <- sl_c2()
  sel <- sl_select(f$tables, nfac = 2,
                   gamma = 10^seq(-1, 1, by = 0.5), K = 4, seed = 5)
  expect_true(sel$scale_matched)
  expect_true(is.numeric(sel$gamma) && length(sel$gamma) == 1L)
  expect_true("gamma" %in% names(sel$grid))
  expect_equal(unname(sel$lambda), unname(sel$gamma * sel$s), tolerance = 1e-10)
  expect_length(sel$s, length(f$tables))
  expect_equal(sum(sel$grid$selected), 1L)
})


## ---- contracts and refusals --------------------------------------------------

test_that("mgcca() default is untouched: penalized still demands an explicit lambda", {
  tabs <- readRDS(fixture_path("mgcca_subset.rds"))
  expect_error(mgcca(tabs, filename = tempfile(fileext = ".h5"),
                     method = "penalized", collect = FALSE),
               "requires 'lambda'")
  # and mgcca() gained no lambda-selection argument
  expect_false(any(c("select_lambda", "lambda_grid") %in% names(formals(mgcca))))
})

test_that("bad input is refused, not guessed", {
  f <- sl_c2()
  expect_error(mgcca_select_lambda(f$tables[1], nfac = 2, lambda = c(1, 10),
                                   scale_matched = FALSE),
               "at least two")
  expect_error(mgcca_select_lambda(f$tables, nfac = 2, lambda = c(-1, 1),
                                   scale_matched = FALSE),
               "positive")
  expect_error(mgcca_select_lambda(f$tables, nfac = 0, lambda = c(1, 10),
                                   scale_matched = FALSE),
               "nfac")
  unnamed <- unname(lapply(f$tables, function(m) { rownames(m) <- NULL; m }))
  expect_error(mgcca_select_lambda(unnamed, nfac = 2, lambda = c(1, 10),
                                   scale_matched = FALSE),
               "rownames")
})

test_that("a boundary optimum is reported, never selected silently", {
  f <- sl_c2()
  # a grid that stops before the CV curve turns over: the argmax lands on the
  # lower edge, which cv_crossblock_gamma's contract calls a boundary optimum.
  expect_warning(
    sel <- mgcca_select_lambda(f$tables, nfac = 2, lambda = c(1e5, 1e6, 1e7),
                               scale_matched = FALSE, K = 4, seed = 11),
    class = "mgcca_lambda_boundary")
  expect_identical(sel$status, "boundary_optimum")
})

test_that("the 1-SE rule is available and never less regularised than argmax", {
  f <- sl_c2()
  a <- sl_select(f$tables, nfac = 2, lambda = sl_grid,
                 scale_matched = FALSE, K = 4, seed = 11)
  b <- sl_select(f$tables, nfac = 2, lambda = sl_grid, rule = "1se",
                 scale_matched = FALSE, K = 4, seed = 11)
  expect_identical(b$rule, "1se")
  expect_gte(b$lambda[1], a$lambda[1])
  expect_true(b$grid$in_1se[b$grid$selected])
  # same CV surface either way: only the pick differs
  expect_equal(a$grid$cv_mean, b$grid$cv_mean)
})
