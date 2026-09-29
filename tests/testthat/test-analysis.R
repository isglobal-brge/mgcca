# Downstream analysis + object methods on the shipped cardiovascular data.
#
# The fixture is a deliberately small slice of `cardiovascular`: 60 individuals
# and 10 / 6 / 4 features. The slice is chosen so that some individuals are
# ABSENT from some tables -- that is the structure these functions exist for,
# and a complete-case slice would not exercise it. Every assertion below is the
# same property the full-data version asserted; only the fixture is smaller, so
# that the file runs on the Bioconductor builders instead of being skipped.
#
# The fit itself is the expensive part, so the tests that only read a fit share
# one (`an_fit()`), rather than refitting the same model each time.

# Every temporary HDF5 file this file creates is removed when it is done.
h5_cleanup_on_exit()

make_X <- function() {
  data(cardiovascular, package = "mgcca", envir = environment())
  X1 <- as.matrix(get("X1"))
  X2 <- as.matrix(get("X2"))
  X3 <- as.matrix(get("X3"))
  u <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))
  # individuals that some table does not carry, kept on purpose
  gaps <- c(head(setdiff(u, rownames(X1)), 8L),
            head(setdiff(u, rownames(X2)), 6L),
            setdiff(u, rownames(X3)))
  ids <- unique(c(gaps, u))[seq_len(60L)]
  keep <- function(m, cols) {
    m <- m[rownames(m) %in% ids, cols, drop = FALSE]
    storage.mode(m) <- "double"
    m
  }
  list(methylation = keep(X1, seq_len(10L)),
       clinical    = keep(X2, colnames(X2)),
       other       = keep(X3, colnames(X3)))
}

# One fit, built on first use and reused by every test that only reads it.
.an_cache <- new.env(parent = emptyenv())

an_fit <- function() {
  if (!is.null(.an_cache$fit)) return(.an_cache$fit)
  X <- make_X()
  fit <- mgcca(X, filename = h5_tmp(), method = "solve",
               scores = TRUE)
  try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
  .an_cache$X <- X
  .an_cache$fit <- fit
  fit
}

an_X <- function() { an_fit(); .an_cache$X }

test_that("mgcca() returns an mgcca object by default (collect = TRUE)", {
  fit <- an_fit()
  expect_s3_class(fit, "mgcca")
  expect_false(is.null(attr(fit, "desc")))
  expect_false(is.null(fit$weights))
  expect_false(is.null(fit$scaling))
  # collect = FALSE gives the descriptor
  d <- mgcca(an_X(), filename = h5_tmp(), method = "solve",
             collect = FALSE)
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  expect_true(is.list(d) && !is.null(d$eig_values) && is.null(attr(d, "class")))
})

test_that("mgcca_associate() reports R2/p per component", {
  fit <- an_fit()
  X <- an_X()
  g <- cut(X$clinical[, "LDL"][rownames(fit$Y)], c(-Inf, 130, Inf),
           c("normal", "high"))
  a <- mgcca_associate(fit, list(LDL = g, glucose = X$clinical[, "Gluc"][rownames(fit$Y)]))
  expect_true(all(c("phenotype", "component", "R2", "p", "n", "p.adj") %in% names(a)))
  expect_equal(nrow(a), 4L)                       # 2 phenotypes x 2 components
  expect_true(all(a$R2 >= 0 & a$R2 <= 1))
  expect_true(all(a$p[!is.na(a$p)] >= 0 & a$p[!is.na(a$p)] <= 1))
})

test_that("predict() reproduces training scores and projects new individuals", {
  fit <- an_fit()
  X <- an_X()
  pr  <- predict(fit, newdata = list(clinical = X$clinical))
  ids <- rownames(X$clinical)
  expect_equal(unname(pr$clinical[ids, ]), unname(fit$scores$clinical[ids, ]),
               tolerance = 1e-8)
  pr2 <- predict(fit, newdata = list(clinical = X$clinical[1:5, , drop = FALSE]))
  expect_equal(dim(pr2$clinical), c(5L, 2L))
  expect_error(predict(fit, newdata = list(nope = X$clinical)), "used at fit time")
})

test_that("mgcca_permtest() returns observed eigenvalues and p-values", {
  # nperm is the number of REFITS, so it is the one knob that has to stay small;
  # 4 is enough for the shape and the 1/(nperm+1) floor the assertions are about.
  pt <- mgcca_permtest(an_X(), filename = h5_tmp(), nperm = 4, method = "solve")
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  expect_length(pt$eigenvalues, 2L)
  expect_length(pt$p.value, 2L)
  expect_true(all(pt$p.value >= 1 / (pt$nperm + 1) & pt$p.value <= 1))
  expect_equal(dim(pt$null), c(4L, 2L))
})

test_that("plot.mgcca / plotScores / plotBiplot return ggplots", {
  skip_if_not_installed("ggplot2")
  fit <- an_fit()
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plotScores(fit, table = "clinical"), "ggplot")
  expect_s3_class(plotBiplot(fit, table = "clinical", top = 4), "ggplot")
})
