# --- Bioconductor build time -------------------------------------------------
# This file is one of the four heaviest in the suite. Bioconductor's builders
# cap `R CMD check` at 10 minutes and this suite is the long pole, so the file
# is skipped THERE only (`IS_BIOC_BUILD_MACHINE`); it runs in full everywhere
# else, including on CRAN-style checks and in development.
testthat::skip_on_bioc()

# Downstream analysis + object methods on the shipped cardiovascular data.

make_X <- function() {
  data(cardiovascular, package = "mgcca", envir = environment())
  X <- list(methylation = as.matrix(get("X1")),
            clinical    = as.matrix(get("X2")),
            other       = as.matrix(get("X3")))
  lapply(X, function(m) { storage.mode(m) <- "double"; m })
}

test_that("mgcca() returns an mgcca object by default (collect = TRUE)", {
  X   <- make_X()
  fit <- mgcca(X, filename = tempfile(fileext = ".h5"), method = "solve",
               scores = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
  expect_s3_class(fit, "mgcca")
  expect_false(is.null(attr(fit, "desc")))
  expect_false(is.null(fit$weights))
  expect_false(is.null(fit$scaling))
  # collect = FALSE gives the descriptor
  d <- mgcca(X, filename = tempfile(fileext = ".h5"), method = "solve",
             collect = FALSE)
  expect_true(is.list(d) && !is.null(d$eig_values) && is.null(attr(d, "class")))
})

test_that("mgcca_associate() reports R2/p per component", {
  X   <- make_X()
  fit <- mgcca(X, filename = tempfile(fileext = ".h5"), method = "solve",
               scores = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
  g <- cut(X$clinical[, "LDL"][rownames(fit$Y)], c(-Inf, 130, Inf),
           c("normal", "high"))
  a <- mgcca_associate(fit, list(LDL = g, glucose = X$clinical[, "Gluc"][rownames(fit$Y)]))
  expect_true(all(c("phenotype", "component", "R2", "p", "n", "p.adj") %in% names(a)))
  expect_equal(nrow(a), 4L)                       # 2 phenotypes x 2 components
  expect_true(all(a$R2 >= 0 & a$R2 <= 1))
  expect_true(all(a$p[!is.na(a$p)] >= 0 & a$p[!is.na(a$p)] <= 1))
})

test_that("predict() reproduces training scores and projects new individuals", {
  X   <- make_X()
  fit <- mgcca(X, filename = tempfile(fileext = ".h5"), method = "solve",
               scores = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
  pr  <- predict(fit, newdata = list(clinical = X$clinical))
  ids <- rownames(X$clinical)
  expect_equal(unname(pr$clinical[ids, ]), unname(fit$scores$clinical[ids, ]),
               tolerance = 1e-8)
  pr2 <- predict(fit, newdata = list(clinical = X$clinical[1:5, , drop = FALSE]))
  expect_equal(dim(pr2$clinical), c(5L, 2L))
  expect_error(predict(fit, newdata = list(nope = X$clinical)), "used at fit time")
})

test_that("mgcca_permtest() returns observed eigenvalues and p-values", {
  X  <- make_X()
  pt <- mgcca_permtest(X, nperm = 9, method = "solve")
  expect_length(pt$eigenvalues, 2L)
  expect_length(pt$p.value, 2L)
  expect_true(all(pt$p.value >= 1 / (pt$nperm + 1) & pt$p.value <= 1))
  expect_equal(dim(pt$null), c(9L, 2L))
})

test_that("plot.mgcca / plotScores / plotBiplot return ggplots", {
  skip_if_not_installed("ggplot2")
  X   <- make_X()
  fit <- mgcca(X, filename = tempfile(fileext = ".h5"), method = "solve",
               scores = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plotScores(fit, table = "clinical"), "ggplot")
  expect_s3_class(plotBiplot(fit, table = "clinical", top = 4), "ggplot")
})
