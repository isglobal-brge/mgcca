# mgcca_results() / mgcca(collect=TRUE) and the post-analysis helpers, on the
# shipped cardiovascular data (three tables with missing individuals). This runs
# anywhere -- it needs no external reference fixtures.

make_X <- function() {
  data(cardiovascular, package = "mgcca", envir = environment())
  X <- list(methylation = as.matrix(get("X1")),
            clinical    = as.matrix(get("X2")),
            other       = as.matrix(get("X3")))
  lapply(X, function(m) { storage.mode(m) <- "double"; m })
}

test_that("mgcca(collect=TRUE) yields a usable 'mgcca' object", {
  X   <- make_X()
  h5  <- tempfile(fileext = ".h5")
  fit <- mgcca(X, filename = h5, nfac = 2, method = "solve",
                   scores = TRUE, collect = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  expect_s3_class(fit, "mgcca")
  m <- length(Reduce(union, lapply(X, rownames)))
  expect_equal(dim(fit$Y), c(m, 2L))
  expect_false(is.null(rownames(fit$Y)))
  expect_true(all(is.finite(fit$Y)))

  expect_named(fit$corsY, names(X))
  expect_named(fit$scores, names(X))
  expect_named(fit$pval.cor, names(X))
  for (nm in names(X)) {
    expect_equal(nrow(fit$corsY[[nm]]), ncol(X[[nm]]))
    expect_identical(rownames(fit$corsY[[nm]]), colnames(X[[nm]]))
  }
  expect_length(fit$AVE$AVE_inner_model, 2L)
})

test_that("mgcca_results() reads the same object from the descriptor and the path", {
  X   <- make_X()
  h5  <- tempfile(fileext = ".h5")
  desc <- mgcca(X, filename = h5, nfac = 2, method = "solve", scores = TRUE,
                collect = FALSE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  a <- mgcca_results(desc)
  b <- mgcca_results(h5, datasets = desc$datasets)
  expect_s3_class(a, "mgcca")
  expect_equal(a$Y, b$Y)
  expect_equal(a$corsY, b$corsY)
})

test_that("plot/summary helpers accept the collected object without error", {
  X   <- make_X()
  h5  <- tempfile(fileext = ".h5")
  fit <- mgcca(X, filename = h5, nfac = 2, method = "solve",
                   scores = TRUE, collect = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plotInds(fit))
  expect_no_error(plotVars(fit, nlab = 3))

  sig <- getSignif(fit, pval.cut = 1e-3)
  expect_true(all(c("variable", "table") %in% names(sig)))
  tv <- topVars(fit, axis = 1, topN = 3)
  expect_named(tv, names(X))
})

test_that("modern ggplot result plots build for an mgcca object", {
  skip_if_not_installed("ggplot2")
  X   <- make_X()
  h5  <- tempfile(fileext = ".h5")
  fit <- mgcca(X, filename = h5, nfac = 2, method = "solve",
                   scores = TRUE, collect = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  expect_s3_class(plotIndividuals(fit), "ggplot")
  expect_s3_class(plotVariables(fit, table = "clinical", top = 4), "ggplot")
  expect_s3_class(plotVariables(fit, table = 1), "ggplot")       # index also works
  expect_s3_class(plotAVE(fit), "ggplot")
  expect_s3_class(plotScree(fit), "ggplot")
  expect_s3_class(plotLoadings(fit, table = "clinical", comp = 1, top = 4), "ggplot")

  # a grouping factor with names is matched to the individuals
  g <- factor(stats::setNames(rep(c("a", "b"), length.out = nrow(fit$Y)),
                              rownames(fit$Y)))
  expect_s3_class(plotIndividuals(fit, group = g), "ggplot")

  # they actually render without error
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_no_error(print(plotIndividuals(fit, group = g)))
  expect_no_error(print(plotVariables(fit, table = "clinical")))

  # input validation
  expect_error(plotIndividuals(list()), "mgcca")
})

test_that("print/summary methods work on an mgcca object", {
  X   <- make_X()
  h5  <- tempfile(fileext = ".h5")
  fit <- mgcca(X, filename = h5, nfac = 2, method = "solve",
               scores = TRUE, collect = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  expect_output(print(fit), "Generalized Canonical Correlation")
  expect_output(print(fit), "tables")
  expect_output(summary(fit, top = 2), "AVE per table")
  expect_output(summary(fit, top = 2), "Top 2 variables")
  expect_invisible(print(fit))
})
