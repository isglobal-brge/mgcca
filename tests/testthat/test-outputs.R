# --- Bioconductor build time -------------------------------------------------
# This file is one of the four heaviest in the suite. Bioconductor's builders
# cap `R CMD check` at 10 minutes and this suite is the long pole, so the file
# is skipped THERE only (`IS_BIOC_BUILD_MACHINE`); it runs in full everywhere
# else, including on CRAN-style checks and in development.
testthat::skip_on_bioc()

# =============================================================================
# Selectable outputs (wave 2, feature 5a).
#
# `outputs =` lets a caller ask for a SUBSET of the components an mgcca fit
# collects, so the HDF5 datasets behind the rest are never read. It must be
# purely SUBTRACTIVE:
#
#   * with the default, the object is bit-for-bit what 1.2.0 returned;
#   * with a subset, each requested component is identical() to its
#     full-collection counterpart -- selecting outputs changes what is READ,
#     never what is COMPUTED.
#
# The skipping itself is proved directly: a dataset that a subset does not ask
# for is DELETED from the file, and the subset still collects while the full
# collection fails on it.
# =============================================================================

out_fit <- function(fixture = "mgcca_subset.rds", method = "penalized",
                    lambda = rep(0.1, 3), scores = TRUE) {
  tabs <- readRDS(fixture_path(fixture))
  if (!is.null(lambda)) lambda <- rep(lambda[1], length(tabs))
  h5   <- tempfile(fileext = ".h5")
  desc <- mgcca(tabs, filename = h5, nfac = 2, method = method,
                lambda = lambda, scores = scores, collect = FALSE)
  list(desc = desc, h5 = h5)
}

OUT_ALL <- c("Y", "corsY", "scores", "pval", "weights", "scaling",
             "AVE", "eigen", "overlap")
OBJ_ALL <- c("Y", "corsY", "scores", "pval.cor", "weights", "scaling",
             "AVE", "eigen", "overlap")


## ---- the argument exists -----------------------------------------------------

test_that("mgcca() and mgcca_results() take an outputs argument", {
  expect_true("outputs" %in% names(formals(mgcca)))
  expect_true("outputs" %in% names(formals(mgcca_results)))
  expect_null(eval(formals(mgcca)$outputs))          # default = everything
  expect_null(eval(formals(mgcca_results)$outputs))
})


## ---- the default is bit-for-bit the old behaviour ---------------------------

test_that("outputs = NULL collects exactly what 1.2.0 collected", {
  f <- out_fit()
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)

  legacy <- mgcca_results(f$desc)                       # no outputs argument
  full   <- mgcca_results(f$desc, outputs = NULL)
  expect_identical(full, legacy)
  expect_identical(names(full), OBJ_ALL)
  expect_null(attr(full, "outputs"))                   # no new attribute either

  # asking for everything by name is the same object as asking for nothing
  every <- mgcca_results(f$desc, outputs = OUT_ALL)
  for (nm in OBJ_ALL) expect_identical(every[[nm]], legacy[[nm]])
})

test_that("mgcca(collect = TRUE) is unchanged when outputs is not given", {
  tabs <- readRDS(fixture_path("mgcca_subset.rds"))
  h5a  <- tempfile(fileext = ".h5"); h5b <- tempfile(fileext = ".h5")
  a <- mgcca(tabs, filename = h5a, nfac = 2, method = "penalized",
             lambda = rep(0.1, length(tabs)), scores = TRUE, collect = TRUE)
  b <- mgcca(tabs, filename = h5b, nfac = 2, method = "penalized",
             lambda = rep(0.1, length(tabs)), scores = TRUE, collect = TRUE,
             outputs = NULL)
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  attr(a, "desc") <- NULL; attr(b, "desc") <- NULL     # filenames differ
  expect_identical(a, b)
  expect_identical(class(a), "mgcca")
})


## ---- subsets are subtractive -------------------------------------------------

test_that("a subset returns the requested components, identical to the full ones", {
  f <- out_fit()
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  full <- mgcca_results(f$desc)

  sub <- mgcca_results(f$desc, outputs = c("Y", "weights", "eig"))
  expect_s3_class(sub, "mgcca")
  expect_identical(names(sub), OBJ_ALL)                # shape is preserved
  expect_identical(sub$Y,       full$Y)
  expect_identical(sub$weights, full$weights)
  expect_identical(sub$eigen,   full$eigen)
  for (nm in setdiff(OBJ_ALL, c("Y", "weights", "eigen")))
    expect_null(sub[[nm]])
  expect_identical(attr(sub, "outputs"), c("Y", "weights", "eigen"))
})

test_that("every single component can be asked for on its own and matches", {
  f <- out_fit()
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  full <- mgcca_results(f$desc)
  for (i in seq_along(OUT_ALL)) {
    one <- mgcca_results(f$desc, outputs = OUT_ALL[i])
    expect_identical(one[[OBJ_ALL[i]]], full[[OBJ_ALL[i]]])
    expect_true(all(vapply(setdiff(OBJ_ALL, OBJ_ALL[i]),
                           function(nm) is.null(one[[nm]]), logical(1))))
  }
})

test_that("mgcca(outputs = ...) subsets the collected object the same way", {
  tabs <- readRDS(fixture_path("mgcca_subset.rds"))
  h5   <- tempfile(fileext = ".h5")
  obj  <- mgcca(tabs, filename = h5, nfac = 2, method = "penalized",
                lambda = rep(0.1, length(tabs)), scores = TRUE, collect = TRUE,
                outputs = c("Y", "AVE"))
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  expect_false(is.null(obj$Y)); expect_false(is.null(obj$AVE))
  expect_null(obj$corsY); expect_null(obj$scores)
  # the descriptor is untouched by the choice of outputs
  d <- attr(obj, "desc")
  expect_false(is.null(d$eigen))
  expect_equal(as.character(d$datasets), names(tabs))
})


## ---- the reads really are skipped -------------------------------------------

test_that("a component that is not requested is not read from the file", {
  f <- out_fit()
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  ds <- f$desc$datasets[1]
  path <- paste0("FINAL_RESULTS/corsY/", ds)

  full_before <- mgcca_results(f$desc)                 # works while it is there
  expect_false(is.null(full_before$corsY))

  BigDataStatMeth::hdf5_close_all()
  BigDataStatMeth::hdf5_remove(f$h5, path)             # take corsY away

  # the full collection now fails on exactly that dataset ...
  expect_error(mgcca_results(f$desc))
  # ... while a subset that never asks for it still succeeds, and returns the
  # same Y as before.
  sub <- mgcca_results(f$desc, outputs = "Y")
  expect_identical(sub$Y, full_before$Y)
})

test_that("the read plan lists exactly the components asked for", {
  p_all <- mgcca:::.mgcca_output_plan(NULL, scores = TRUE, pval = TRUE)
  expect_identical(sort(names(p_all)[p_all]), sort(OUT_ALL))

  p <- mgcca:::.mgcca_output_plan(c("Y", "eig"), scores = TRUE, pval = TRUE)
  expect_identical(sort(names(p)[p]), sort(c("Y", "eigen")))

  # outputs never RESURRECTS something the fit did not compute: scores = FALSE
  # still means no scores, even if they are asked for by name.
  p0 <- mgcca:::.mgcca_output_plan(c("Y", "scores"), scores = FALSE, pval = TRUE)
  expect_false(p0[["scores"]])
  expect_true(p0[["Y"]])
})


## ---- refusals ----------------------------------------------------------------

test_that("an unknown output name is an error, not a silent drop", {
  f <- out_fit()
  on.exit(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE), add = TRUE)
  expect_error(mgcca_results(f$desc, outputs = c("Y", "Ynot")), "Ynot")
  expect_error(mgcca_results(f$desc, outputs = character(0)), "at least one")
  expect_error(mgcca_results(f$desc, outputs = 1:3), "character")
})
