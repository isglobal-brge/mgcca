# Provenance manifest persistence (native HDF5 attributes) and mgcca_load().
# mgcca() now writes a self-describing manifest; mgcca_load() rebuilds the full
# object from the file alone. Also covers the no-manifest fallback path.

test_that("mgcca_load rebuilds the object and provenance from the file alone", {
  skip_if_not(have_reference("mgcca_subset.rds", "penalized"),
              "reference fixtures not shipped with the built package")

  o <- run_fixture("mgcca_subset.rds", "penalized", lambda = c(0.75, 1))
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
  ds <- o$res$datasets

  fit  <- mgcca_load(o$h5)
  desc <- attr(fit, "desc")

  # object shape
  expect_s3_class(fit, "mgcca")
  expect_equal(dim(fit$Y), c(80L, 2L))
  expect_named(fit$corsY, ds)

  # numeric parity with the collector on the same file (sign-tolerant for Y)
  ref <- mgcca_results(o$res)
  expect_true(eqs(fit$Y, ref$Y))
  for (d in ds) expect_true(eqs(fit$corsY[[d]], ref$corsY[[d]]), info = d)

  # provenance recovered from the manifest
  expect_equal(desc$method, "penalized")
  expect_equal(desc$nfac, 2L)
  expect_equal(desc$m, 80L)
  expect_equal(as.character(desc$datasets), ds)
  expect_equal(length(desc$eig_values), 2L)
  expect_equal(as.numeric(desc$eig_values), as.numeric(o$res$eig_values),
               tolerance = 1e-8)

  # per-table facts: lambda + route_dual, named by dataset
  expect_named(desc$lambda, ds)
  expect_equal(unname(desc$lambda), c(0.75, 1), tolerance = 1e-8)
  expect_named(desc$route_dual, ds)
  expect_type(desc$route_dual, "logical")
  # mgcca_subset is p < n -> covariance route -> not dual
  expect_equal(unname(desc$route_dual), c(FALSE, FALSE))

  # version/date stamped
  expect_true(nchar(desc$mgcca_version) > 0)
  expect_match(desc$mgcca_date, "^[0-9]{4}-[0-9]{2}-[0-9]{2}")
})

test_that("mgcca_load falls back to corsY discovery when no manifest exists", {
  skip_if_not(have_reference("mgcca_subset.rds", "penalized"),
              "reference fixtures not shipped with the built package")

  o <- run_fixture("mgcca_subset.rds", "penalized", lambda = c(0.75, 1))
  ref <- mgcca_results(o$res)
  BigDataStatMeth::hdf5_close_all()
  ds <- o$res$datasets

  # Build an equivalent results file WITHOUT any provenance attributes, by
  # re-writing only the numeric blocks that mgcca_results() consumes.
  h5 <- tempfile(fileext = ".h5")
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
  cm <- function(path, m)
    BigDataStatMeth::hdf5_create_matrix(h5, path, data = as.matrix(m),
                                        overwrite = TRUE)
  cm("FINAL_RESULTS/Y", ref$Y)
  for (d in ds) {
    cm(paste0("FINAL_RESULTS/corsY/", d), ref$corsY[[d]])
    cm(paste0("FINAL_RESULTS/pval/", d),  ref$pval.cor[[d]])
  }
  cm("FINAL_RESULTS/AVE/AVE_X",     ref$AVE$AVE_X)
  cm("FINAL_RESULTS/AVE/AVE_outer", matrix(ref$AVE$AVE_outer_model, ncol = 1))
  cm("FINAL_RESULTS/AVE/AVE_inner", matrix(ref$AVE$AVE_inner_model, ncol = 1))
  BigDataStatMeth::hdf5_close_all()

  fit  <- mgcca_load(h5)
  desc <- attr(fit, "desc")

  expect_s3_class(fit, "mgcca")
  expect_setequal(as.character(desc$datasets), ds)   # discovered from corsY
  expect_true(eqs(fit$Y, ref$Y))
  # provenance is unknown without a manifest
  expect_true(is.na(desc$method))
  expect_true(is.na(desc$route))
})
