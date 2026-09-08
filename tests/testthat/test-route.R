# The covariance (X'X) and dual (Gram XX') routes are two algebraically
# equivalent ways to build M_j; on a well-conditioned fixture they must produce
# the same Mi, eigenvalues and shared components Y.

test_that("cov and dual routes agree on Mi, eigenvalues and Y", {
  skip_if_not(file.exists(fixture_path("mgcca_subset_solve.rds")),
              "fixtures not shipped with the built package")

  cov  <- run_fixture("mgcca_subset_solve.rds", "geninv", route = "cov")
  dual <- run_fixture("mgcca_subset_solve.rds", "geninv", route = "dual")
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  expect_false(any(as.logical(cov$res$route_dual)))
  expect_true(all(as.logical(dual$res$route_dual)))
  expect_equal(sort(as.numeric(cov$res$eig_values)),
               sort(as.numeric(dual$res$eig_values)), tolerance = 1e-5)
  for (ds in cov$res$datasets)
    expect_true(eqm(cov$read(paste0("MGCCA_TMP/Mi/", ds)),
                    dual$read(paste0("MGCCA_TMP/Mi/", ds)), tol = 1e-5), info = ds)
  expect_true(eqs(cov$read("FINAL_RESULTS/Y"), dual$read("FINAL_RESULTS/Y")))
})

test_that("route must be one of auto/cov/dual", {
  skip_if_not(file.exists(fixture_path("mgcca_subset_solve.rds")),
              "fixtures not shipped with the built package")
  tabs <- readRDS(fixture_path("mgcca_subset_solve.rds"))
  expect_error(
    mgcca(tabs, filename = tempfile(fileext = ".h5"),
              method = "geninv", route = "bogus"))
})
