# Inversion methods (solve / geninv) vs the oracle, plus the route dispatch.
# For a well-conditioned p<n fixture with distinct eigenvalues we validate the
# full output; for the p>n fixture the top eigenvalue is degenerate (multiplicity),
# so Y is only defined up to a rotation and we validate the deterministic
# quantities (route-invariant Mi, Msum, MKsum05, eigenvalues) instead.

# Compare the intermediates + (optionally) the sign-ambiguous outputs against the
# oracle. `route`-aware: Mgram/XKX are covariance-route artefacts, absent under
# the dual (Gram) route, so they are only checked for cov-routed tables.
check_case <- function(fixture, method, lambda = NULL, route = "auto") {
  o   <- run_fixture(fixture, method, lambda = lambda, route = route)
  REF <- reference_dir(fixture, method)
  rd  <- function(n) readRDS(file.path(REF, paste0(n, ".rds")))
  read <- o$read; res <- o$res
  ev  <- as.numeric(res$eig_values)

  for (i in seq_along(res$datasets)) {
    ds <- res$datasets[i]
    if (!isTRUE(res$route_dual[[i]])) {
      expect_true(eqm(read(paste0("MGCCA_TMP/M/", ds)),   rd(sprintf("t%d_Mgram", i))), info = ds)
      expect_true(eqm(read(paste0("MGCCA_TMP/XKX/", ds)), rd(sprintf("t%d_xkx", i))),   info = ds)
    }
    expect_true(eqm(read(paste0("MGCCA_TMP/Mi/", ds)), rd(sprintf("t%d_Mi", i))), info = ds)
  }
  expect_true(eqm(read("MGCCA_TMP/Msum"),    rd("M")))
  expect_true(eqm(read("MGCCA_TMP/MKsum05"), rd("MKsum05")))
  expect_true(eqm(sort(ev), sort(as.numeric(rd("eig_values")))))

  if (abs(ev[1] - ev[2]) >= 1e-8) {              # non-degenerate: full check
    expect_true(eqs(read("FINAL_RESULTS/Y"), rd("Y")))
    for (i in seq_along(res$datasets)) {
      ds <- res$datasets[i]
      expect_true(eqs(read(paste0("FINAL_RESULTS/corsY/", ds)),  rd(sprintf("t%d_corsY", i))),  info = ds)  # sign-tolerant: corsY inherits Y's sign
      expect_true(eqs(read(paste0("FINAL_RESULTS/scores/", ds)), rd(sprintf("t%d_scores", i))), info = ds)
    }
  }
  BigDataStatMeth::hdf5_close_all()
  invisible(res)
}

test_that("solve matches the oracle (well-conditioned p<n)", {
  skip_if_not(have_reference("mgcca_subset_solve.rds", "solve"), "reference not shipped")
  res <- check_case("mgcca_subset_solve.rds", "solve")
  expect_false(any(as.logical(res$route_dual)))   # p<n -> covariance route
})

test_that("geninv matches the oracle (well-conditioned p<n)", {
  skip_if_not(have_reference("mgcca_subset_solve.rds", "geninv"), "reference not shipped")
  check_case("mgcca_subset_solve.rds", "geninv")
})

test_that("geninv on p>n: auto picks the dual route, deterministics match", {
  skip_if_not(have_reference("mgcca_subset_p100.rds", "geninv"), "reference not shipped")
  res <- check_case("mgcca_subset_p100.rds", "geninv", route = "auto")
  expect_true(all(as.logical(res$route_dual)))    # p>n -> dual route for every table
})

test_that("geninv on p>n via forced covariance route matches the oracle", {
  skip_if_not(have_reference("mgcca_subset_p100.rds", "geninv"), "reference not shipped")
  res <- check_case("mgcca_subset_p100.rds", "geninv", route = "cov")
  expect_false(any(as.logical(res$route_dual)))   # forced covariance
})
