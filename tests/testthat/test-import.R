# mgcca_import_hdf5: raw round-trip of values, orientation and dimnames.

test_that("mgcca_import_hdf5 round-trips values and dimnames", {
  skip_if_not(file.exists(fixture_path("mgcca_subset.rds")),
              "fixtures not shipped with the built package")

  fix  <- readRDS(fixture_path("mgcca_subset.rds"))
  h5   <- tempfile(fileext = ".h5")
  desc <- mgcca_import_hdf5(fix, filename = h5, group = "MGCCA_IN",
                            overwriteFile = TRUE, overwriteDataset = TRUE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  expect_setequal(desc$datasets, names(fix))
  for (nm in names(fix)) {
    hm <- BigDataStatMeth::hdf5_matrix(h5, paste0("MGCCA_IN/", nm))
    dn <- dimnames(hm)                                # /1 = rownames, /2 = colnames
    expect_true(eqm(as.matrix(hm), fix[[nm]], tol = 1e-10),
                info = paste("values", nm))
    expect_identical(as.character(dn[[1]]), rownames(fix[[nm]]))
    expect_identical(as.character(dn[[2]]), colnames(fix[[nm]]))
  }
})
