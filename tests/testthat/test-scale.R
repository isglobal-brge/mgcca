# Real-scale smoke test of the dual route: full TCGA-ACC methylation (~485k
# features x 80) + RNASeq. Heavy and depends on the 243 MB imputed MAE, which is
# not shipped -- skipped unless it is present locally. Confirms route="auto"
# picks the Gram route for both p>>n tables and the pipeline completes with a
# finite, well-shaped Y / scores (no p x p ever formed).

test_that("dual route scales to full methylation (p >> n)", {
  skip_on_cran()
  mae <- testthat::test_path("..", "data", "bd", "multiassayexperiment.imputed.rda")
  skip_if_not(file.exists(mae), "imputed MAE not present")

  suppressMessages({
    library(MultiAssayExperiment); library(SummarizedExperiment)
  })
  e <- new.env(); load(mae, envir = e)
  obj <- get(ls(e)[1], envir = e)
  tabs <- list()
  for (j in seq_along(experiments(obj))) {
    m <- t(as.matrix(assay(obj[[j]]))); storage.mode(m) <- "double"
    rownames(m) <- substr(rownames(m), 1, 12)
    keep <- apply(m, 2, function(v) all(is.finite(v)) && sd(v) > 1e-8)
    tabs[[names(experiments(obj))[j]]] <- m[, keep, drop = FALSE]
  }
  rm(obj, e); invisible(gc())

  h5  <- tempfile(fileext = ".h5")
  res <- mgcca(tabs, filename = h5, nfac = 2, method = "geninv",
               scores = TRUE, collect = FALSE)
  on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)

  expect_true(all(as.logical(res$route_dual)))     # both tables p >> n -> dual
  Y <- as.matrix(BigDataStatMeth::hdf5_matrix(h5, "FINAL_RESULTS/Y"))
  expect_equal(dim(Y), c(res$m, 2L))
  expect_true(all(is.finite(Y)))
  ids <- Reduce(union, lapply(tabs, rownames))
  for (ds in res$datasets) {
    S <- as.matrix(BigDataStatMeth::hdf5_matrix(h5, paste0("FINAL_RESULTS/scores/", ds)))
    expect_equal(nrow(S), res$m)
    # every row of the union is there, but the individuals a table does not hold
    # carry NaN in it (they have no score in that table), not a number
    absent <- !stats::complete.cases(S)
    expect_equal(sum(absent), length(setdiff(ids, rownames(tabs[[ds]]))))
    expect_true(all(is.nan(S[absent, , drop = FALSE])))
    expect_true(all(is.finite(S[!absent, , drop = FALSE])))
  }
})
