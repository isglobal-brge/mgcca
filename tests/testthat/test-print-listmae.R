# print() on what getTables() actually returns. The class the method is
# registered for and the class getTables() assigns used to differ in the case of
# one letter, so the method never fired; the check below is on dispatch, not on
# the text it writes.

test_that("print() dispatches to print.ListMAE on a getTables() result", {
  skip_if_not_installed("MultiAssayExperiment")

  data(cardiovascular, package = "mgcca", envir = environment())
  sel <- rownames(get("X2"))[1:20]
  a1 <- t(as.matrix(get("X1")[rownames(get("X1")) %in% sel, 1:8, drop = FALSE]))
  a2 <- t(as.matrix(get("X2")[rownames(get("X2")) %in% sel, , drop = FALSE]))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
      MultiAssayExperiment::ExperimentList(
          list(methylation = a1, clinical = a2)))

  tabs <- getTables(mae)
  expect_s3_class(tabs, "ListMAE")

  # The method is registered for the class getTables() assigns, so it is the
  # one UseMethod("print") selects.
  expect_identical(
      utils::getS3method("print", "ListMAE"),
      utils::getFromNamespace("print.ListMAE", "mgcca"))

  out <- utils::capture.output(print(tabs))
  expect_true(any(grepl("Object of class ListMAE", out, fixed = TRUE)))
  expect_true(any(grepl("2 assays", out, fixed = TRUE)))
  expect_true(any(grepl("methylation", out, fixed = TRUE)))

  # Auto-printing takes the same route.
  expect_identical(out, utils::capture.output(tabs))
})
