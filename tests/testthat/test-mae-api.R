# The MultiAssayExperiment-facing half of the API -- getTables(), mgccaImpute(),
# matrixChr2Num(), genCoords() and corTestP() -- on a tiny object built here
# from scratch, plus the regressions for the argument-passing fixes. Everything
# is a known-value assertion: no fixture files, no reference dumps.

# A small MultiAssayExperiment: two assays, features in rows, individuals in
# columns, overlapping on s03..s05 only.
mae_fixture <- function(na = FALSE) {
  a1 <- matrix(as.double(1:40), nrow = 8, ncol = 5,
               dimnames = list(paste0("cg", sprintf("%02d", 1:8)),
                               paste0("s", sprintf("%02d", 1:5))))
  a2 <- matrix(as.double(101:118), nrow = 3, ncol = 6,
               dimnames = list(c("age", "bmi", "sbp"),
                               paste0("s", sprintf("%02d", 3:8))))
  if (na) {
    a1[2, 2] <- NA
    a1[5, 4] <- NA
  }
  MultiAssayExperiment::MultiAssayExperiment(
    MultiAssayExperiment::ExperimentList(
      list(methylation = a1, clinical = a2)))
}

test_that("getTables() transposes every assay and labels the list", {
  skip_if_not_installed("MultiAssayExperiment")
  mae  <- mae_fixture()
  tabs <- getTables(mae)

  expect_s3_class(tabs, "ListMAE")
  expect_named(tabs, c("methylation", "clinical"))
  expect_equal(dim(tabs$methylation), c(5L, 8L))   # individuals x variables
  expect_equal(dim(tabs$clinical),    c(6L, 3L))
  expect_identical(rownames(tabs$methylation),
                   paste0("s", sprintf("%02d", 1:5)))
  expect_identical(colnames(tabs$methylation),
                   paste0("cg", sprintf("%02d", 1:8)))
  expect_identical(colnames(tabs$clinical), c("age", "bmi", "sbp"))
  # values, not just shapes: column-major fill means cg01/s02 is 9
  expect_equal(tabs$methylation["s02", "cg01"], 9)
  expect_equal(tabs$clinical["s03", "age"], 101)
  expect_true(all(vapply(tabs, is.matrix, logical(1))))

  expect_error(getTables(list(1, 2)), "MultiAssayExperiment")
})

test_that("matrixChr2Num() converts one table and keeps the dimnames", {
  skip_if_not_installed("MultiAssayExperiment")
  mae <- mae_fixture()
  tabs <- getTables(mae)
  storage.mode(tabs$methylation) <- "character"
  expect_true(is.character(tabs$methylation))

  out <- matrixChr2Num(tabs, 1)

  expect_s3_class(out, "ListMAE")
  expect_true(is.numeric(out$methylation))
  expect_equal(dim(out$methylation), c(5L, 8L))
  expect_identical(dimnames(out$methylation), dimnames(tabs$methylation))
  expect_equal(unname(out$methylation), unname(t(matrix(as.double(1:40), 8, 5))))
  # the other table is untouched
  expect_identical(out$clinical, tabs$clinical)

  expect_error(matrixChr2Num(unclass(tabs), 1), "ListMAE")
  expect_error(matrixChr2Num(tabs), "Index")
})

test_that("mgccaImpute() fills the missing entries and leaves the rest alone", {
  skip_if_not_installed("MultiAssayExperiment")
  mae <- mae_fixture(na = TRUE)
  expect_true(anyNA(mae[[1]]))

  imp <- suppressMessages(mgccaImpute(mae, method = "knn", k = 3))

  expect_s4_class(imp, "MultiAssayExperiment")
  expect_false(anyNA(imp[[1]]))
  expect_true(all(is.finite(imp[[1]])))
  expect_equal(sum(!is.na(imp[[1]])), 40L)
  expect_identical(dimnames(imp[[1]]), dimnames(mae[[1]]))
  # the observed entries come back untouched
  obs <- !is.na(mae[[1]])
  expect_equal(imp[[1]][obs], mae[[1]][obs])
  expect_equal(imp[[2]], mae[[2]])
})

test_that("mgccaImpute() hands k, rowmax and colmax to impute.knn by name", {
  skip_if_not_installed("MultiAssayExperiment")
  mae <- mae_fixture(na = TRUE)

  # Same call, made directly: identical output means k/rowmax/colmax landed on
  # the formals of the same name (they used to be passed positionally, so
  # rowmax arrived as the neighbour count k).
  expected <- suppressMessages(
    impute::impute.knn(mae[[1]], k = 3, rowmax = 0.5, colmax = 0.8)$data)
  got <- suppressMessages(
    mgccaImpute(mae, method = "knn", k = 3, rowmax = 0.5, colmax = 0.8))
  expect_equal(got[[1]], expected)

  # k really is the neighbour count: a different k gives a different answer.
  other <- suppressMessages(
    mgccaImpute(mae, method = "knn", k = 5, rowmax = 0.5, colmax = 0.8))
  expect_false(isTRUE(all.equal(got[[1]], other[[1]])))

  # colmax really is colmax: impute.knn refuses a column that is too empty.
  bad <- mae
  a <- bad[[1]]
  a[, 2] <- NA_real_
  bad[[1]] <- a
  expect_error(suppressMessages(
    mgccaImpute(bad, method = "knn", k = 3, colmax = 0.5)),
    "a column has more than 50 % missing values", fixed = TRUE)
  # rowmax really is rowmax: a row too empty for kNN falls back to the column
  # mean, which is a different number from the kNN estimate.
  lax <- suppressWarnings(suppressMessages(
    mgccaImpute(mae, method = "knn", k = 3, rowmax = 0.001, colmax = 0.8)))
  expect_equal(lax[[1]], suppressWarnings(suppressMessages(
    impute::impute.knn(mae[[1]], k = 3, rowmax = 0.001, colmax = 0.8)$data)))
  expect_false(isTRUE(all.equal(got[[1]], lax[[1]])))
})

test_that("mgccaImpute() only accepts the method it implements", {
  skip_if_not_installed("MultiAssayExperiment")
  mae <- mae_fixture(na = TRUE)

  expect_error(mgccaImpute(mae, method = "hmisc"), "knn")
  expect_error(mgccaImpute(mae, method = "mean"), "knn")
  # partial matching still works, and the default is "knn"
  expect_false(anyNA(
    suppressMessages(mgccaImpute(mae, method = "kn", k = 3))[[1]]))
  expect_false(anyNA(suppressMessages(mgccaImpute(mae, k = 3))[[1]]))
  expect_error(mgccaImpute(list(1, 2), method = "knn"), "MultiAssayExperiment")
})

# genCoords() ------------------------------------------------------------

# A non-ranged assay (methylation-like) whose rowData carries the position and
# the chromosome, plus a ranged one. cg03 sits inside the region but on chr2, so
# it must be dropped; cg05 is on chr1 but outside the region.
gen_fixture <- function(pos.col = "position") {
  m1 <- matrix(as.double(1:20), nrow = 5, ncol = 4,
               dimnames = list(paste0("cg0", 1:5), paste0("s0", 1:4)))
  rd <- data.frame(pos = c(1200, 2200, 3200, 4200, 9000),
                   chr = c("chr1", "chr1", "chr2", "chr1", "chr1"),
                   row.names = rownames(m1))
  names(rd)[1] <- pos.col
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(beta = m1), rowData = rd)

  m2 <- matrix(as.double(1:12), nrow = 3, ncol = 4,
               dimnames = list(paste0("g0", 1:3), paste0("s0", 1:4)))
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1500, 2500, 9500), width = 50))
  names(gr) <- rownames(m2)
  rse2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = m2), rowRanges = gr)

  MultiAssayExperiment::MultiAssayExperiment(
    MultiAssayExperiment::ExperimentList(
      list(methylation = se1, cells = rse2)))
}

test_that("genCoords() honours col.name and the chromosome", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  mae <- gen_fixture(pos.col = "position")
  region <- GenomicRanges::GRanges(
    seqnames = "chr1", ranges = IRanges::IRanges(start = 1000, end = 5000))

  sub <- genCoords(mae, region, meth.index = 1, col.name = "position")

  # cg03 is inside the interval but on chr2; cg05 is past the end
  expect_identical(rownames(sub[[1]]), c("cg01", "cg02", "cg04"))
  # the ranged assay keeps the two features contained in the region
  expect_identical(rownames(sub[[2]]), c("g01", "g02"))

  # naming the chromosome column explicitly gives the same answer
  sub2 <- genCoords(mae, region, meth.index = 1, col.name = "position",
                    seqname.col = "chr")
  expect_identical(rownames(sub2[[1]]), rownames(sub[[1]]))

  # a region on the other chromosome picks up cg03 only
  region2 <- GenomicRanges::GRanges(
    seqnames = "chr2", ranges = IRanges::IRanges(start = 1000, end = 5000))
  # (the ranged assay has nothing on chr2, which GenomicRanges warns about)
  sub3 <- suppressWarnings(
    genCoords(mae, region2, meth.index = 1, col.name = "position"))
  expect_identical(rownames(sub3[[1]]), "cg03")
})

test_that("genCoords() validates its arguments", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  mae <- gen_fixture(pos.col = "position")
  region <- GenomicRanges::GRanges(
    seqnames = "chr1", ranges = IRanges::IRanges(start = 1000, end = 5000))

  expect_error(genCoords(list(1), region, 1, "position"),
               "MultiAssayExperiment")
  expect_error(genCoords(mae, "chr1:1000-5000", 1, "position"), "GRanges")

  # start()/end() would recycle silently over several ranges
  two <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1000, 6000), end = c(5000, 9500)))
  expect_error(genCoords(mae, two, 1, "position"), "exactly one range")

  # an unknown position column is reported by name, not as a NULL comparison
  expect_error(genCoords(mae, region, 1, "Genomic_Coordinate"),
               "Genomic_Coordinate")
  expect_error(genCoords(mae, region, 1), "col.name")
  expect_error(genCoords(mae, region, 1, "position", seqname.col = "nope"),
               "nope")

  # without any chromosome column the position filter stands alone, with a
  # warning rather than a silent wrong answer
  mae.nochr <- gen_fixture(pos.col = "position")
  rd <- SummarizedExperiment::rowData(mae.nochr[[1]])
  SummarizedExperiment::rowData(mae.nochr[[1]]) <- rd[, "position", drop = FALSE]
  expect_warning(
    sub <- genCoords(mae.nochr, region, meth.index = 1, col.name = "position"),
    "chromosome")
  expect_identical(rownames(sub[[1]]), c("cg01", "cg02", "cg03", "cg04"))
})

# corTestP() ------------------------------------------------------------

test_that("corTestP() reproduces the p-value of stats::cor.test()", {
  x <- c(2.1, 3.4, 1.8, 5.6, 4.2, 6.9, 3.3, 7.1, 5.5, 2.8)
  y <- c(1.4, 2.9, 2.2, 4.9, 4.4, 6.1, 2.6, 6.8, 5.2, 3.1)
  r <- stats::cor(x, y)
  n <- length(x)

  expect_equal(corTestP(r, n), stats::cor.test(x, y)$p.value)

  # vectorised in r, and the obvious anchors
  ps <- corTestP(c(0.1, 0.5, 0.9), 50)
  expect_length(ps, 3L)
  expect_true(all(diff(ps) < 0))
  expect_equal(corTestP(0, 30), 1)
  # and the closed form it documents, on a coefficient with no data behind it
  expect_equal(corTestP(0.35, 100),
               2 * stats::pt(0.35 * sqrt(98) / sqrt(1 - 0.35^2), 98,
                             lower.tail = FALSE))
})

# getSignif() ------------------------------------------------------------

# A hand-built 'mgcca' object: getSignif() reads x$pval.cor and nothing else,
# so a fit is not needed to pin down its table selection.
signif_fixture <- function() {
  mk <- function(v, nms) matrix(v, nrow = length(nms), ncol = 2,
                                dimnames = list(nms, c("Comp1", "Comp2")))
  structure(list(pval.cor = list(
    methylation = mk(c(0.001, 0.4, 0.2, 0.9), c("cg01", "cg02")),
    clinical    = mk(c(0.02, 0.5, 0.6, 0.001), c("age", "bmi")),
    other       = mk(c(0.7, 0.8, 0.9, 0.95), c("v1", "v2")))),
    class = "mgcca")
}

test_that("getSignif() accepts a vector of tables", {
  fit <- signif_fixture()

  # the default NA means every table
  all.sig <- getSignif(fit, pval.cut = 0.05)
  expect_identical(all.sig$variable, c("cg01", "age", "bmi"))
  expect_identical(all.sig$table,
                   c("methylation", "clinical", "clinical"))

  # a vector of positions used to abort with "the condition has length > 1"
  two <- getSignif(fit, df = c(1, 2), pval.cut = 0.05)
  expect_identical(two$variable, c("cg01", "age", "bmi"))
  expect_identical(unique(two$table), c("methylation", "clinical"))

  # a vector of names, and a single name
  by.name <- getSignif(fit, df = c("clinical", "other"), pval.cut = 0.05)
  expect_identical(by.name$variable, c("age", "bmi"))
  expect_identical(getSignif(fit, df = "clinical", pval.cut = 0.05)$variable,
                   c("age", "bmi"))
  expect_equal(nrow(getSignif(fit, df = 3, pval.cut = 0.05)), 0L)

  expect_error(getSignif(fit, df = c(1, 99)), "not a valid name")
  expect_error(getSignif(list(), df = 1), "mgcca object expected")
})

test_that("getSignif() points at the call that collects p-values", {
  fit <- signif_fixture()
  fit$pval.cor <- NULL
  expect_error(getSignif(fit), "mgcca_results")
  expect_error(getSignif(fit), "outputs")
})

# mgcca_permtest() -------------------------------------------------------

test_that("mgcca_permtest() removes the HDF5 file of every permutation", {
  n <- 14L
  ids <- paste0("i", sprintf("%02d", seq_len(n)))
  mk <- function(off, p) {
    m <- outer(seq_len(n), seq_len(p), function(i, j) sin(i * j + off) + i / n)
    dimnames(m) <- list(ids, paste0("v", off, "_", seq_len(p)))
    m
  }
  X <- list(a = mk(1, 3), b = mk(2, 3), c = mk(3, 3))

  h5 <- tempfile(fileext = ".h5")
  on.exit({ unlink(h5); BigDataStatMeth::hdf5_close_all() }, add = TRUE)

  before <- list.files(tempdir(), pattern = "\\.h5$")
  set.seed(11)
  pt <- suppressMessages(
    mgcca_permtest(X, filename = h5, nperm = 3L, nfac = 2L,
                   method = "penalized", lambda = rep(0.1, 3)))
  after <- list.files(tempdir(), pattern = "\\.h5$")

  expect_length(pt$eigenvalues, 2L)
  expect_equal(dim(pt$null), c(3L, 2L))
  expect_equal(pt$nperm, 3L)
  expect_true(all(pt$p.value >= 1 / 4 & pt$p.value <= 1))

  # only the observed fit's file is left behind; the three permutation files
  # used to stay in tempdir() for the rest of the session
  expect_setequal(setdiff(after, before), basename(h5))
})
