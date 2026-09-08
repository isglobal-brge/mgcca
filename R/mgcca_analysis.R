# Downstream analyses on an mgcca fit: associate the shared components with
# external phenotypes, and test their significance by permutation.

#' Associate the shared components with external phenotypes
#'
#' @description Regresses each shared component on one or more external variables
#'   (a phenotype, a clinical group, a batch, ...) and reports how much of the
#'   component each explains (\eqn{R^2}) with a p-value. This answers "do the
#'   integrated axes capture this variable?" -- useful to interpret an mgcca
#'   result against clinical or technical annotation.
#'
#' @param x an \code{mgcca} object (from \code{\link{mgcca}} /
#'   \code{\link{mgcca_results}}).
#' @param phenotypes a single vector (factor or numeric), a named \code{list} of
#'   such vectors, or a \code{data.frame}. Named vectors / data-frame row names
#'   are matched to the individuals of \code{x}; otherwise they are assumed to be
#'   in the same order.
#' @param comps components to test. Default: all.
#' @return A \code{data.frame} with one row per (phenotype, component):
#'   \code{phenotype}, \code{component}, \code{R2}, \code{p}, \code{n} (non-missing
#'   individuals) and \code{p.adj} (Benjamini-Hochberg across all rows).
#' @seealso \code{\link{plotIndividuals}}, \code{\link{mgcca_permtest}}
#' @examples
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:150]
#' num <- function(d, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% ids, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = num(X1, 1:20), clinical = num(X2), other = num(X3))
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3))
#'
#' ## An external annotation, matched to the individuals by row name. Here it
#' ## is simulated, so the components are not expected to capture it.
#' set.seed(1)
#' pheno <- data.frame(
#'     group = factor(sample(c("case", "control"), nrow(fit$Y), TRUE)),
#'     age   = rnorm(nrow(fit$Y), 55, 8),
#'     row.names = rownames(fit$Y))
#'
#' a <- mgcca_associate(fit, pheno)
#' a
#' a[a$p.adj < 0.05, ]
#' @export
mgcca_associate <- function(x, phenotypes, comps = NULL) {
    .check_mgcca(x)
    Y <- as.matrix(x$Y)
    if (is.null(comps)) comps <- seq_len(ncol(Y))
    if (!all(comps %in% seq_len(ncol(Y)))) stop("'comps' out of range")

    if (is.data.frame(phenotypes)) {
        rn <- rownames(phenotypes)
        ph <- lapply(as.list(phenotypes), function(v) { names(v) <- rn; v })
    } else if (is.list(phenotypes)) {
        ph <- phenotypes
    } else {
        ph <- list(phenotype = phenotypes)
    }
    if (is.null(names(ph)) || any(names(ph) == ""))
        names(ph) <- paste0("phenotype", seq_along(ph))

    align <- function(v) {
        if (!is.null(names(v)) && !is.null(rownames(Y))) v[rownames(Y)] else v
    }
    one <- function(y, g) {
        ok <- !is.na(g) & !is.na(y)
        n  <- sum(ok)
        if (is.character(g) || is.factor(g)) g <- droplevels(factor(g[ok])) else g <- g[ok]
        if ((is.factor(g) && nlevels(g) < 2) || n < 3)
            return(c(R2 = NA_real_, p = NA_real_, n = n))
        m <- stats::lm(y[ok] ~ g); f <- stats::summary.lm(m)$fstatistic
        c(R2 = unname(stats::summary.lm(m)$r.squared),
          p  = unname(stats::pf(f[1L], f[2L], f[3L], lower.tail = FALSE)),
          n  = n)
    }

    rows <- list()
    for (nm in names(ph)) {
        v <- align(ph[[nm]])
        for (cc in comps) {
            r <- one(Y[, cc], v)
            rows[[length(rows) + 1L]] <- data.frame(
                phenotype = nm, component = cc,
                R2 = r[["R2"]], p = r[["p"]], n = r[["n"]],
                stringsAsFactors = FALSE)
        }
    }
    out <- do.call(rbind, rows)
    out$p.adj <- stats::p.adjust(out$p, "BH")
    rownames(out) <- NULL
    out
}

#' Permutation test for the significance of the shared components
#'
#' @description Builds a null distribution for the mgcca eigenvalues by breaking
#'   the correspondence between tables (independently permuting each table's
#'   individuals) and re-fitting \code{nperm} times. A component whose observed
#'   eigenvalue exceeds the permuted ones carries genuine cross-table signal.
#'
#'   Most informative when \eqn{p < n} (eigenvalues below 1); for \eqn{p \gg n}
#'   with penalisation the eigenvalues sit near 1 by construction and the test is
#'   uninformative. Re-fits the model \code{nperm} times, so use it on modest data
#'   (or a feature subset).
#'
#' @param x input tables, exactly as passed to \code{\link{mgcca}} (a named list
#'   of matrices, a \code{MultiAssayExperiment}, ...).
#' @param filename target HDF5 file for the observed fit. Default a temp file.
#' @param nperm number of permutations. Default 99.
#' @param nfac,method,lambda,route,scale passed to \code{\link{mgcca}}.
#' @param ... further arguments to \code{\link{mgcca}}.
#' @return A list with \code{eigenvalues} (observed, descending), \code{p.value}
#'   (per component), the \code{null} matrix (\code{nperm} x \code{nfac}) and
#'   \code{nperm}.
#' @seealso \code{\link{mgcca_associate}}
#' @examples
#' ## Each permutation is a full file-backed re-fit, so even a token number of
#' ## them takes tens of seconds; the example is kept out of the check budget.
#' \donttest{
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:60]
#' num <- function(d, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% ids, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = num(X1, 1:4), clinical = num(X2, 1:4),
#'           other = num(X3, 1:4))
#'
#' h5 <- tempfile(fileext = ".h5")
#' pt <- mgcca_permtest(X, filename = h5, nperm = 9, nfac = 2,
#'                      method = "penalized", lambda = rep(0.1, 3))
#' pt$eigenvalues
#' pt$p.value          # nperm = 9 bounds these below at 0.1
#' unlink(h5)
#' }
#' @export
mgcca_permtest <- function(x, filename = tempfile(fileext = ".h5"), nperm = 99,
                           nfac = 2, method = "penalized", lambda = NULL,
                           route = "auto", scale = TRUE, ...) {
    fit0 <- mgcca(x, filename = filename, nfac = nfac, method = method,
                  lambda = lambda, route = route, scale = scale,
                  collect = FALSE, ...)
    obs <- sort(as.numeric(fit0$eig_values), decreasing = TRUE)

    tabs <- .mgcca_as_tables(x)                 # coerce once to a list of matrices
    null <- matrix(NA_real_, nperm, length(obs))
    for (b in seq_len(nperm)) {
        xb <- lapply(tabs, function(m) {
            rownames(m) <- sample(rownames(m)); m
        })
        fb <- mgcca(xb, filename = tempfile(fileext = ".h5"), nfac = nfac,
                    method = method, lambda = lambda, route = route,
                    scale = scale, collect = FALSE, ...)
        null[b, ] <- sort(as.numeric(fb$eig_values), decreasing = TRUE)
    }
    p <- vapply(seq_along(obs), function(k)
        (1 + sum(null[, k] >= obs[k])) / (nperm + 1), numeric(1))
    list(eigenvalues = obs, p.value = p, null = null, nperm = nperm)
}
