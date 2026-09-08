# S3 methods that make an 'mgcca' object behave like a standard R model object
# (as returned by mgcca_results() / mgcca(collect = TRUE)).

#' Print an mgcca object
#'
#' @description One-screen overview of the fit. The \code{spectral gap} line is
#'   the external gap \eqn{\mu_L - \mu_{L+1}} of the mgcca operator, relative to
#'   \eqn{\mu_L}: it says whether the leading \eqn{L}-dimensional subspace is
#'   numerically identified at all. A gap at machine zero is flagged
#'   \code{NOT IDENTIFIED} -- the solver converged, but the directions it
#'   returned are arbitrary (see \code{\link{mgcca}}).
#'
#' @param x an \code{mgcca} object.
#' @param ... ignored.
#' @return Called for its side effect (the overview is written to the console).
#'   Returns \code{x} invisibly.
#' @examples
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:150]
#' num <- function(d, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% ids, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = num(X1, 1:20), clinical = num(X2), other = num(X3))
#'
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3),
#'              scores = TRUE)
#' print(fit)
#' @export
#' @method print mgcca
print.mgcca <- function(x, ...) {
    # An object collected with mgcca(outputs = ...) carries only what was asked
    # for; every block below is printed if and only if its component is there.
    # A fully collected object prints exactly what it always printed.
    ev <- if (!is.null(x$AVE)) as.numeric(x$AVE$AVE_inner_model) else numeric(0)
    L  <- if (length(ev)) length(ev)
          else if (!is.null(x$Y)) ncol(x$Y) else NA_integer_
    cat("Generalized Canonical Correlation Analysis (mgcca)\n")
    if (!is.null(x$corsY))
        cat(sprintf("  tables      : %d  (%s)\n", length(x$corsY),
                    paste(names(x$corsY), collapse = ", ")))
    if (!is.null(x$Y))
        cat(sprintf("  individuals : %d\n", nrow(x$Y)))
    if (!is.null(x$AVE)) {
        cat(sprintf("  components  : %d\n", length(ev)))
        cat(sprintf("  eigenvalues : %s\n", paste(sprintf("%.4g", ev), collapse = ", ")))
        cat(sprintf("  AVE (outer) : %s\n",
                    paste(sprintf("%.4g", x$AVE$AVE_outer_model), collapse = ", ")))
    }
    eg <- x$eigen
    if (!is.null(eg) && is.finite(eg$gap_rel) && !is.na(L))
        cat(sprintf("  spectral gap: %.4g (relative, external at L = %d)%s\n",
                    eg$gap_rel, L,
                    if (.mgcca_subspace_degenerate(eg))
                        "   <- NOT IDENTIFIED" else ""))
    if (!is.null(x$scores))
        cat("  scores      : available (per table)\n")
    # An in-memory fit (mgcca(filename = NULL)) is not written to disk, so it
    # cannot be reloaded with mgcca_load() and the reliability layer cannot run on
    # it. File-backed fits print exactly as before (their desc has no backend).
    if (identical(attr(x, "desc")$backend, "memory"))
        cat("  backend     : in-memory fit -- not reloadable from a file\n")
    invisible(x)
}

#' Plot an mgcca object
#'
#' @description Default plot for an \code{mgcca} object: the individuals on the
#'   two shared components. A thin wrapper around \code{\link{plotIndividuals}};
#'   all its arguments (\code{group}, \code{comps}, \code{ellipse}, ...) are
#'   forwarded.
#'
#' @param x an \code{mgcca} object.
#' @param ... passed to \code{\link{plotIndividuals}}.
#' @return A \code{ggplot} object: the individuals of \code{x} on two shared
#'   components. Drawn when printed.
#' @examples
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:150]
#' num <- function(d, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% ids, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = num(X1, 1:20), clinical = num(X2), other = num(X3))
#'
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3),
#'              scores = TRUE)
#' p <- plot(fit)
#' class(p)
#' @export
#' @method plot mgcca
plot.mgcca <- function(x, ...) plotIndividuals(x, ...)

#' Summarise an mgcca object
#'
#' @description Prints the model overview, the variance explained per table
#'   (outer-model AVE), the pairwise overlap between tables and, for each table
#'   and component, the \code{top} variables most correlated with the shared
#'   components.
#'
#' @details The overlap block reports, for every pair of tables \eqn{(j,k)}, the
#'   number \eqn{n_{jk}} of individuals both tables hold and its normalised
#'   version \eqn{\alpha_{jk} = n_{jk}/n} over the union, plus the dispersion of
#'   \eqn{\alpha} across pairs. It is descriptive: nothing in the fit depends on
#'   it. It is worth reading because it is the structural quantity that decides
#'   how much the treatment of missing individuals can matter -- with
#'   \eqn{\alpha} flat and near 1 every sensible treatment agrees, and the
#'   choice only starts to bite when the pairs disagree about who they share.
#'   Two notes are printed when they apply: one when the dispersion of
#'   \eqn{\alpha} is large, one when some pair shares very few individuals.
#'   Both are reporting heuristics with no inferential status.
#'
#' @param object an \code{mgcca} object.
#' @param top number of top variables per component to list. Default 5.
#' @param ... ignored.
#' @return Called for its side effect (the overview, the per-table AVE, the
#'   pairwise overlap and the top variables are written to the console). Returns
#'   \code{object} invisibly.
#' @examples
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:150]
#' num <- function(d, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% ids, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = num(X1, 1:20), clinical = num(X2), other = num(X3))
#'
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3),
#'              scores = TRUE)
#' summary(fit, top = 3)
#' @export
#' @method summary mgcca
summary.mgcca <- function(object, top = 5, ...) {
    print(object)
    if (!is.null(object$AVE) && !is.null(object$corsY)) {
        A <- as.matrix(object$AVE$AVE_X)             # components x tables
        dimnames(A) <- list(paste0("comp", seq_len(nrow(A))), names(object$corsY))
        cat("\nAVE per table (outer model):\n")
        print(round(t(A), 4))                        # tables x components
    }

    ov <- object$overlap
    if (!is.null(ov)) {
        cat(sprintf(
            "\nPairwise overlap (n = %d individuals in the union):\n", ov$n))
        tab <- data.frame(
            n_jk     = ov$pairs$n_jk,
            alpha_jk = round(ov$pairs$alpha_jk, 4),
            row.names = paste(ov$pairs$table1, ov$pairs$table2, sep = "-"))
        print(tab)
        cat(sprintf(
            "  n_j: %s\n",
            paste(sprintf("%s = %d", names(ov$n_present), ov$n_present),
                  collapse = ", ")))
        cat(sprintf(
            paste0("  alpha_jk = n_jk / n;  sd(alpha) = %.4f (cv = %.3f),",
                   "  range %.4f - %.4f\n"),
            ov$alpha_sd, ov$alpha_cv, ov$alpha_min, ov$alpha_max))
        if (isTRUE(ov$heterogeneous))
            cat("  Note: the overlap is heterogeneous across pairs -- this is",
                "the regime\n        in which the treatment of missing",
                "individuals matters most.\n")
        if (isTRUE(ov$sparse_pair))
            cat(sprintf(
                paste0("  Note: the sparsest pair shares only %d individuals;",
                       " the association between\n        those two tables",
                       " rests on very few individuals.\n"), ov$n_jk_min))
    }

    if (is.null(object$corsY)) return(invisible(object))
    cat(sprintf("\nTop %d variables per component (correlation with Y):\n", top))
    for (tb in names(object$corsY)) {
        C <- as.matrix(object$corsY[[tb]])
        cat(sprintf("  [%s]\n", tb))
        for (k in seq_len(ncol(C))) {
            o  <- order(abs(C[, k]), decreasing = TRUE)[seq_len(min(top, nrow(C)))]
            vs <- paste(sprintf("%s (%+.2f)", rownames(C)[o], C[o, k]),
                        collapse = ", ")
            cat(sprintf("    comp%d: %s\n", k, vs))
        }
    }
    invisible(object)
}

#' Project new individuals onto the shared canonical space
#'
#' @description Applies the per-table weights learned by \code{\link{mgcca}} to
#'   new individuals, returning their canonical scores per table. The new tables
#'   are centred/scaled with the \emph{training} parameters, so the projection is
#'   consistent with the fit. Requires a fit produced with \code{scores = TRUE}
#'   (weights) and \code{scale = TRUE} (scaling parameters).
#'
#' @param object an \code{mgcca} object.
#' @param newdata a named \code{list} of new tables (individuals x variables) with
#'   the same variable columns used at fit time; tables not seen at fit time are
#'   ignored. Only the tables present are projected.
#' @param ... ignored.
#' @return A named \code{list} of score matrices (new individuals x components),
#'   one per projected table.
#' @seealso \code{\link{mgcca}}, \code{\link{plotScores}}
#' @examples
#' data(cardiovascular)
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))
#' num <- function(d, keep, cols = seq_len(ncol(d))) {
#'     m <- as.matrix(d[rownames(d) %in% keep, cols, drop = FALSE])
#'     storage.mode(m) <- "double"
#'     m
#' }
#' train <- ids[1:150]
#' test  <- ids[151:180]
#' X  <- list(methylation = num(X1, train, 1:20), clinical = num(X2, train),
#'            other = num(X3, train))
#' Xn <- list(methylation = num(X1, test, 1:20), clinical = num(X2, test),
#'            other = num(X3, test))
#'
#' ## weights (scores = TRUE) and scaling (scale = TRUE) are what a
#' ## projection needs; both are on by default except for scores.
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = rep(0.1, 3),
#'              scores = TRUE)
#' proj <- predict(fit, newdata = Xn)
#' vapply(proj, nrow, integer(1))
#' head(proj$clinical)
#' @export
#' @method predict mgcca
predict.mgcca <- function(object, newdata, ...) {
    if (is.null(object$weights))
        stop("this fit has no weights; refit with mgcca(..., scores = TRUE)")
    if (!is.list(newdata) || is.null(names(newdata)))
        stop("'newdata' must be a named list of tables")

    out <- list()
    for (nm in names(newdata)) {
        W <- object$weights[[nm]]
        if (is.null(W)) next                       # table not part of the fit
        vars <- rownames(W)
        Xn <- as.matrix(newdata[[nm]]); storage.mode(Xn) <- "double"
        miss <- setdiff(vars, colnames(Xn))
        if (length(miss))
            stop("newdata table '", nm, "' is missing ", length(miss),
                 " variable(s) used at fit time (e.g. ", miss[1], ")")
        Xn <- Xn[, vars, drop = FALSE]
        cs <- object$scaling[[nm]]
        if (!is.null(cs)) {                        # apply training centre/scale
            cs <- cs[vars, , drop = FALSE]
            sc <- cs[, "scale"]; sc <- ifelse(sc == 0, 1, sc)
            Xn <- sweep(sweep(Xn, 2, cs[, "center"], "-"), 2, sc, "/")
        }
        s <- Xn %*% W
        colnames(s) <- paste0("comp", seq_len(ncol(s)))
        out[[nm]] <- s
    }
    if (!length(out))
        stop("none of the tables in 'newdata' were used at fit time")
    out
}
