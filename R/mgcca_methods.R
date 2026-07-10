# S3 methods that make an 'mgcca' object behave like a standard R model object
# (as returned by mgcca_results() / mgcca(collect = TRUE)).

#' Print an mgcca object
#'
#' @param x an \code{mgcca} object.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
#' @method print mgcca
print.mgcca <- function(x, ...) {
    ev <- as.numeric(x$AVE$AVE_inner_model)
    cat("Generalized Canonical Correlation Analysis (mgcca)\n")
    cat(sprintf("  tables      : %d  (%s)\n", length(x$corsY),
                paste(names(x$corsY), collapse = ", ")))
    cat(sprintf("  individuals : %d\n", nrow(x$Y)))
    cat(sprintf("  components  : %d\n", length(ev)))
    cat(sprintf("  eigenvalues : %s\n", paste(sprintf("%.4g", ev), collapse = ", ")))
    cat(sprintf("  AVE (outer) : %s\n",
                paste(sprintf("%.4g", x$AVE$AVE_outer_model), collapse = ", ")))
    if (!is.null(x$scores))
        cat("  scores      : available (per table)\n")
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
#' @return A \code{ggplot} object.
#' @export
#' @method plot mgcca
plot.mgcca <- function(x, ...) plotIndividuals(x, ...)

#' Summarise an mgcca object
#'
#' @description Prints the model overview, the variance explained per table
#'   (outer-model AVE) and, for each table and component, the \code{top} variables
#'   most correlated with the shared components.
#'
#' @param object an \code{mgcca} object.
#' @param top number of top variables per component to list. Default 5.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @export
#' @method summary mgcca
summary.mgcca <- function(object, top = 5, ...) {
    print(object)
    A <- as.matrix(object$AVE$AVE_X)                 # components x tables
    dimnames(A) <- list(paste0("comp", seq_len(nrow(A))), names(object$corsY))
    cat("\nAVE per table (outer model):\n")
    print(round(t(A), 4))                            # tables x components

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
