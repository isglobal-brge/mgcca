#' Per-block contribution profile of a sensitivity result
#'
#' @description Reformats the block-level decomposition that
#'   \code{\link{mgcca_sensitivity}} already computed: how much of a query's
#'   first-order sensitivity each block carries, and how much of that is carried
#'   by the grouping.
#'
#' @details \strong{This is an extractor, not a second analysis.} It recomputes
#'   nothing. The block-level quantities come back from the same kernel pass that
#'   produced the totals, so there is exactly one route to them; a second
#'   independent implementation would eventually disagree with the first, and the
#'   disagreement would be silent.
#'
#'   \strong{Shares are relative within a query.} A block carrying most of a
#'   negligible sensitivity still carries most of it, so the absolute columns are
#'   returned alongside and should be read with the shares.
#'
#' @param x An object of class \code{"mgcca_sensitivity"}.
#' @param query Optional query name, or names, to restrict to. Default: all.
#' @return An object of class \code{"mgcca_block_profile"}: a data frame with one
#'   row per query and block, plus the originating provenance.
#' @seealso \code{\link{mgcca_sensitivity}}
#' @examples
#' ## mgcca_sensitivity() needs the HDF5-backed source blocks of a fit, so the
#' ## example works from a minimal stand-in carrying the one element this
#' ## extractor reads, `by_block` -- exactly the shape mgcca_sensitivity()
#' ## returns for an ungrouped analysis of three blocks and two queries.
#' sens <- structure(
#'     list(by_block = data.frame(
#'              query = rep(c("age", "bmi"), each = 3),
#'              block = rep(c("methylation", "clinical", "other"), 2),
#'              S_total = c(0.62, 0.25, 0.13, 0.31, 0.44, 0.25),
#'              share_total = c(0.62, 0.25, 0.13, 0.31, 0.44, 0.25),
#'              stringsAsFactors = FALSE),
#'          groups = NULL,
#'          settings = list(grouped = FALSE),
#'          provenance = list(file = NA_character_, datasets = NULL),
#'          call = NULL),
#'     class = "mgcca_sensitivity")
#'
#' bp <- mgcca_block_profile(sens)
#' bp$profile
#'
#' ## restrict to one query
#' mgcca_block_profile(sens, query = "age")$profile
#' @export
mgcca_block_profile <- function(x, query = NULL) {
    if (!inherits(x, "mgcca_sensitivity"))
        stop("`x` must be the result of mgcca_sensitivity(). ",
             "This function extracts a decomposition that has already been computed; ",
             "it does not compute one.", call. = FALSE)
    d <- x$by_block
    if (is.null(d) || !nrow(d))
        stop("this sensitivity result carries no per-block decomposition", call. = FALSE)
    if (!is.null(query)) {
        unknown <- setdiff(query, unique(d$query))
        if (length(unknown))
            stop("unknown query name(s): ", paste(unknown, collapse = ", "),
                 ". Available: ", paste(unique(d$query), collapse = ", "), call. = FALSE)
        d <- d[d$query %in% query, , drop = FALSE]
    }
    res <- list(profile = d, groups = x$groups, grouped = isTRUE(x$settings$grouped),
                provenance = x$provenance, from = x$call)
    class(res) <- "mgcca_block_profile"
    res
}

#' Print a block contribution profile
#'
#' @description Prints the per-query, per-block table built by
#'   \code{\link{mgcca_block_profile}}, preceded by a one-line header with the
#'   number of queries and blocks and the grouping the sensitivity was run with.
#'   When no grouping was supplied, a closing note says so: the between-group
#'   columns are undefined there, not missing.
#'
#' @param x An object of class \code{"mgcca_block_profile"}.
#' @param ... ignored.
#' @return Called for its side effect (the profile is written to the console).
#'   Returns \code{x} invisibly.
#' @seealso \code{\link{mgcca_block_profile}}
#' @examples
#' ## A minimal stand-in for an ungrouped mgcca_sensitivity() result; see
#' ## mgcca_block_profile() for why the example does not fit one here.
#' sens <- structure(
#'     list(by_block = data.frame(
#'              query = rep(c("age", "bmi"), each = 3),
#'              block = rep(c("methylation", "clinical", "other"), 2),
#'              S_total = c(0.62, 0.25, 0.13, 0.31, 0.44, 0.25),
#'              share_total = c(0.62, 0.25, 0.13, 0.31, 0.44, 0.25),
#'              stringsAsFactors = FALSE),
#'          groups = NULL,
#'          settings = list(grouped = FALSE),
#'          provenance = list(file = NA_character_, datasets = NULL),
#'          call = NULL),
#'     class = "mgcca_sensitivity")
#'
#' print(mgcca_block_profile(sens))
#' @export
#' @method print mgcca_block_profile
print.mgcca_block_profile <- function(x, ...) {
    cat("mgcca block profile\n")
    cat(sprintf("  queries : %d ; blocks : %d ; grouping : %s\n",
                length(unique(x$profile$query)), length(unique(x$profile$block)),
                if (x$grouped) paste(x$groups, collapse = ", ") else "none"))
    cat("\n")
    print(format(x$profile, digits = 4), row.names = FALSE)
    if (!x$grouped)
        cat("\n  No grouping was supplied to mgcca_sensitivity(), so the between-group\n",
            "  columns are not shown: they are undefined here, not missing.\n", sep = "")
    invisible(x)
}

#' Plot a block contribution profile
#'
#' @description Bar chart of \code{share_total}, the share of a query's
#'   first-order sensitivity each block carries, with one panel per query and the
#'   blocks ordered by their share. The shares are relative within a query, so
#'   read them next to the absolute \code{S_total} column of
#'   \code{x$profile}: a block carrying most of a negligible sensitivity still
#'   carries most of it.
#'
#' @param x An object of class \code{"mgcca_block_profile"}.
#' @param ... ignored.
#' @return A \code{ggplot} object. Drawn when printed.
#' @seealso \code{\link{mgcca_block_profile}}
#' @examples
#' ## A minimal stand-in for an ungrouped mgcca_sensitivity() result; see
#' ## mgcca_block_profile() for why the example does not fit one here.
#' sens <- structure(
#'     list(by_block = data.frame(
#'              query = rep(c("age", "bmi"), each = 3),
#'              block = rep(c("methylation", "clinical", "other"), 2),
#'              S_total = c(0.62, 0.25, 0.13, 0.31, 0.44, 0.25),
#'              share_total = c(0.62, 0.25, 0.13, 0.31, 0.44, 0.25),
#'              stringsAsFactors = FALSE),
#'          groups = NULL,
#'          settings = list(grouped = FALSE),
#'          provenance = list(file = NA_character_, datasets = NULL),
#'          call = NULL),
#'     class = "mgcca_sensitivity")
#'
#' p <- plot(mgcca_block_profile(sens))
#' class(p)
#' @export
#' @method plot mgcca_block_profile
plot.mgcca_block_profile <- function(x, ...) {
    d <- x$profile
    ggplot2::ggplot(d, ggplot2::aes(x = stats::reorder(.data$block, .data$share_total),
                                    y = .data$share_total)) +
        ggplot2::geom_col(width = 0.6, fill = .mgcca_palette[1]) +
        ggplot2::facet_wrap(~ .data$query) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
        ggplot2::labs(x = NULL, y = "share of the total sensitivity",
                      title = "Block contribution",
                      subtitle = "shares are within a query; read them with the absolute values") +
        .mgcca_theme()
}
