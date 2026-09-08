#' @export
print.mgcca_sensitivity <- function(x, ...) {
    s <- x$settings
    cat("mgcca sensitivity\n")
    cat(sprintf("  queries   : %d\n", nrow(x$overall)))
    cat(sprintf("  blocks    : %d\n", length(unique(x$by_block$block))))
    cat(sprintf("  grouping  : %s\n",
                if (isTRUE(s$grouped))
                    paste0(length(x$groups), " levels (",
                           paste(utils::head(x$groups, 4), collapse = ", "),
                           if (length(x$groups) > 4) ", ..." else "", ")")
                else "none -- total sensitivity only"))
    # Results produced before the per-individual decomposition existed carry no
    # `by_individual`; they print exactly as they always did.
    if (!is.null(x[["by_individual"]]))
        cat(sprintf("  individuals: %d per query (by_individual, exploratory)\n",
                    length(unique(x[["by_individual"]]$id))))
    cat(sprintf("  backend   : %s ; L = %d ; ridge = %s\n",
                s$backend, s$L, s$lambda_source))
    cat("\n")
    print(format(x$overall, digits = 4), row.names = FALSE)
    if (!isTRUE(x$validity$valid))
        cat("\n  ! the retained and discarded subspaces are not separated (",
            x$validity$reason, "). The first-order expansion is not meaningful\n",
            "    here, and these sensitivities must not be read as if it were.\n",
            sep = "")
    invisible(x)
}

#' @export
summary.mgcca_sensitivity <- function(object, ...) {
    print(object)
    cat("\nPer-block decomposition\n")
    print(format(object$by_block, digits = 4), row.names = FALSE)
    # Older results carry no per-individual table; the section is skipped in
    # silence rather than announced as missing.
    bi <- object[["by_individual"]]
    if (!is.null(bi) && nrow(bi)) {
        cat("\nPer-individual sensitivity contributions (exploratory)\n")
        for (q in unique(bi$query)) {
            d <- bi[bi$query == q, , drop = FALSE]
            s <- rowsum(d$S_total, d$id)
            # ⚠️ NOT the table's `share_total`, which is a share WITHIN a block
            # and so does not add up across blocks. The displayed share is
            # computed here against the query's own total, which is the only
            # denominator that makes a summed-over-blocks number a share at all.
            St <- object$overall$S_total[object$overall$query == q]
            sh <- if (length(St) == 1L && is.finite(St) && St > 0)
                s[, 1L] / St else rep(NA_real_, nrow(s))
            ord <- order(s[, 1L], decreasing = TRUE)
            take <- utils::head(ord, 5L)
            top <- data.frame(id = rownames(s)[take], S = s[take, 1L],
                              share_of_query = sh[take], stringsAsFactors = FALSE)
            cat(sprintf("  query '%s' -- top %d of %d participants by S_total\n",
                        q, nrow(top), nrow(s)))
            print(format(top, digits = 4), row.names = FALSE)
        }
        cat("\n  S is summed over blocks; share_of_query is its share of the query's\n")
        cat("  total sensitivity (the by_individual table's own share_total is a\n")
        cat("  share WITHIN one block, and sums to one there).\n")
        cat("  EXPLORATORY: the block-level parent has MIXTO calibration evidence,\n")
        cat("  and no individual-level calibration has been established. These are\n")
        cat("  sensitivity contributions -- not variance, uncertainty, causal\n")
        cat("  influence, or reliability guarantees -- conditioned on the grouped\n")
        cat("  participant set.\n")
    }
    cat("\n  T is the share of the query lying in the retained subspace.\n")
    cat("  S_total is its first-order sensitivity; S_between the part carried by the\n")
    cat("  grouping. Read R_between WITH its components: a high ratio can accompany a\n")
    cat("  negligible absolute sensitivity, and a low one a large signal.\n")
    invisible(object)
}

#' Plot a sensitivity result
#'
#' @description Four views. \code{"map"} places the individuals on the layer's
#'   reference axes and draws each query as an arrow, so a query that points along
#'   the direction separating your groups is visible rather than inferred from
#'   \code{R_between}. \code{"blocks"} (default) shows how much each block
#'   contributes to a query's sensitivity and how much of that is carried by the
#'   grouping; \code{"overall"} compares the queries; \code{"individuals"} cuts
#'   the same sensitivity the other way and shows which participants carry it.
#'
#' @details The absolute sensitivities are drawn, not only the shares, because a
#'   share cannot be interpreted without the magnitude it is a share of: a block
#'   can carry most of a sensitivity that is itself negligible.
#'
#'   \strong{Why the between-group part is overlaid rather than placed beside it.}
#'   Side by side it is routinely two orders of magnitude shorter than the total
#'   and renders as an invisible sliver, so the plot fails to show the one
#'   comparison it exists for. The between-group part is a component of the total,
#'   so it is drawn inside the same bar with its share printed above -- a sliver,
#'   however honest, cannot be read off an axis.
#'
#'   \strong{Why the \code{"individuals"} view always draws an "all others" bar.}
#'   A top-N chart without the remainder silently overstates concentration: fifteen
#'   tall bars look like the whole story whether they carry 80\% of the total or
#'   8\%. The remainder is therefore drawn as one further bar, stacked by block
#'   like the rest, so the panel accounts for the query's entire sensitivity and
#'   the reader can see how much the top N leaves out. It is the same rule that
#'   puts the shares next to the absolutes everywhere else in this package.
#'
#'   The per-participant quantities it draws are an exploratory decomposition:
#'   the block-level parent has empirical calibration evidence with a MIXTO
#'   verdict, no individual-level calibration has been established, and the bars
#'   are sensitivity contributions -- not variance, uncertainty, causal
#'   influence, or reliability guarantees. They are conditioned on the grouped
#'   participant set, so a different grouping is a different question.
#'
#' @param x An object of class \code{"mgcca_sensitivity"}.
#' @param type \code{"blocks"} (default), \code{"overall"}, \code{"map"} or
#'   \code{"individuals"}.
#' @param top_n For \code{type = "individuals"}, how many participants to draw
#'   per query, ranked by their sensitivity contribution summed over blocks.
#'   Default 15. Everyone else is aggregated into a single "all others" bar
#'   rather than dropped.
#' @param noise For \code{type = "overall"}, the names of the control queries --
#'   columns of random numbers put through the same machinery. Their range is
#'   shaded as the level chance alone reaches, which is the only scale
#'   \code{T} has.
#' @param comps For \code{type = "map"}, which two components to draw. With more
#'   than two components, call it repeatedly on different pairs: three dimensions
#'   need rotating to be read and do not survive a printed page, and beyond three
#'   there is nothing to rotate.
#' @param ... Ignored.
#' @return A \code{ggplot} object, which draws when printed.
#' @seealso \code{\link{mgcca_sensitivity}}, which produces the object plotted
#'   here.
#' @examples
#' ## A small three-block slice of the shipped cohort; the fit must be
#' ## file-backed for the reliability layer to reach the source blocks.
#' data(cardiovascular, package = "mgcca")
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:120]
#' mk  <- function(d, j = NULL) {
#'     m <- as.matrix(d[rownames(d) %in% ids, , drop = FALSE])
#'     if (!is.null(j)) m <- m[, j, drop = FALSE]
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = mk(X1, 1:10), clinical = mk(X2), cells = mk(X3))
#'
#' fit <- mgcca(X, filename = tempfile(fileext = ".h5"), nfac = 2,
#'              method = "penalized", lambda = rep(0.1, 3), scores = TRUE)
#'
#' bmi <- X$clinical[, "BMI", drop = FALSE]
#' cd4 <- X$cells[, "CD4T"]
#' grp <- stats::setNames(
#'     ifelse(cd4 > stats::median(cd4), "CD4T high", "CD4T low"),
#'     rownames(X$cells))
#' sens <- mgcca_sensitivity(fit, query = bmi, group = grp)
#'
#' ## Where the sensitivity sits, block by block, with the between-group part
#' ## drawn inside the total rather than beside it.
#' plot(sens, type = "blocks")
#'
#' ## The same fit as a map: individuals on the layer's own reference axes,
#' ## the query as an arrow, and the dashed line showing how the groups split.
#' plot(sens, type = "map")
#'
#' ## And cut the other way: which participants carry it. The "all others" bar
#' ## is always drawn, so the panel accounts for the whole total.
#' plot(sens, type = "individuals", top_n = 8)
#' @export
plot.mgcca_sensitivity <- function(x, type = c("blocks", "overall", "map",
                                              "individuals"),
                                  comps = c(1, 2), noise = NULL, top_n = 15, ...) {
    type <- match.arg(type)
    grouped <- isTRUE(x$settings$grouped)

    if (identical(type, "individuals")) {
        # ⚠️ `[[`, not `$`. On a result that predates the per-individual table,
        # `x$by_individual` PARTIALLY MATCHES `by_individual_check` and hands
        # back a list instead of NULL, so the refusal below never fires and the
        # failure surfaces much later as an unreadable type error.
        bi <- x[["by_individual"]]
        if (is.null(bi) || !nrow(bi))
            stop("this result carries no per-individual decomposition; it was ",
                 "produced by an older version of mgcca_sensitivity()", call. = FALSE)
        if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
            top_n < 1)
            stop("`top_n` must be a single number of participants to draw, ",
                 "at least 1", call. = FALSE)
        top_n <- as.integer(top_n)
        gmap <- x$reference$group
        blocks <- unique(bi$block)
        # ⚠️ ONE DISCRETE AXIS, DIFFERENT ORDERS PER PANEL. Participants rank
        # differently for different queries, and a single factor cannot hold two
        # orders at once. The axis value is therefore query-prefixed so the levels
        # can be ordered per query, `scales = "free_y"` drops the levels that do
        # not belong to a panel, and the prefix is stripped again by the label
        # function. The separator is a control character, so it cannot collide
        # with a participant identifier.
        SEP <- ""                 # unit separator: not a legal id character
        built <- lapply(split(bi, bi$query), function(d) {
            tot <- rowsum(d$S_total, d$id)[, 1L]
            keep <- utils::head(names(sort(tot, decreasing = TRUE)), top_n)
            other <- setdiff(names(tot), keep)
            named <- if (is.null(gmap)) d$id else {
                g <- unname(gmap[d$id])
                ifelse(is.na(g), d$id, paste0(d$id, " (", g, ")"))
            }
            # The remainder bar carries no group suffix: it is not one
            # participant, so a group label on it would be a false statement.
            oth_lab <- sprintf("all others (%d)", length(other))
            d$label <- ifelse(d$id %in% other, oth_lab, named)
            a <- stats::aggregate(list(S_total = d$S_total),
                                  list(label = d$label, block = d$block), sum)
            a$query <- d$query[1L]
            a$is_other <- a$label == oth_lab & length(other) > 0L
            # Largest at the top, and the remainder pinned to the bottom, where it
            # reads as the baseline the top N is measured against.
            ltot <- vapply(split(a$S_total, a$label), sum, numeric(1))
            if (length(other)) ltot[oth_lab] <- -Inf
            list(d = a, lev = paste0(a$query[1L], SEP, names(sort(ltot))))
        })
        dd <- do.call(rbind, lapply(built, `[[`, "d"))
        dd$key <- factor(paste0(dd$query, SEP, dd$label),
                         levels = unlist(lapply(built, `[[`, "lev"),
                                         use.names = FALSE))
        # The remainder is drawn lighter so the eye separates it from the named
        # participants without spending a second colour scale on it.
        dd$alpha <- ifelse(dd$is_other, 0.55, 1)
        # Stacked in the legend's order rather than against it: ggplot2 stacks
        # the last level first by default, so a reader matching colours to the
        # key reads the bar backwards.
        dd$block <- factor(dd$block, levels = blocks)
        return(ggplot2::ggplot(dd,
                ggplot2::aes(x = .data$S_total, y = .data$key,
                             fill = .data$block, alpha = .data$alpha)) +
            ggplot2::geom_col(width = 0.72,
                              position = ggplot2::position_stack(reverse = TRUE)) +
            ggplot2::scale_alpha_identity() +
            # A single query does not pay for a two-column layout: ncol = 2
            # would leave half the width blank.
            ggplot2::facet_wrap(~ .data$query, scales = "free_y",
                                ncol = min(2L, length(unique(dd$query)))) +
            ggplot2::scale_y_discrete(
                labels = function(v) sub(paste0("^.*", SEP), "", v)) +
            ggplot2::scale_fill_manual(
                values = stats::setNames(
                    rep(.mgcca_palette, length.out = length(blocks)), blocks),
                breaks = blocks, name = NULL) +
            ggplot2::scale_x_continuous(
                expand = ggplot2::expansion(mult = c(0, 0.06))) +
            ggplot2::labs(x = "first-order sensitivity contribution", y = NULL,
                          title = "Which participants carry the sensitivity",
                          subtitle = paste0(
                              "exploratory decomposition, conditioned on the ",
                              "grouped participant set\n",
                              "no individual-level calibration has been established"),
                          caption = sprintf(paste("top %d per query, stacked by",
                              "block; the 'all others' bar is everyone else, so",
                              "each panel accounts for the whole total"), top_n)) +
            .mgcca_theme() +
            ggplot2::theme(legend.position = "bottom",
                           plot.caption = ggplot2::element_text(
                               size = 7.5, colour = "grey40", hjust = 0)))
    }

    if (identical(type, "blocks")) {
        d <- x$by_block
        p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$block)) +
            ggplot2::geom_col(ggplot2::aes(y = .data$S_total, fill = "total"),
                              width = 0.68)
        if (grouped && "S_between" %in% names(d)) {
            p <- p +
                ggplot2::geom_col(ggplot2::aes(y = .data$S_between,
                                               fill = "between groups"),
                                  width = 0.34) +
                ggplot2::geom_text(
                    ggplot2::aes(y = .data$S_total,
                                 label = sprintf("%.1f%%",
                                                 100 * .data$S_between / .data$S_total)),
                    vjust = -0.45, size = 3, colour = "grey25")
        }
        p <- p +
            ggplot2::facet_wrap(~ .data$query, scales = "free_y", ncol = 2) +
            ggplot2::scale_fill_manual(
                values = stats::setNames(.mgcca_palette[seq_len(2)],
                                         c("total", "between groups")),
                breaks = c("total", "between groups"), name = NULL) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0, 0.16))) +
            ggplot2::labs(x = NULL, y = "first-order sensitivity",
                          title = "Sensitivity by block",
                          subtitle = if (grouped)
                              "bar = total; inner bar and label = the part carried by the grouping"
                          else "no grouping supplied: total sensitivity only") +
            .mgcca_theme() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
        return(p)
    }

    if (identical(type, "map")) {
        R <- x$reference
        if (is.null(R) || is.null(R$Y))
            stop("this result carries no reference subspace; it was produced by an ",
                 "older version of mgcca_sensitivity()", call. = FALSE)
        L <- ncol(R$Y)
        if (max(comps) > L || min(comps) < 1 || length(comps) != 2)
            stop("`comps` must be two component numbers between 1 and ", L, call. = FALSE)
        Y <- as.data.frame(R$Y[, comps, drop = FALSE])
        names(Y) <- c("a", "b")
        Y$group <- if (is.null(R$group)) NA_character_ else
            unname(R$group[match(rownames(R$Y), names(R$group))])
        keep <- !is.na(Y$group) | is.null(R$group)
        Y <- Y[keep, , drop = FALSE]

        # Arrows are scaled to the cloud, not to their own length: their DIRECTION
        # is the message, and their length in the original units would be
        # meaningless next to individual coordinates.
        W <- R$directions[comps, , drop = FALSE]
        nrmW <- sqrt(colSums(W^2)); nrmW[nrmW == 0] <- 1
        sc <- 0.80 * max(abs(unlist(Y[, c("a", "b")])), na.rm = TRUE)
        A <- data.frame(lab = colnames(W),
                        xend = W[1, ] / nrmW * sc, yend = W[2, ] / nrmW * sc,
                        stringsAsFactors = FALSE)

        p <- ggplot2::ggplot(Y, ggplot2::aes(.data$a, .data$b))
        if (!is.null(R$group)) {
            p <- p +
                ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.9, alpha = 0.75) +
                ggplot2::stat_ellipse(ggplot2::aes(colour = .data$group), level = 0.68,
                                      linewidth = 0.4) +
                ggplot2::scale_colour_manual(values = .mgcca_palette, name = NULL)
            cen <- stats::aggregate(cbind(a, b) ~ group, Y, mean)
            cen <- cen[order(cen$group), , drop = FALSE]
            # ⚠️ DRAW THE THING THE ARROWS ARE COMPARED AGAINST, AND FROM THE SAME
            # ORIGIN THEY START AT. A path joining the centroids is the obvious
            # choice and it fails: the centroids sit almost on top of each other, so
            # the path is a few pixels long while the arrows span the panel, and the
            # eye ends up connecting an arrow tip to the nearest coloured cloud
            # instead -- which reverses the reading entirely. Instead the direction
            # from the first to the last centroid is drawn through the origin at the
            # same scale as the arrows, so the comparison is one of ANGLE, which is
            # what the quantity actually is.
            gdir <- as.numeric(cen[nrow(cen), c("a", "b")] - cen[1L, c("a", "b")])
            gn <- sqrt(sum(gdir^2))
            if (is.finite(gn) && gn > 0) {
                gdir <- gdir / gn * sc
                p <- p + ggplot2::annotate("segment",
                        x = -gdir[1], y = -gdir[2], xend = gdir[1], yend = gdir[2],
                        linetype = "22", linewidth = 0.6, colour = "grey45") +
                    # Label at the end that points RIGHT, right-aligned so the text
                    # runs back INTO the panel. Centring it on an endpoint clips it
                    # against whichever edge that endpoint happens to fall near, and
                    # which end that is depends on the data.
                    # ...and on a white plate, because wherever it lands it lands on
                    # top of the individuals it is describing.
                    ggplot2::annotate("label",
                        x = if (gdir[1] >= 0) gdir[1] else -gdir[1],
                        y = if (gdir[1] >= 0) gdir[2] else -gdir[2],
                        label = "how the groups separate", hjust = 1, vjust = 1.25,
                        size = 2.8, colour = "grey45", fill = "white",
                        alpha = 0.85, label.size = 0, label.padding =
                            ggplot2::unit(0.10, "lines"))
            }
            p <- p + ggplot2::geom_point(data = cen, ggplot2::aes(colour = .data$group),
                                         size = 4, shape = 18)
        } else {
            p <- p + ggplot2::geom_point(size = 1.9, alpha = 0.7,
                                         colour = .mgcca_palette[1])
        }
        return(p +
            ggplot2::geom_segment(data = A,
                ggplot2::aes(x = 0, y = 0, xend = .data$xend, yend = .data$yend),
                arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "closed"),
                linewidth = 0.7, colour = "grey15", inherit.aes = FALSE) +
            ggplot2::geom_text(data = A,
                ggplot2::aes(x = .data$xend, y = .data$yend, label = .data$lab),
                hjust = -0.05, vjust = -0.4, size = 3.1, colour = "grey15",
                inherit.aes = FALSE) +
            ggplot2::coord_equal(clip = "off") +
            ggplot2::labs(x = paste("component", comps[1]),
                          y = paste("component", comps[2]),
                          title = "Groups, and where each query points",
                          subtitle = if (is.null(R$group))
                              "arrows = the direction each query represents"
                          else "arrow along the dashed line = it moves with the groups") +
            .mgcca_theme() +
            # Legend underneath: the arrow labels live on the right of the panel and
            # collided with it. Moving the legend is cheaper than shortening the
            # labels, which are what makes the figure readable in the first place.
            ggplot2::theme(legend.position = "bottom",
                           legend.margin = ggplot2::margin(t = -4),
                           plot.margin = ggplot2::margin(6, 30, 6, 6)) +
            ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)))
    }

    o <- x$overall
    # ⚠️ The noise band is what makes T readable. Without it the reader has a
    # number with no scale: T has no conventional threshold, and its floor depends
    # on the sample size and on how many axes were kept. Naming the control
    # queries shades that floor, so "is this above chance?" becomes a look rather
    # than a calculation.
    band <- NULL
    if (!is.null(noise)) {
        unknown <- setdiff(noise, o$query)
        if (length(unknown))
            stop("`noise` names queries that are not in this result: ",
                 paste(unknown, collapse = ", "), call. = FALSE)
        band <- range(o$T[o$query %in% noise])
        o <- o[!o$query %in% noise, , drop = FALSE]
        if (!nrow(o)) stop("every query was named as noise", call. = FALSE)
    }
    ord <- o$query[order(o$T)]
    o$query <- factor(o$query, levels = ord)
    p <- ggplot2::ggplot(o, ggplot2::aes(x = .data$query, y = .data$T))
    if (!is.null(band))
        p <- p + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                                   ymin = band[1], ymax = band[2],
                                   fill = "grey70", alpha = 0.35)
    p <- p +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$query, y = 0, yend = .data$T),
                              colour = "grey75", linewidth = 0.4) +
        ggplot2::geom_point(size = 2.6, colour = .mgcca_palette[1]) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "alignment with the shared axes",
                      title = "How much of each query the shared signal accounts for",
                      subtitle = if (is.null(band))
                          "no control queries named: pass `noise =` to shade the chance level"
                      else "grey band = the level reached by the control queries, i.e. chance") +
        .mgcca_theme()
    p
}
