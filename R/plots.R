# Modern ggplot2-based result plots for mgcca objects (as returned by
# mgcca_results() / mgcca(collect = TRUE)). These replace the ageing base-R
# plotInds()/plotVars(): a single consistent, colourblind-safe look, ggplot
# objects you can further customise, and sensible handling of omics-scale tables
# (thousands of variables) via top-N labelling.

# Okabe-Ito colourblind-safe qualitative palette.
.mgcca_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                    "#56B4E9", "#F0E442", "#999999", "#000000")

.mgcca_theme <- function(base_size = 12) {
    ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = "grey90"),
            plot.title       = ggplot2::element_text(face = "bold"),
            legend.position  = "right",
            axis.title       = ggplot2::element_text(colour = "grey20"))
}

.check_mgcca <- function(x) {
    if (!inherits(x, "mgcca"))
        stop("'x' must be an 'mgcca' object (see mgcca_results() or ",
             "mgcca(collect = TRUE)).", call. = FALSE)
}

# Axis label for a single component. When the eigenvalue is available it is
# appended via plotmath (so the lambda symbol renders on every device, unlike a
# literal unicode character which warns on some graphics devices).
.comp_label <- function(comp, eig = NULL) {
    if (!is.null(eig) && length(eig) >= comp)
        bquote("Shared component" ~ .(comp) ~ (lambda == .(signif(eig[comp], 3))))
    else
        paste("Shared component", comp)
}

#' Plot individuals on the shared components
#'
#' @description Scatter plot of the shared canonical components \code{Y} (one
#'   point per individual), optionally coloured by a grouping factor. A modern,
#'   colourblind-safe replacement for \code{\link{plotInds}} that returns a
#'   \code{ggplot} object you can further customise.
#'
#' @param x an \code{mgcca} object (from \code{\link{mgcca_results}} or
#'   \code{mgcca(collect = TRUE)}).
#' @param group optional factor to colour individuals by. If it carries names
#'   they are matched to the row names of \code{Y}; otherwise it is assumed to be
#'   in the same order as the individuals.
#' @param comps the two components to display. Default \code{c(1, 2)}.
#' @param label logical; label points with individual IDs. Default FALSE.
#' @param ellipse logical; when a \code{group} is given, overlay a 95\% confidence
#'   ellipse per group (groups with fewer than 4 individuals, and the \code{NA}
#'   group, are skipped). Makes the separation legible when points overlap.
#'   Default TRUE.
#' @param level confidence level for the ellipses. Default 0.95.
#' @param point_size point size. Default 2.2.
#' @param title plot title.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotVariables}}, \code{\link{plotAVE}}
#' @export
plotIndividuals <- function(x, group = NULL, comps = c(1, 2), label = FALSE,
                            ellipse = TRUE, level = 0.95, point_size = 2.2,
                            title = "Individuals") {
    .check_mgcca(x)
    Y <- as.matrix(x$Y)
    if (max(comps) > ncol(Y)) stop("'comps' exceed the number of components")
    df <- data.frame(comp1 = Y[, comps[1]], comp2 = Y[, comps[2]],
                     id = rownames(Y), stringsAsFactors = FALSE)
    if (!is.null(group)) {
        g <- group
        if (!is.null(names(g)) && !is.null(rownames(Y)))
            g <- g[rownames(Y)]
        df$group <- factor(g)
    }
    eig <- x$AVE$AVE_inner_model
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$comp1, y = .data$comp2)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70")
    if (is.null(group)) {
        p <- p + ggplot2::geom_point(size = point_size, colour = .mgcca_palette[1],
                                     alpha = 0.85)
    } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$group),
                                     size = point_size, alpha = 0.75)
        if (ellipse) {
            # only non-NA groups with enough points support an ellipse
            dfe <- df[!is.na(df$group), , drop = FALSE]
            keep <- names(which(table(droplevels(dfe$group)) >= 4L))
            dfe <- dfe[as.character(dfe$group) %in% keep, , drop = FALSE]
            if (nrow(dfe) > 0) {
                dfe$group <- droplevels(dfe$group)
                p <- p + ggplot2::stat_ellipse(
                    data = dfe, ggplot2::aes(colour = .data$group),
                    type = "t", level = level, linewidth = 0.8)
            }
        }
        p <- p + ggplot2::scale_colour_manual(values = .mgcca_palette, name = NULL,
                                              na.value = "grey75")
    }
    if (label)
        p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$id),
                                    size = 2.6, vjust = -0.6, colour = "grey30")
    p + ggplot2::labs(title = title, x = .comp_label(comps[1], eig),
                      y = .comp_label(comps[2], eig)) +
        .mgcca_theme()
}

#' Correlation circle of variables for one table
#'
#' @description For a chosen table, draws the correlation of each variable with
#'   the two shared components as a correlation circle (arrows within the unit
#'   circle). Only the \code{top} variables with the largest correlation radius
#'   are drawn as labelled arrows, so the plot stays readable even for
#'   omics-scale tables; the rest are shown as faint points.
#'
#' @param x an \code{mgcca} object.
#' @param table table to display: a name or an index into \code{x$corsY}.
#'   Default 1.
#' @param comps the two components to display. Default \code{c(1, 2)}.
#' @param top number of top variables (by correlation radius) to draw as
#'   labelled arrows. Default 10.
#' @param title plot title. Default: the table name.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotIndividuals}}, \code{\link{getSignif}}
#' @export
plotVariables <- function(x, table = 1, comps = c(1, 2), top = 10,
                          title = NULL) {
    .check_mgcca(x)
    nm <- names(x$corsY)
    tb <- if (is.character(table)) table else nm[table]
    if (is.na(tb) || is.null(x$corsY[[tb]]))
        stop("table '", table, "' not found in x$corsY")
    C <- as.matrix(x$corsY[[tb]])
    if (max(comps) > ncol(C)) stop("'comps' exceed the number of components")
    df <- data.frame(c1 = C[, comps[1]], c2 = C[, comps[2]],
                     var = rownames(C), stringsAsFactors = FALSE)
    df$r <- sqrt(df$c1^2 + df$c2^2)
    df <- df[order(-df$r), ]
    top <- min(top, nrow(df))
    lead <- utils::head(df, top)

    circ <- data.frame(t = seq(0, 2 * pi, length.out = 200))
    circ$x <- cos(circ$t); circ$y <- sin(circ$t)
    eig <- x$AVE$AVE_inner_model
    if (is.null(title)) title <- tb

    ggplot2::ggplot() +
        ggplot2::geom_path(data = circ, ggplot2::aes(x = .data$x, y = .data$y),
                           colour = "grey80", linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey85") +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey85") +
        ggplot2::geom_point(data = df, ggplot2::aes(x = .data$c1, y = .data$c2),
                            colour = "grey75", alpha = 0.5, size = 0.8) +
        ggplot2::geom_segment(data = lead,
            ggplot2::aes(x = 0, y = 0, xend = .data$c1, yend = .data$c2),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.16, "cm")),
            colour = .mgcca_palette[1], linewidth = 0.5) +
        ggplot2::geom_text(data = lead,
            ggplot2::aes(x = .data$c1, y = .data$c2, label = .data$var),
            size = 2.8, colour = "grey15", vjust = -0.4) +
        ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
        ggplot2::labs(title = title, subtitle = sprintf("top %d variables", top),
                      x = .comp_label(comps[1], eig),
                      y = .comp_label(comps[2], eig)) +
        .mgcca_theme()
}

#' Barplot of the Average Variance Explained (AVE) per table
#'
#' @description The outer-model AVE (\code{AVE_X}) measures how much of each
#'   table's variance is captured by each shared component. This draws it as a
#'   grouped bar chart (component on the x axis, one bar per table).
#'
#' @param x an \code{mgcca} object.
#' @param title plot title. Default \code{"Average Variance Explained"}.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotScree}}
#' @export
plotAVE <- function(x, title = "Average Variance Explained") {
    .check_mgcca(x)
    A <- as.matrix(x$AVE$AVE_X)                 # components x tables
    tbl <- names(x$corsY)
    if (is.null(tbl) || length(tbl) != ncol(A)) tbl <- paste0("table", seq_len(ncol(A)))
    df <- data.frame(
        component = factor(rep(paste0("comp", seq_len(nrow(A))), times = ncol(A)),
                           levels = paste0("comp", seq_len(nrow(A)))),
        table     = factor(rep(tbl, each = nrow(A)), levels = tbl),
        ave       = as.vector(A))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$component, y = .data$ave,
                                     fill = .data$table)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                          width = 0.7) +
        ggplot2::scale_fill_manual(values = .mgcca_palette, name = NULL) +
        ggplot2::labs(title = title, x = NULL, y = "AVE (outer model)") +
        .mgcca_theme()
}

#' Scree plot of the shared-component eigenvalues
#'
#' @description Bar plot of the eigenvalues (inner-model AVE) associated with
#'   each shared component.
#'
#' @param x an \code{mgcca} object.
#' @param title plot title. Default \code{"Eigenvalues"}.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotAVE}}
#' @export
plotScree <- function(x, title = "Eigenvalues") {
    .check_mgcca(x)
    ev <- as.numeric(x$AVE$AVE_inner_model)
    df <- data.frame(component = factor(paste0("comp", seq_along(ev)),
                                        levels = paste0("comp", seq_along(ev))),
                     eigenvalue = ev)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$component, y = .data$eigenvalue)) +
        ggplot2::geom_col(width = 0.6, fill = .mgcca_palette[1]) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3g", .data$eigenvalue)),
                           vjust = -0.4, size = 3, colour = "grey30") +
        ggplot2::labs(title = title, x = NULL, y = "Eigenvalue (AVE inner)") +
        .mgcca_theme()
}

#' Barplot of the top variable loadings for one component
#'
#' @description For a chosen table and shared component, shows the variables most
#'   correlated with that component as a horizontal bar chart (signed
#'   correlation, coloured by sign). Complements the correlation circle of
#'   \code{\link{plotVariables}} and scales to omics tables via \code{top}.
#'
#' @param x an \code{mgcca} object.
#' @param table table name or index into \code{x$corsY}. Default 1.
#' @param comp component to display. Default 1.
#' @param top number of top variables (by absolute correlation). Default 15.
#' @param title plot title.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotVariables}}, \code{\link{getSignif}}
#' @export
plotLoadings <- function(x, table = 1, comp = 1, top = 15, title = NULL) {
    .check_mgcca(x)
    nm <- names(x$corsY)
    tb <- if (is.character(table)) table else nm[table]
    if (is.na(tb) || is.null(x$corsY[[tb]]))
        stop("table '", table, "' not found in x$corsY")
    C <- as.matrix(x$corsY[[tb]])
    if (comp > ncol(C)) stop("'comp' exceeds the number of components")
    v <- C[, comp]
    o <- order(abs(v), decreasing = TRUE)[seq_len(min(top, length(v)))]
    o <- rev(o)                                     # largest at the top of the bar chart
    df <- data.frame(var = factor(rownames(C)[o], levels = rownames(C)[o]),
                     cor = v[o],
                     sign = ifelse(v[o] >= 0, "positive", "negative"),
                     stringsAsFactors = FALSE)
    if (is.null(title))
        title <- sprintf("%s -- component %d loadings", tb, comp)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$cor, y = .data$var,
                                     fill = .data$sign)) +
        ggplot2::geom_col(width = 0.75) +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
        ggplot2::scale_fill_manual(
            values = c(positive = .mgcca_palette[1], negative = .mgcca_palette[2]),
            guide = "none") +
        ggplot2::labs(title = title, x = "correlation with component", y = NULL) +
        .mgcca_theme()
}

#' Scatter of one table's canonical scores
#'
#' @description Plots the per-table scores (\code{X_j} projected onto its
#'   canonical weights) for one table on two components -- each table's own view
#'   of the individuals. Requires the fit to have been computed with
#'   \code{scores = TRUE}.
#'
#' @param x an \code{mgcca} object.
#' @param table table name or index into \code{x$scores}. Default 1.
#' @param comps the two components to display. Default \code{c(1, 2)}.
#' @param group optional factor to colour by (matched to individuals as in
#'   \code{\link{plotIndividuals}}).
#' @param ellipse,level passed through to draw per-group confidence ellipses.
#' @param title plot title. Default: the table name.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotIndividuals}}, \code{\link{plotBiplot}}
#' @export
plotScores <- function(x, table = 1, comps = c(1, 2), group = NULL,
                       ellipse = TRUE, level = 0.95, title = NULL) {
    .check_mgcca(x)
    if (is.null(x$scores))
        stop("no scores in 'x'; refit with mgcca(..., scores = TRUE)")
    nm <- names(x$scores)
    tb <- if (is.character(table)) table else nm[table]
    if (is.na(tb) || is.null(x$scores[[tb]]))
        stop("table '", table, "' not found in x$scores")
    S <- as.matrix(x$scores[[tb]])
    # reuse the individuals plotter by wrapping the scores as a pseudo-Y
    proxy <- x; proxy$Y <- S
    if (is.null(title)) title <- sprintf("%s scores", tb)
    p <- plotIndividuals(proxy, group = group, comps = comps, ellipse = ellipse,
                         level = level, title = title)
    p + ggplot2::labs(x = paste("Component", comps[1]),
                      y = paste("Component", comps[2]))
}

#' Biplot: individuals and variable loadings together
#'
#' @description Overlays the individuals (shared components \code{Y}) and, for one
#'   table, the loadings of its top variables as arrows -- a compact view of which
#'   variables pull the individuals where. Loadings are scaled to the span of the
#'   individual cloud for readability.
#'
#' @param x an \code{mgcca} object.
#' @param table table name or index into \code{x$corsY} for the variable arrows.
#'   Default 1.
#' @param comps the two components to display. Default \code{c(1, 2)}.
#' @param top number of variables (by correlation radius) to draw. Default 8.
#' @param group optional factor to colour individuals by.
#' @param title plot title.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotIndividuals}}, \code{\link{plotVariables}}
#' @export
plotBiplot <- function(x, table = 1, comps = c(1, 2), top = 8, group = NULL,
                       title = "Biplot") {
    .check_mgcca(x)
    p <- plotIndividuals(x, group = group, comps = comps, ellipse = FALSE,
                         title = title)
    nm <- names(x$corsY)
    tb <- if (is.character(table)) table else nm[table]
    C  <- as.matrix(x$corsY[[tb]])[, comps, drop = FALSE]
    r  <- sqrt(rowSums(C^2))
    o  <- order(r, decreasing = TRUE)[seq_len(min(top, nrow(C)))]
    Y  <- as.matrix(x$Y)[, comps, drop = FALSE]
    # scale unit-ish correlations to ~80% of the individual cloud span
    span <- 0.8 * max(abs(Y))
    lead <- data.frame(x = C[o, 1] * span, y = C[o, 2] * span,
                       var = rownames(C)[o], stringsAsFactors = FALSE)
    p +
        ggplot2::geom_segment(data = lead,
            ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.16, "cm")),
            colour = "grey25", linewidth = 0.5, inherit.aes = FALSE) +
        ggplot2::geom_text(data = lead,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$var),
            colour = "grey10", size = 2.9, vjust = -0.4, inherit.aes = FALSE)
}
