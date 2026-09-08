#' Stability of the shared subspace under resampling
#'
#' @description Refits the shared subspace on resamples of the participants and
#'   measures how much it moves. \code{"subsample"} draws random subsets;
#'   \code{"loco"} leaves out one level of \code{group} at a time -- the generic
#'   form of leave-one-cohort-out, for any grouping the user has.
#'
#' @details
#' Two quantities are returned per resample. \code{overlap} is the agreement
#' between the resampled subspace and the reference one, compared \emph{on the
#' participants both contain}. \code{rank_margin} is the smallest singular value
#' of the reference frame restricted to those shared participants.
#'
#' \strong{The rank margin is not a measure of component separation.} It says
#' whether the comparison is well posed on the people both fits have: if the
#' reference frame becomes near-degenerate once restricted, the overlap is
#' comparing directions that are barely present in the shared set, and a high
#' overlap would mean little. An orthonormal frame has every singular value equal
#' to one, so a value well below one is the warning, not the reassurance.
#'
#'   \strong{Why the ridge cannot be zero, and cannot be the fit's.} A block Gram
#'   is \eqn{n \times n} in \emph{participants}, so its rank is at most
#'   \eqn{\min(n, p_j)}, and centring costs one dimension more. A wide block
#'   (\eqn{p_j \gg n}, say methylation) therefore returns rank \eqn{n-1}; a
#'   narrow one (\eqn{p_j < n}, say a clinical table of four columns sitting
#'   beside it) returns rank \eqn{p_j}, which can be far below \eqn{n}. Measured:
#'   \eqn{n=50} with \eqn{p=500} gives rank 49, and \eqn{n=620} with \eqn{p=4}
#'   gives rank 4. So \eqn{(G_j + \lambda I)^{-1}} does not exist at
#'   \eqn{\lambda = 0} in \emph{either} case -- the two simply fail by very
#'   different margins. Inheriting the fit's own \code{lambda} is no better: it
#'   was chosen for a different operator on a different scale. The default rule is
#'   proportional to each block's own mean positive eigenvalue, which is what lets
#'   one number serve a 485,000-column block and a 4-column one in the same fit.
#'
#' \strong{What moves and what does not.} The blocks' standardisation is that of
#' the full fit and is held fixed; only the participant set changes. So this
#' measures stability of the estimated subspace to the sample, not to the whole
#' preprocessing pipeline.
#'
#' @param x A fitted \code{"mgcca"} object; see \code{\link{mgcca_sensitivity}}
#'   for what it must carry.
#' @param method \code{"subsample"} (default), \code{"loco"}, or \code{"both"}.
#'   \code{"loco"} requires \code{group}; \code{"subsample"} does not, so the
#'   function is usable with no grouping at all.
#' @param group Grouping over participants, required for \code{"loco"}.
#' @param strata Optional grouping to stratify the subsampling by. Deliberately
#'   separate from \code{group}: leaving out study centres while stratifying
#'   subsamples by sex is a reasonable thing to want, and one argument serving
#'   both roles would make it impossible.
#' @param B Number of subsamples.
#' @param fraction Fraction of participants retained in each subsample.
#' @param seed Optional seed. The realised resampling plan is archived in the
#'   result either way, so a run can be reproduced from what it returns.
#' @param lambda,gamma Ridge for the layer's operator; see
#'   \code{\link{mgcca_sensitivity}}.
#' @param backend,block_size,threads As in \code{\link{mgcca_sensitivity}}.
#' @param store If \code{TRUE}, persist under \code{RELIABILITY/}. Default
#'   \code{FALSE}: an analysis function should not write to a user's file
#'   unasked.
#' @param store_name Name of the group under \code{RELIABILITY/} when
#'   \code{store = TRUE}. Default \code{"stability"}. An existing name is
#'   refused rather than overwritten, so a second run with different settings
#'   cannot leave a manifest that no longer describes the data beside it.
#'
#' @return An object of class \code{"mgcca_stability"}: a list with
#'   \item{resamples}{one row per resample, with \code{method}, \code{id},
#'     \code{left_out}, the retained \code{n}, the \code{overlap} with the
#'     reference subspace, the \code{rank_margin}, and \code{ok} / \code{reason}
#'     for resamples whose refit failed -- reported, never dropped.}
#'   \item{summary_metrics}{one row per method: number of resamples, number
#'     failed, minimum and median overlap, and the minimum rank margin.}
#'   \item{reference}{the reference spectrum \code{mu}, its eigenvalue
#'     \code{gap}, and \code{L}.}
#'   \item{validity}{the gap check on the reference fit.}
#'   \item{plan}{the realised resampling plan, so a run reproduces from what it
#'     returns rather than from a seed plus matching RNG settings.}
#'   \item{settings, groups, provenance, call}{the record of how it was
#'     produced. When \code{store = TRUE} a \code{storage} element records the
#'     HDF5 group written to.}
#' @seealso \code{\link{mgcca_sensitivity}} for the query-based companion, and
#'   \code{\link{mgcca_select_lambda}} -- instability is a diagnostic at a
#'   chosen penalty, never the way to choose one.
#' @examples
#' ## A small three-block slice of the shipped cohort. As for
#' ## mgcca_sensitivity(), the fit must be file-backed: the resampling refits
#' ## the operator from the SOURCE blocks in HDF5.
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
#' ## B is tiny here so the example runs quickly; a real run uses the default.
#' ## "both" adds leave-one-level-out on top of the subsampling, so the two
#' ## methods appear side by side in the summary.
#' cd4 <- X$cells[, "CD4T"]
#' grp <- stats::setNames(
#'     ifelse(cd4 > stats::median(cd4), "CD4T high", "CD4T low"),
#'     rownames(X$cells))
#'
#' stab <- mgcca_stability(fit, method = "both", group = grp, B = 3, seed = 1)
#' stab
#'
#' ## Read the overlap next to the rank margin: a high overlap on a low margin
#' ## compares directions barely present in the shared participants.
#' stab$summary_metrics
#' stab$resamples
#' @export
mgcca_stability <- function(x, method = c("subsample", "loco", "both"),
                            group = NULL, strata = NULL, B = 200L,
                            fraction = 0.632, seed = NULL,
                            lambda = NULL, gamma = 1,
                            store_name = "stability",
                            backend = c("auto", "memory", "hdf5"),
                            block_size = 0L, store = FALSE, threads = NULL) {
    cl <- match.call()
    .mgcca_refuse_memory(x)          # in-memory fits have no source blocks on disk
    method <- match.arg(method); backend <- match.arg(backend)
    if (method %in% c("loco", "both") && is.null(group))
        stop("method '", method, "' leaves out one level of `group` at a time, so ",
             "`group` is required. Use method = \"subsample\" if you have no grouping.",
             call. = FALSE)
    if (!(fraction > 0 && fraction < 1))
        stop("`fraction` must lie strictly between 0 and 1", call. = FALSE)

    ctx <- .mgcca_rel_context(x); ids <- ctx$ids; n <- length(ids)
    G  <- if (is.null(group))  NULL else .mgcca_rel_group(group, ids)
    ST <- if (is.null(strata)) NULL else .mgcca_rel_group(strata, ids, min_size = 1L)

    chosen <- .mgcca_rel_backend(backend, ctx)
    blocks <- .mgcca_rel_block_grams(ctx, chosen, block_size, threads)
    lam <- if (is.null(lambda)) .mgcca_rel_lambda(blocks$Glist, ctx$L, gamma)
           else if (length(lambda) == 1L) rep(as.numeric(lambda), length(ctx$datasets))
           else as.numeric(lambda)

    ref <- reliability_reference_fit(blocks$Glist, blocks$present, lam, ctx$L)
    Vref <- ref$V[, seq_len(ctx$L), drop = FALSE]
    rownames(Vref) <- ids

    # The plan is drawn ONCE and archived, so the result carries what produced it
    # rather than a seed that only reproduces it under the same RNG settings.
    if (!is.null(seed)) set.seed(seed)
    rng <- list(kind = RNGkind(), seed = seed)
    plan <- list()
    if (method %in% c("subsample", "both")) {
        k <- max(2L, floor(fraction * n))
        for (b in seq_len(B)) {
            keep <- if (is.null(ST)) sort(sample.int(n, k)) else
                sort(unlist(lapply(split(seq_len(n), ST$code[match(seq_len(n), ST$index + 1L)]),
                                   function(idx) {
                                       idx <- idx[!is.na(idx)]
                                       if (!length(idx)) return(integer(0))
                                       sample(idx, max(1L, floor(fraction * length(idx))))
                                   })))
            plan[[length(plan) + 1L]] <- list(method = "subsample", id = b,
                                              keep = keep, left_out = NA_character_)
        }
    }
    if (method %in% c("loco", "both")) {
        gi <- G$index + 1L
        for (l in seq_along(G$levels)) {
            drop_pos <- gi[G$code == (l - 1L)]
            plan[[length(plan) + 1L]] <- list(method = "loco", id = l,
                                              keep = setdiff(seq_len(n), drop_pos),
                                              left_out = G$levels[l])
        }
    }

    rows <- vector("list", length(plan))
    for (i in seq_along(plan)) {
        pl <- plan[[i]]; keep <- pl$keep
        Gs <- lapply(blocks$Glist,  function(g) g[keep, keep, drop = FALSE])
        Ps <- lapply(blocks$present, function(p) p[keep])
        fitb <- try(reliability_reference_fit(Gs, Ps, lam, ctx$L), silent = TRUE)
        if (inherits(fitb, "try-error")) {
            rows[[i]] <- data.frame(method = pl$method, id = pl$id,
                                    left_out = pl$left_out, n = length(keep),
                                    overlap = NA_real_, rank_margin = NA_real_,
                                    ok = FALSE,
                                    reason = trimws(conditionMessage(attr(fitb, "condition"))),
                                    stringsAsFactors = FALSE)
            next
        }
        Vb <- fitb$V[, seq_len(ctx$L), drop = FALSE]
        Vr <- Vref[keep, , drop = FALSE]
        # Is the comparison well posed on the shared participants at all?
        margin <- min(svd(Vr, nu = 0, nv = 0)$d)
        ov <- reliability_subspace_overlap(Vr, Vb)
        rows[[i]] <- data.frame(method = pl$method, id = pl$id,
                                left_out = pl$left_out, n = length(keep),
                                overlap = ov, rank_margin = margin, ok = TRUE,
                                reason = NA_character_, stringsAsFactors = FALSE)
    }
    res_tab <- do.call(rbind, rows)

    smry <- do.call(rbind, lapply(split(res_tab, res_tab$method), function(d) {
        o <- d$overlap[d$ok]
        data.frame(method = d$method[1], n_resamples = nrow(d), n_failed = sum(!d$ok),
                   overlap_min = if (length(o)) min(o) else NA_real_,
                   overlap_median = if (length(o)) stats::median(o) else NA_real_,
                   rank_margin_min = if (length(o)) min(d$rank_margin[d$ok]) else NA_real_,
                   stringsAsFactors = FALSE)
    }))

    out <- list(resamples = res_tab, summary_metrics = smry,
                reference = list(mu = ref$mu, gap = ref$gap, L = ctx$L),
                validity = .mgcca_rel_validity(list(gap = ref$gap, mu = ref$mu)),
                settings = list(method = method, B = B, fraction = fraction,
                                backend = chosen, block_size = block_size,
                                lambda = lam, gamma = gamma,
                                lambda_source = if (is.null(lambda)) "guarded rule" else "user",
                                strata = if (is.null(ST)) NULL else ST$levels,
                                rng = rng),
                plan = plan, groups = if (is.null(G)) NULL else G$levels,
                provenance = ctx$provenance, call = cl)
    class(out) <- "mgcca_stability"
    if (isTRUE(store))
        out$storage <- .mgcca_rel_store(out, ctx$file, store_name,
                                        list(resamples = out$resamples, summary_metrics = out$summary_metrics))
    out
}

#' @export
print.mgcca_stability <- function(x, ...) {
    s <- x$settings
    cat("mgcca subspace stability\n")
    cat(sprintf("  method   : %s ; L = %d ; backend = %s\n",
                s$method, x$reference$L, s$backend))
    if (s$method %in% c("subsample", "both"))
        cat(sprintf("  subsample: B = %d, fraction = %.3f%s\n", s$B, s$fraction,
                    if (is.null(s$strata)) "" else
                        paste0(", stratified by ", length(s$strata), " levels")))
    if (!is.null(x$groups))
        cat(sprintf("  groups   : %s\n", paste(x$groups, collapse = ", ")))
    cat("\n"); print(format(x$summary_metrics, digits = 4), row.names = FALSE)
    nf <- sum(!x$resamples$ok)
    if (nf)
        cat(sprintf("\n  ! %d resample(s) did not produce a fit; see $resamples$reason.\n", nf))
    cat("\n  overlap = agreement with the reference subspace on the SHARED participants.\n")
    cat("  rank_margin = smallest singular value of the reference frame restricted to\n")
    cat("  them: it says whether that comparison is well posed, NOT whether the\n")
    cat("  components are separated. An orthonormal frame has all singular values 1,\n")
    cat("  so a value well below 1 is a warning about the comparison itself.\n")
    invisible(x)
}

#' @export
summary.mgcca_stability <- function(object, ...) {
    print(object)
    cat("\nPer-resample\n")
    d <- object$resamples
    print(format(utils::head(d[order(d$overlap), ], 10), digits = 4), row.names = FALSE)
    if (nrow(d) > 10) cat("  ... ", nrow(d) - 10, " more (worst shown first)\n", sep = "")
    invisible(object)
}

#' @export
plot.mgcca_stability <- function(x, ...) {
    d <- x$resamples[x$resamples$ok, , drop = FALSE]
    if (!nrow(d)) stop("no resample produced a fit, so there is nothing to plot",
                       call. = FALSE)
    d$label <- ifelse(is.na(d$left_out), paste0("#", d$id), d$left_out)
    ggplot2::ggplot(d, ggplot2::aes(x = .data$rank_margin, y = .data$overlap,
                                    colour = .data$method)) +
        ggplot2::geom_point(size = 2, alpha = 0.8) +
        ggplot2::scale_colour_manual(values = .mgcca_palette, name = NULL) +
        ggplot2::labs(x = "rank margin (is the comparison well posed?)",
                      y = "subspace overlap with the reference",
                      title = "Subspace stability under resampling",
                      subtitle = "points at a low rank margin compare directions barely present in the shared set") +
        .mgcca_theme()
}
