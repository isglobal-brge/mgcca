#' Sensitivity of the shared subspace to a participant-level query
#'
#' @description Measures how much of a participant-level quantity -- a phenotype,
#'   an exposure, a score -- aligns with the shared subspace estimated by
#'   \code{\link{mgcca}}, and how much of that alignment's sensitivity is carried
#'   by a grouping variable such as batch, centre, plate or cohort.
#'
#' @details
#' Two quantities are returned for each query. \code{T} is the share of the query
#' that lies in the retained subspace, bounded in \eqn{[0,1]}. \code{S_total} is
#' the first-order sensitivity of that alignment to a structured perturbation of
#' the blocks; when \code{group} is supplied, \code{S_between} is the part of that
#' sensitivity attributable to differences between groups, and \code{R_between}
#' their ratio.
#'
#' \strong{Read the ratio with its components, never alone.} A high
#' \code{R_between} can accompany a negligible absolute sensitivity, and a low one
#' can accompany a large signal; the ratio says how the sensitivity is distributed,
#' not how much there is.
#'
#' \strong{Participants are matched by identifier, never by position.} A query
#' without names is an error rather than an assumption, because a silent
#' misalignment produces a plausible number instead of a failure and cannot be
#' detected from the output.
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
#' \strong{What this is not.} These are first-order sensitivity coefficients of a
#' fitted operator. They are not causal, and they are not a psychometric
#' reliability coefficient.
#'
#' @param x A fitted object of class \code{"mgcca"}. It must carry a descriptor
#'   pointing at an HDF5 file that still contains the input blocks: the source
#'   data is required, and a results-only object is refused rather than
#'   half-processed. Use \code{\link{mgcca_load}} to obtain one from a file.
#' @param query Numeric vector or matrix over participants. A vector must be
#'   named; a matrix must have row names, and column names that are unique
#'   because they become the query identifiers. Participants absent from the
#'   query contribute nothing and are counted in \code{n_used}.
#' @param group Optional grouping over participants: factor, character or
#'   integer. Supplying it adds the between-group decomposition; omitting it
#'   returns the total quantities only, with the group-dependent fields absent
#'   rather than filled with zeroes. Levels with fewer than two participants are
#'   refused, since a between-group quantity is not defined for them.
#' @param covariates Optional named, participant-aligned numeric matrix or data
#'   frame. The query is residualised on these before the sensitivity is taken.
#'   The realised design, its rank and the retained participants are archived in
#'   the result.
#' @param lambda Ridge for the layer's own operator: \code{NULL} (default) uses a
#'   scale-aware guarded rule, \code{gamma} times each block's mean positive
#'   eigenvalue. Supply one value, or one per block, to override it. The fit's own
#'   \code{lambda} is deliberately NOT inherited: it was chosen for a different
#'   operator on a different scale.
#' @param gamma Multiplier for the guarded ridge rule. Ignored when
#'   \code{lambda} is supplied.
#' @param backend One of \code{"auto"} (default), \code{"memory"} or
#'   \code{"hdf5"}. \code{"auto"} chooses by block size, as
#'   \code{mgcca(route = "auto")} does. The other two force the path; the choice
#'   is recorded in the result so a run is never ambiguous about how it was
#'   produced.
#' @param block_size Block size for the out-of-core steps. \code{0} lets
#'   BigDataStatMeth choose. A positive value forces it, which exists so the
#'   blocked path can be exercised on small data.
#' @param store If \code{TRUE}, persist the result into the fit's HDF5 file under
#'   \code{RELIABILITY/}. Default \code{FALSE}: an analysis function should not
#'   write to a user's file unasked.
#' @param store_name Name of the group under \code{RELIABILITY/} when
#'   \code{store = TRUE}. An existing name is refused rather than overwritten.
#' @param threads Optional thread count.
#'
#' @return An object of class \code{"mgcca_sensitivity"}: a list with
#'   \code{overall} (a data frame, one row per query, carrying \code{n_used},
#'   \code{T}, \code{S_total} and -- only when \code{group} is supplied --
#'   \code{S_between} and \code{R_between}), \code{by_block} (one row per query
#'   and block), \code{by_individual}, \code{by_individual_check},
#'   \code{validity} (the eigenvalue-gap check on the first-order expansion),
#'   \code{reference} (the layer's own participant coordinates \code{Y}, the
#'   query \code{directions}, and the grouping actually used),
#'   \code{groups} (the group levels, or \code{NULL}), \code{settings},
#'   \code{provenance} and the matched \code{call}. When \code{store = TRUE} a
#'   \code{storage} element records the HDF5 group the result was written to.
#'
#'   \code{by_individual} decomposes each block's sensitivity over participants,
#'   in long form with columns \code{query}, \code{id}, \code{block},
#'   \code{S_total}, \code{share_total} and -- only when \code{group} is
#'   supplied -- \code{S_between} and \code{share_between}. Writing \eqn{s_{ij}}
#'   for a participant's contribution and \eqn{S_j} for the block quantity
#'   already reported in \code{by_block}, the decomposition is
#'   \eqn{\sum_i s_{ij} = S_j}.
#'
#'   \strong{Rows.} \code{by_individual} contains exactly the ordered
#'   participant set used by the grouped sensitivity calculation, including
#'   participants numerically absent from a particular block. Key it by the
#'   \code{query}, \code{id} and \code{block} columns; row order is not part of
#'   the contract. A grouping that covers fewer participants therefore yields
#'   different contributions: that is the estimand being conditioned on the
#'   grouped set, not an instability.
#'
#'   \strong{Shares.} \code{share_total} is
#'   \eqn{s_{ij} / S_j}, so the shares of a given query and block sum to one;
#'   \code{share_between} is defined the same way over the between-group parts.
#'   Note that \code{S_between} is built from group means, so it is constant
#'   across the members of a group: it is the group's between-group contribution
#'   carried in each member's row, and ranking participants by it ranks their
#'   groups. \code{S_total} is the per-participant quantity.
#'   When a denominator is not finite and strictly positive the share is
#'   \code{NA}, never a silent zero. The absolute contributions inherit the
#'   non-identified scale of the parent sensitivity quantity; the shares remove
#'   this common multiplicative scale and are the preferred dimensionless
#'   representation, but they remain conditional on the fitted model, the ridge
#'   specification, the query, the block and the grouped participant set. No
#'   invariant ranking of participants is implied. Like \eqn{S_j}, the
#'   contributions are invariant to rescaling the query or a block and are NOT
#'   invariant to the layer's ridge. A participant absent from a block
#'   contributes a numerical zero there (of order 1e-34), not an exact zero,
#'   which is why a per-entry relative comparison is the wrong way to check them.
#'
#'   \code{by_individual_check} reports the decomposition identity actually
#'   achieved: \code{sum_rel_total}, and \code{sum_rel_between} when grouped, are
#'   the largest \eqn{|\sum_i s_{ij} - S_j| / \max(1, |S_j|)} over queries and
#'   blocks. The identity is definitional but not bit for bit, because a trace
#'   and a sum over rows accumulate in different orders.
#'
#'   \strong{What these values are.} \code{by_individual} is an exploratory
#'   decomposition of the first-order sensitivity diagnostic. Its block-level
#'   parent has empirical calibration evidence with a MIXTO verdict; no
#'   individual-level calibration has been established. The reported values are
#'   sensitivity contributions, not variance, uncertainty, causal influence, or
#'   reliability guarantees.
#'
#' @seealso \code{\link{mgcca_block_profile}} to reformat the per-block
#'   decomposition, \code{\link{mgcca_stability}} for the resampling companion,
#'   and \code{\link{mgcca}} for the fit itself.
#' @examples
#' ## A small three-block slice of the shipped cohort. The reliability layer
#' ## reads the SOURCE blocks back from HDF5, so the fit must be file-backed:
#' ## an in-memory fit (filename = NULL) is refused, not silently approximated.
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
#' ## The query is BMI; the grouping splits the cohort on CD4T cell proportion.
#' ## A one-column matrix is used so the query carries its own name into the
#' ## result; a named vector would work too.
#' bmi <- X$clinical[, "BMI", drop = FALSE]
#' cd4 <- X$cells[, "CD4T"]
#' grp <- stats::setNames(
#'     ifelse(cd4 > stats::median(cd4), "CD4T high", "CD4T low"),
#'     rownames(X$cells))
#'
#' sens <- mgcca_sensitivity(fit, query = bmi, group = grp)
#' sens
#'
#' ## Read the ratio with its components, never alone.
#' sens$overall
#' sens$by_block
#' @export
mgcca_sensitivity <- function(x, query, group = NULL, covariates = NULL,
                              lambda = NULL, gamma = 1,
                              store_name = "sensitivity",
                              backend = c("auto", "memory", "hdf5"),
                              block_size = 0L, store = FALSE, threads = NULL) {
    cl <- match.call()
    .mgcca_refuse_memory(x)          # in-memory fits have no source blocks on disk
    backend <- match.arg(backend)
    ctx <- .mgcca_rel_context(x)
    ids <- ctx$ids

    Q <- .mgcca_rel_query_matrix(query, ids)
    G <- if (is.null(group)) NULL else .mgcca_rel_group(group, ids)
    CV <- if (is.null(covariates)) NULL else .mgcca_rel_covariates(covariates, ids)

    chosen <- .mgcca_rel_backend(backend, ctx)
    blocks <- .mgcca_rel_block_grams(ctx, chosen, block_size, threads)

    # The ridge is the LAYER's, not the fit's: see .mgcca_rel_lambda for why
    # neither zero nor the fit's lambda is the right default.
    lam <- if (is.null(lambda)) .mgcca_rel_lambda(blocks$Glist, ctx$L, gamma)
           else { if (length(lambda) == 1L) rep(as.numeric(lambda), length(ctx$datasets))
                  else as.numeric(lambda) }
    if (length(lam) != length(ctx$datasets))
        stop("`lambda` must be one value, or one per block (", length(ctx$datasets), ")",
             call. = FALSE)
    fit <- reliability_reference_fit(blocks$Glist, blocks$present, lam, ctx$L)
    val <- .mgcca_rel_validity(list(gap = fit$gap, mu = fit$mu))

    D <- fit$D
    n <- length(ids)
    # Everything the sealed kernel needs that does not depend on the query.
    prm  <- if (is.null(G)) which(rep(TRUE, n)) - 1L else G$index
    grpc <- if (is.null(G)) rep(0L, n) else G$code
    ginfo <- if (is.null(G)) list(levels = NULL, n_used = n, n_dropped = 0L) else G

    overall <- vector("list", ncol(Q))
    byblock <- vector("list", ncol(Q))
    byind   <- vector("list", ncol(Q))
    dirs    <- vector("list", ncol(Q))
    # The decomposition identity holds definitionally but NOT bit for bit: the
    # kernel's trace and a sum over rows accumulate in different orders. It is
    # therefore checked with a tolerance rather than with `identical()`, and the
    # worst departure is reported so a caller can gate on it.
    ind_rel_tot <- 0
    ind_rel_btw <- 0
    for (k in seq_len(ncol(Q))) {
        z <- Q[, k]
        used <- is.finite(z) & z != 0
        if (!is.null(CV)) z <- .mgcca_rel_residualise(z, CV$design)
        zs <- z / D
        # ⚠️ A NUMERICALLY zero query is not an exactly zero one. Residualising a
        # query on a design that already spans it leaves a residual of order 1e-16
        # times the original scale -- not `== 0`, but carrying no information. An
        # exact test lets that through and returns a meaningless T computed from
        # rounding error. The comparison is therefore RELATIVE to the query's own
        # magnitude before residualisation.
        scale0 <- max(abs(Q[, k] / D))
        if (max(abs(zs)) <= 1e-10 * max(scale0, .Machine$double.eps))
            stop("query '", colnames(Q)[k], "' has no variation left after ",
                 "residualisation on the covariates: it is spanned by them, so its ",
                 "alignment is not defined.", call. = FALSE)
        qk <- reliability_query(fit$V, fit$mu, zs, ctx$L)

        # --- per-block sensitivity through the sealed kernel ---
        # It reads each block's PRESENT-ONLY Gram, which is why the Gram is kept
        # rather than only its embedded copy: the kernel indexes participants in
        # the block's own order and maps them to fit rows through `mids_to_fit`.
        dirs[[k]] <- qk$w
        bb <- vector("list", length(ctx$datasets))
        ii <- vector("list", length(ctx$datasets))
        for (j in seq_along(ctx$datasets)) {
            bi <- reliability_block_inputs(fit$Rlist[[j]], D, fit$V, ctx$L)
            m2f <- match(blocks$block_ids[[j]], ids) - 1L
            kk <- reliability_sensitivity_gram(
                blocks$gram_file, blocks$gram_group, ctx$datasets[j],
                as.integer(m2f), bi$Rr, bi$Sm, qk$C, lam[j],
                blocks$present[[j]], as.integer(prm), as.integer(grpc))
            bb[[j]] <- data.frame(query = colnames(Q)[k], block = ctx$datasets[j],
                                  S_total = kk$trtH, S_between = kk$trbH,
                                  stringsAsFactors = FALSE)
            # The kernel's accumulators are already per participant: one row per
            # member of the GROUPED set, in `prm` order. The rows are squared and
            # summed here rather than in C++ so that src/ stays untouched.
            ii[[j]] <- data.frame(query = colnames(Q)[k], id = ids[prm + 1L],
                                  block = ctx$datasets[j],
                                  S_total = rowSums(kk$accX^2),
                                  S_between = rowSums(kk$accB^2),
                                  stringsAsFactors = FALSE)
        }
        bb <- do.call(rbind, bb)
        ii <- do.call(rbind, ii)
        # ⚠️ THE SHARE DENOMINATOR IS THE BLOCK'S OWN TOTAL, not the whole
        # participant-by-block grid. A share answers "how much of THIS block's
        # sensitivity does this participant carry", so it sums to one within each
        # query x block. Dividing by the grid instead would answer a different
        # question with the same column name, which is the worse failure.
        #
        # The divisor is the published `by_block` value, so the share divides by
        # the number the caller can see rather than by a second sum of its own.
        #
        # A denominator that is not a positive finite number leaves the share
        # UNDEFINED, and undefined is reported as NA -- never as a zero, which
        # would read as a real answer to a question that has none.
        ii$share_total   <- NA_real_
        ii$share_between <- NA_real_
        for (j in seq_along(ctx$datasets)) {
            m <- ii$block == ctx$datasets[j]
            St_j <- bb$S_total[j]
            ind_rel_tot <- max(ind_rel_tot,
                abs(sum(ii$S_total[m]) - St_j) / max(1, abs(St_j)))
            if (is.finite(St_j) && St_j > 0)
                ii$share_total[m] <- ii$S_total[m] / St_j
            if (!is.null(G)) {
                Sb_j <- bb$S_between[j]
                ind_rel_btw <- max(ind_rel_btw,
                    abs(sum(ii$S_between[m]) - Sb_j) / max(1, abs(Sb_j)))
                if (is.finite(Sb_j) && Sb_j > 0)
                    ii$share_between[m] <- ii$S_between[m] / Sb_j
            }
        }
        if (is.null(G)) ii$S_between <- NA_real_
        byind[[k]] <- ii
        # Shares are reported next to the absolutes, never instead of them.
        bb$share_total   <- bb$S_total   / sum(bb$S_total)
        bb$share_between <- if (is.null(G)) NA_real_ else bb$S_between / sum(bb$S_between)
        if (is.null(G)) bb$S_between <- NA_real_
        byblock[[k]] <- bb

        St <- sum(bb$S_total); Sb <- if (is.null(G)) NA_real_ else sum(bb$S_between)
        overall[[k]] <- data.frame(
            query = colnames(Q)[k], n_used = sum(used), T = qk$T,
            S_total = St, S_between = Sb,
            R_between = if (is.null(G)) NA_real_ else Sb / St,
            stringsAsFactors = FALSE)
    }
    overall <- do.call(rbind, overall)
    byblock <- do.call(rbind, byblock)
    byind   <- do.call(rbind, byind)
    rownames(byind) <- NULL
    # ⚠️ THIS DOES NOT GO INSIDE `validity`. That element must stay exactly what
    # the previous version produced, so a regression bridge can compare it with
    # identical(); a new field inside it would break that comparison for every
    # caller at once. The decomposition check gets its own documented element.
    indchk <- list(sum_rel_total = ind_rel_tot)
    if (!is.null(G)) indchk$sum_rel_between <- ind_rel_btw
    # ⚠️ ABSENT MEANS ABSENT. Without a grouping the between-group quantities are
    # not defined, so the columns are DROPPED rather than filled with NA. A column
    # of NAs is a column that should not be there: it invites the reader to wonder
    # what went wrong, when in fact nothing did -- they simply asked a question
    # that has no group in it.
    if (is.null(G)) {
        overall <- overall[, setdiff(names(overall), c("S_between", "R_between")), drop = FALSE]
        byblock <- byblock[, setdiff(names(byblock), c("S_between", "share_between")), drop = FALSE]
        byind   <- byind[, setdiff(names(byind), c("S_between", "share_between")), drop = FALSE]
    }

    # The map view needs two things and neither is worth recomputing later: the
    # coordinates of the individuals on the layer's OWN reference axes, and where
    # each query points on them.
    Yref <- (fit$D * fit$V)[, seq_len(ctx$L), drop = FALSE]
    dimnames(Yref) <- list(ids, paste0("comp", seq_len(ctx$L)))
    W <- do.call(cbind, dirs); dimnames(W) <- list(paste0("comp", seq_len(ctx$L)),
                                                   colnames(Q))

    res <- list(overall = overall, by_block = byblock, by_individual = byind,
                by_individual_check = indchk, validity = val,
                reference = list(Y = Yref, directions = W,
                                 group = if (is.null(G)) NULL else
                                     stats::setNames(G$levels[G$code + 1L],
                                                     ids[G$index + 1L])),
                groups = ginfo$levels,
                settings = list(backend = chosen, block_size = block_size,
                                L = ctx$L, lambda = lam, gamma = gamma,
                                lambda_source = if (is.null(lambda)) "guarded rule" else "user",
                                grouped = !is.null(G), covariates = !is.null(CV)),
                provenance = ctx$provenance, call = cl)
    class(res) <- "mgcca_sensitivity"
    if (isTRUE(store))
        res$storage <- .mgcca_rel_store(res, ctx$file, store_name,
                                        list(overall = res$overall, by_block = res$by_block,
                                             by_individual = res$by_individual))
    res
}

# ---- context: what the layer needs from a fit, or a precise refusal ---------
# The reliability layer needs the SOURCE blocks, not just the results. An
# "mgcca" object carries only final results; its descriptor carries the file and
# the input group. Anything missing is named, because "it did not work" is not
# something a user can act on.
.mgcca_rel_context <- function(x) {
    if (!inherits(x, "mgcca"))
        stop("`x` must be a fitted object of class \"mgcca\". ",
             "If you have a file, load it first with mgcca_load().", call. = FALSE)
    d <- attr(x, "desc")
    if (is.null(d))
        stop("this \"mgcca\" object carries no provenance descriptor, so the input ",
             "blocks cannot be located. Refit, or reload the fit with mgcca_load().",
             call. = FALSE)
    if (is.null(d$filename) || !nzchar(d$filename) || !file.exists(d$filename))
        stop("the fit's HDF5 file is not available (",
             if (is.null(d$filename)) "no filename recorded" else d$filename,
             "). The reliability analyses need the source blocks, not only the results.",
             call. = FALSE)
    ig <- d$input_group
    if (is.null(ig) || is.na(ig) || !nzchar(ig))
        stop("this fit does not record which HDF5 group holds its INPUT blocks, so the ",
             "source data cannot be located. Fits written by older versions of mgcca ",
             "do not carry it. Refit, or supply the input group explicitly.",
             call. = FALSE)
    ds <- as.character(d$datasets)
    if (!length(ds)) stop("the fit records no datasets", call. = FALSE)

    ids <- rownames(x$Y)
    if (is.null(ids) || anyDuplicated(ids))
        stop("the fit's participant identifiers are missing or duplicated", call. = FALSE)
    L <- if (!is.null(d$nfac)) as.integer(d$nfac) else ncol(x$Y)
    lam <- if (!is.null(d$lambda) && length(d$lambda) == length(ds))
        as.numeric(d$lambda) else rep(0, length(ds))

    list(file = d$filename, input_group = ig, datasets = ds, ids = ids,
         L = L, lambda = lam,
         provenance = list(file = d$filename, input_group = ig, datasets = ds,
                           mgcca_version = d$mgcca_version, method = d$method,
                           route = d$route, scale = d$scale))
}

# ---- backend choice ---------------------------------------------------------
# "auto" mirrors mgcca(route = "auto"): decide by size, but let the caller force
# either path. Automatic behaviour with no override turns into a support problem
# the user cannot escape; one argument prevents a whole class of them.
.mgcca_rel_backend <- function(backend, ctx, max_cells = 5e7) {
    if (backend != "auto") return(backend)
    cells <- 0
    for (ds in ctx$datasets) {
        h <- try(BigDataStatMeth::hdf5_matrix(ctx$file,
                     paste0(ctx$input_group, "/", ds)), silent = TRUE)
        if (!inherits(h, "try-error")) {
            cells <- max(cells, prod(dim(h)))
            try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
        }
    }
    if (cells > max_cells) "hdf5" else "memory"
}

# ---- per-block participant Grams -------------------------------------------
# The stored input block holds exactly its own individuals, so standardising it
# IS the present-only convention. The padded copies in the temporary group must
# never be read here: doing so would reintroduce the zero-padding defect that
# this layer exists to prevent, and would do it silently.
.mgcca_rel_block_grams <- function(ctx, backend, block_size, threads) {
    n <- length(ctx$ids)
    Glist <- vector("list", length(ctx$datasets))
    present <- vector("list", length(ctx$datasets))
    block_ids <- vector("list", length(ctx$datasets))
    # The sealed kernel reads each block's PRESENT-ONLY Gram from HDF5, so it is
    # written there on both backends. It is n_pr x n_pr in PARTICIPANTS -- small
    # whatever the block's width -- so this costs nothing even for the memory path.
    gram_file <- ctx$file; gram_group <- "RELIABILITY_TMP/G" 
    for (j in seq_along(ctx$datasets)) {
        ds <- ctx$datasets[j]
        path <- paste0(ctx$input_group, "/", ds)
        h <- BigDataStatMeth::hdf5_matrix(ctx$file, path)
        bids <- rownames(h)
        if (is.null(bids))
            stop("block '", ds, "' has no participant row names, so its individuals ",
                 "cannot be identified", call. = FALSE)
        unknown <- setdiff(bids, ctx$ids)
        if (length(unknown))
            stop("block '", ds, "' contains ", length(unknown),
                 " participant(s) absent from the fit", call. = FALSE)
        Gsub <- if (identical(backend, "memory")) {
            Xj <- as.matrix(h)
            g <- tcrossprod(scale(Xj, center = TRUE, scale = TRUE))
            dimnames(g) <- list(bids, bids)
            BigDataStatMeth::hdf5_close_all()
            BigDataStatMeth::hdf5_create_matrix(ctx$file,
                paste0(gram_group, "/", ds), data = g, overwrite = TRUE)
            g
        } else {
            BigDataStatMeth::hdf5_close_all()
            r <- reliability_gram_block(ctx$file, ctx$input_group, ds,
                                        gram_group, as.integer(block_size), threads)
            as.matrix(BigDataStatMeth::hdf5_matrix(ctx$file, r$path))
        }
        BigDataStatMeth::hdf5_close_all()
        block_ids[[j]] <- bids
        Gfull <- matrix(0, n, n, dimnames = list(ctx$ids, ctx$ids))
        i <- match(bids, ctx$ids)
        Gfull[i, i] <- Gsub
        Glist[[j]] <- Gfull
        present[[j]] <- ctx$ids %in% bids
    }
    list(Glist = Glist, present = present, block_ids = block_ids,
         gram_file = gram_file, gram_group = gram_group)
}

# ---- query, group and covariate contracts ----------------------------------
.mgcca_rel_query_matrix <- function(query, ids) {
    if (is.null(dim(query))) {
        i <- .mgcca_align_ids(query, ids, "query")
        Q <- matrix(0, length(ids), 1L, dimnames = list(ids, "query"))
        ok <- !is.na(i)
        v <- as.numeric(query)[i[ok]]
        fin <- is.finite(v)
        Q[ids[ok][fin], 1L] <- v[fin]
        return(Q)
    }
    Q0 <- as.matrix(query)
    if (is.null(colnames(Q0)))
        stop("a query matrix must have column names: they become the query identifiers.",
             call. = FALSE)
    if (anyDuplicated(colnames(Q0)))
        stop("the query matrix has duplicated column names, so its results could not ",
             "be told apart.", call. = FALSE)
    i <- .mgcca_align_ids(Q0, ids, "query")
    Q <- matrix(0, length(ids), ncol(Q0), dimnames = list(ids, colnames(Q0)))
    ok <- !is.na(i)
    sub <- Q0[i[ok], , drop = FALSE]
    sub[!is.finite(sub)] <- 0
    Q[ids[ok], ] <- sub
    Q
}

.mgcca_rel_covariates <- function(covariates, ids) {
    M <- if (is.data.frame(covariates)) as.matrix(covariates) else as.matrix(covariates)
    if (!is.numeric(M))
        stop("`covariates` must be numeric in this version. Expand factors yourself, ",
             "so the realised design is the one you intended rather than one this ",
             "function guessed.", call. = FALSE)
    i <- .mgcca_align_ids(M, ids, "covariates", allow_missing = FALSE)
    D <- cbind(`(Intercept)` = 1, M[i, , drop = FALSE])
    q <- qr(D)
    list(design = D, rank = q$rank, ncol = ncol(D), ids = ids)
}

.mgcca_rel_residualise <- function(z, design) {
    q <- qr(design)
    as.numeric(z - qr.fitted(q, z))
}
