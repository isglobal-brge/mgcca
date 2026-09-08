# =============================================================================
# Internal algebra for the mgcca reliability layer.
#
# Ported from the sealed research implementation, which contains no study-specific
# code: the maths is generic and was verified to be so before porting (no cohort,
# no omic, no dataset names anywhere in it).
#
# WHY THIS LAYER EXISTS AND WHAT IT IS NOT. The reliability layer builds its own
# PARTICIPANT-SPACE reference from the same HDF5 blocks the fit used. It is a
# second method over the same substrate, not a re-reading of `mgcca()`'s output:
# the quantities it needs -- the weighted operator, its eigenbasis and the
# per-block resolvents -- are not part of an `"mgcca"` result object.
#
# HOW IT SCALES. Everything here is n x n in participant space, where n is the
# number of individuals. The only step whose cost grows with the number of
# variables is the per-block accumulation, and that is done by the sealed C++
# kernel streaming from HDF5, never in R. So a block with hundreds of thousands
# of variables costs no more R memory than a block with thirty.
# =============================================================================

# ---- identifier discipline --------------------------------------------------
# Participant alignment is BY IDENTIFIER, never by position. This is not caution
# for its own sake: silent misalignment produces a plausible number rather than
# an error, and is undetectable from the output. Every entry point that takes a
# participant-indexed object goes through here.
.mgcca_align_ids <- function(x, ids, what, allow_missing = TRUE) {
    nm <- if (is.matrix(x) || is.data.frame(x)) rownames(x) else names(x)
    if (is.null(nm))
        stop(what, " must carry participant identifiers (",
             if (is.matrix(x) || is.data.frame(x)) "rownames" else "names",
             "). Positional matching is not supported: a silent misalignment ",
             "would produce a plausible result rather than a failure.", call. = FALSE)
    if (anyDuplicated(nm))
        stop(what, " has duplicated identifiers (",
             paste(utils::head(unique(nm[duplicated(nm)]), 3), collapse = ", "),
             "): alignment would be ambiguous.", call. = FALSE)
    unknown <- setdiff(nm, ids)
    if (length(unknown))
        stop(what, " contains ", length(unknown), " identifier(s) absent from the fit (",
             paste(utils::head(unknown, 3), collapse = ", "),
             if (length(unknown) > 3) ", ..." else "", ").", call. = FALSE)
    i <- match(ids, nm)
    if (!allow_missing && anyNA(i))
        stop(what, " does not cover every participant in the fit (", sum(is.na(i)),
             " missing).", call. = FALSE)
    i
}

# ---- the participant-space reference fit ------------------------------------
# G_j are the per-block participant Grams (n x n, block-absent rows and columns
# zero); present[[j]] the block-presence masks; lambda the per-block ridge.
#
# D = 1/sqrt(K) with K the number of blocks each participant is present in. That
# weighting is what makes the operator well defined when individuals are missing
# from whole blocks, and it is the reason a participant present in no block at
# all is a hard error rather than an Inf.
.mgcca_rel_fit <- function(Glist, present, lambda, L, ids) {
    J <- length(Glist); n <- length(ids)
    if (J < 1L) stop("no blocks supplied", call. = FALSE)
    if (length(present) != J || length(lambda) != J)
        stop("Glist, present and lambda must have the same length", call. = FALSE)
    if (anyDuplicated(ids)) stop("duplicated participant identifiers in the fit", call. = FALSE)
    if (!all(vapply(Glist, function(G) nrow(G) == n && ncol(G) == n, logical(1))))
        stop("every block Gram must be n x n over the participant universe", call. = FALSE)
    # A Gram that carries dimnames must agree with the universe IN ORDER; one that
    # carries none cannot be checked and is taken on the caller's contract.
    if (any(vapply(Glist, function(G) {
            rn <- rownames(G); cn <- colnames(G)
            (!is.null(rn) && !identical(rn, ids)) || (!is.null(cn) && !identical(cn, ids))
        }, logical(1))))
        stop("a block Gram's dimnames are not aligned to the participant universe", call. = FALSE)
    if (!all(vapply(Glist, function(G) all(is.finite(G)), logical(1))))
        stop("non-finite entries in a block Gram", call. = FALSE)

    Kcount <- Reduce(`+`, lapply(present, as.numeric))
    if (any(Kcount < 1))
        stop(sum(Kcount < 1), " participant(s) are present in no block at all, so the ",
             "availability weight 1/sqrt(K) is not defined for them. Remove them from ",
             "the universe or check the presence masks.", call. = FALSE)
    D <- 1 / sqrt(Kcount)

    Rlist <- lapply(seq_len(J), function(j) solve(Glist[[j]] + lambda[j] * diag(n)))
    M <- Reduce(`+`, lapply(seq_len(J), function(j) Glist[[j]] %*% Rlist[[j]]))
    S <- outer(D, D) * M
    S <- (S + t(S)) / 2                      # symmetrise: M is only symmetric up to rounding
    eig <- eigen(S, symmetric = TRUE)
    if (L < 1L || L >= n) stop("L must satisfy 1 <= L < n", call. = FALSE)

    structure(list(
        ids = ids, n = n, J = J, L = L, lambda = lambda, present = present,
        D = D, Glist = Glist, Rlist = Rlist, S = S,
        V = eig$vectors, mu = eig$values,
        gap = eig$values[L] - eig$values[L + 1L]),
        class = "mgcca_rel_fit")
}

# ---- a query in the fit's metric --------------------------------------------
# Participants absent from the query contribute zero, which is the same convention
# the estimator uses for block-absent participants: an absent value is not a value
# of zero, it is an absence, and the metric already accounts for it through D.
.mgcca_rel_zstar <- function(fit, z) {
    i <- .mgcca_align_ids(z, fit$ids, "query")
    zc <- stats::setNames(numeric(fit$n), fit$ids)
    ok <- !is.na(i)
    zv <- as.numeric(z)[i[ok]]
    fin <- is.finite(zv)
    zc[fit$ids[ok][fin]] <- zv[fin]
    if (all(zc == 0))
        stop("the query has no finite value on any participant of the fit", call. = FALSE)
    list(zstar = zc / fit$D, n_used = sum(fin))
}

# ---- the alignment functional T ---------------------------------------------
# T = ||V_L' z*||^2 / ||z*||^2 : the share of the query that lies in the retained
# subspace. Bounded in [0, 1] by construction.
.mgcca_rel_T <- function(fit, z) {
    zs <- .mgcca_rel_zstar(fit, z)$zstar
    w <- crossprod(fit$V[, seq_len(fit$L), drop = FALSE], zs)
    as.numeric(sum(w^2) / sum(zs^2))
}

# ---- the first-order perturbation coefficients ------------------------------
# C is the retained-by-discarded block of the eigenbasis response to a perturbation
# in the direction A, divided by the eigenvalue gaps. It is shared across blocks;
# the per-block parts are Rr_j and Sm_j. This is the standard first-order
# expansion of an eigenprojector, applied -- not derived -- here.
#
# The eigenvalue gap in the denominator is why a fit whose L-th and (L+1)-th
# eigenvalues are nearly equal has an ill-conditioned sensitivity: the coefficients
# blow up because the subspace itself is not well separated. That is reported as a
# validity flag rather than silently returned as a large number.
.mgcca_rel_C <- function(fit, A) {
    L <- fit$L; n <- fit$n
    idxL <- seq_len(L); idxR <- (L + 1L):n
    VtAV <- crossprod(fit$V, A %*% fit$V)
    gaps <- outer(fit$mu[idxL], fit$mu[idxR], `-`)
    2 * VtAV[idxL, idxR, drop = FALSE] / gaps
}

# Per-block inputs for the streaming kernel: Rr_j = R_j D V_L, Sm_j = R_j D V_R.
.mgcca_rel_block_inputs <- function(fit, j) {
    DV <- fit$D * fit$V
    idxL <- seq_len(fit$L); idxR <- (fit$L + 1L):fit$n
    list(Rr = fit$Rlist[[j]] %*% DV[, idxL, drop = FALSE],
         Sm = fit$Rlist[[j]] %*% DV[, idxR, drop = FALSE],
         alpha = fit$lambda[j])
}

# ---- validity of the first-order expansion ----------------------------------
# The expansion is only meaningful where the retained subspace is separated from
# the discarded one. Reported, never silently absorbed.
.mgcca_rel_validity <- function(fit, tol = 1e-8) {
    rel_gap <- fit$gap / max(abs(fit$mu[1L]), .Machine$double.eps)
    list(gap = fit$gap, rel_gap = rel_gap,
         valid = is.finite(rel_gap) && rel_gap > tol,
         reason = if (!is.finite(rel_gap)) "non-finite eigenvalue gap"
                  else if (rel_gap <= tol) "retained and discarded subspaces are not separated"
                  else NA_character_)
}

# ---- grouping ----------------------------------------------------------------
# A user's grouping is a factor, character or integer vector over participants.
# It is converted to contiguous zero-based integers for the kernel, and the
# original labels are kept for the result: a user must never have to guess which
# of their groups an integer refers to.
.mgcca_rel_group <- function(group, ids, min_size = 2L) {
    i <- .mgcca_align_ids(group, ids, "group")
    g <- as.character(group)[i]
    keep <- !is.na(g)
    if (sum(keep) < 2L)
        stop("the grouping covers fewer than two participants of the fit", call. = FALSE)
    f <- factor(g[keep])
    tb <- table(f)
    small <- names(tb)[tb < min_size]
    if (length(small))
        stop(length(small), " group level(s) have fewer than ", min_size,
             " participants (", paste(utils::head(small, 3), collapse = ", "),
             "). A between-group decomposition is not defined for them.", call. = FALSE)
    list(index = which(keep) - 1L,          # 0-based participant positions
         code = as.integer(f) - 1L,          # 0-based group ids
         levels = levels(f),
         n_used = sum(keep),
         n_dropped = sum(!keep))
}

# ---- the layer's own ridge --------------------------------------------------
# ⚠️ THE RELIABILITY LAYER CHOOSES ITS OWN RIDGE, and neither obvious shortcut is
# right.
#
# ZERO IS ALWAYS WRONG, and for a reason that is easy to state incorrectly. A block
# Gram is n x n in PARTICIPANTS, so its rank is at most min(n, p_j). Centring costs
# one further dimension. Therefore:
#   * a WIDE block (p_j >> n, e.g. methylation) still comes back at rank n - 1;
#   * a NARROW block (p_j < n, e.g. a clinical table of four columns beside it)
#     comes back at rank p_j, which can be far below n.
# Measured: n = 50 with p = 500 gives rank 49; n = 620 with p = 4 gives rank 4.
# So the resolvent (G + lambda I)^-1 does not exist at lambda = 0 in EITHER case,
# and the two cases fail by wildly different margins.
#
# INHERITING THE FIT'S LAMBDA IS ALSO WRONG: it was chosen for a different operator
# on a different scale.
#
# The rule is scale-aware by construction: lambda = gamma * (sum of the positive
# eigenvalues) / (numerical rank), i.e. gamma times the MEAN positive eigenvalue.
# It therefore tracks each block's own magnitude instead of imposing an absolute
# number, which is what makes it safe across blocks of wildly different size.
#
# The guards are part of the rule, not decoration: a Gram that is not symmetric,
# has a non-positive trace, is not positive semidefinite, or whose numerical rank
# is below L, produces a named error rather than a lambda that would hide it.
.mgcca_rel_lambda <- function(Glist, L, gamma = 1,
                              tol_sym = 1e-10, tol_rank = 1e-7, tol_pos = 1e-8) {
    vapply(seq_along(Glist), function(j) {
        G <- Glist[[j]]
        nm <- if (!is.null(names(Glist))) names(Glist)[j] else as.character(j)
        if (max(abs(G - t(G))) > tol_sym * max(1, max(abs(G))))
            stop("block '", nm, "': the Gram is not symmetric", call. = FALSE)
        tr <- sum(diag(G))
        if (!is.finite(tr) || tr <= 0)
            stop("block '", nm, "': the Gram's trace is not positive", call. = FALSE)
        ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
        if (min(ev) < -tol_rank * max(ev))
            stop("block '", nm, "': the Gram is not positive semidefinite (min/max = ",
                 format(min(ev) / max(ev), digits = 3), ")", call. = FALSE)
        pos <- sum(ev > tol_pos * max(ev))
        if (pos < L)
            stop("block '", nm, "': numerical rank ", pos, " is below L = ", L,
                 ", so an L-dimensional subspace is not identified in it", call. = FALSE)
        gamma * sum(pmax(ev, 0)) / pos
    }, numeric(1))
}

# ---- optional persistence ---------------------------------------------------
# Writes a result under RELIABILITY/<name>/ in the fit's own file, with the same
# attribute-manifest convention the package already uses for FINAL_RESULTS.
#
# ⚠️ REFUSES TO OVERWRITE. Two runs with different settings under the same name
# would leave a file whose contents no longer match its manifest, and nothing
# downstream could detect it. The caller picks a new name or deletes the old one
# deliberately.
.mgcca_rel_store <- function(res, file, name, tables) {
    grp <- paste0("RELIABILITY/", name)
    existing <- tryCatch(mgcca_list_group_rcpp(file, grp), error = function(e) character(0))
    BigDataStatMeth::hdf5_close_all()
    if (length(existing))
        stop("'", grp, "' already exists in ", file,
             ". Choose another `store_name`, or remove it deliberately: overwriting ",
             "would leave a manifest that no longer describes the data beside it.",
             call. = FALSE)
    for (nm in names(tables)) {
        d <- tables[[nm]]
        if (is.null(d) || !nrow(d)) next
        num <- vapply(d, is.numeric, logical(1))
        if (any(num))
            BigDataStatMeth::hdf5_create_matrix(
                file, paste0(grp, "/", nm),
                data = as.matrix(d[, num, drop = FALSE]), overwrite = TRUE)
        BigDataStatMeth::hdf5_close_all()
    }
    man <- list(mgcca_version = as.character(utils::packageVersion("mgcca")),
                mgcca_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                analysis = class(res)[1],
                backend = as.character(res$settings$backend),
                L = as.integer(res$settings$L %||% res$reference$L),
                lambda = as.numeric(res$settings$lambda),
                lambda_source = as.character(res$settings$lambda_source),
                tables = names(tables))
    mgcca_write_attrs_rcpp(file, grp, "", man)
    BigDataStatMeth::hdf5_close_all()
    list(group = grp, tables = names(tables))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
