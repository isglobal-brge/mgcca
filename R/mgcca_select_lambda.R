## =============================================================================
## Choosing lambda: the held-out cross-block agreement criterion.
##
## PROVENANCE. This file is a port of the selection criterion developed and
## validated in the thesis reliability calibration:
##
##   reliability/05_lambda_selection.md          -- the criterion and its result
##   reliability/calibration/estimator.R:33-77   -- fit_pgcca (the estimator)
##   reliability/calibration/estimator.R:96-124  -- lambda_scale_matched (s_j)
##   reliability/calibration/estimator.R:445-476 -- cv_crossblock (the criterion)
##   reliability/calibration/estimator.R:501-671 -- cv_crossblock_gamma (the
##                                                  frozen-contract version over
##                                                  the dimensionless ridge)
##   MGCCA_Simulations/R/metrics.R:23-30         -- subspace_rv
##
## WHY IT EXISTS, IN ONE PARAGRAPH. A penalized fit needs a lambda, and the
## tempting way to pick one -- take the lambda whose fit is most stable under
## resampling -- is wrong, and not mildly. Instability is minimised exactly where
## the operator has been regularised into inertia: each block latches onto its
## own high-variance, NON-shared axis, the fit stops moving because there is
## nothing left in it that could move, and the shared structure GCCA exists to
## find is gone. In the specific-dominant regime of 05_lambda_selection.md, the
## stability-optimal lambda costs 0.666 of subspace recovery; the criterion
## implemented here costs 0.002. That is the whole argument for this function.
##
## THE CRITERION. K-fold over individuals. Fit on the training folds; project
## each held-out individual through every block's weights; measure how well the
## per-block held-out scores AGREE across blocks (mean pairwise RV). It never
## sees the truth, and it penalises over-regularisation BY CONSTRUCTION: when
## each block has fallen onto its own non-shared axis, the blocks disagree about
## where a held-out individual sits, and the score collapses.
##
## TWO DELIBERATE DEPARTURES from the calibration script, both documented in the
## .Rd and neither of them a change of criterion:
##   (1) the per-block weights are formed in the GRAM (dual) form,
##       B_j = X_j'(X_jX_j' + lambda_j I)^-1 Y, instead of the p x p primal form
##       of estimator.R:461. The push-through identity makes these the same
##       number; only the dual one survives p >> n, which is the regime mgcca is
##       for. Tested as an identity (test-select-lambda.R).
##   (2) columns with (near-)zero training-fold variance are dropped inside the
##       fold rather than producing NaN. The calibration's synthetic tables never
##       contained one; real tables do.
## =============================================================================


## ---- RV coefficient (rotation / sign / scale invariant) --------------------
## metrics.R:23-30. Two score configurations of the same individuals agree
## perfectly (1) when one is a rotation+rescaling of the other -- which is all a
## subspace is defined up to.
.mgcca_subspace_rv <- function(A, B) {
    A <- scale(as.matrix(A), center = TRUE, scale = FALSE)
    B <- scale(as.matrix(B), center = TRUE, scale = FALSE)
    SA <- tcrossprod(A); SB <- tcrossprod(B)
    den <- sqrt(sum(SA * SA) * sum(SB * SB))
    if (!is.finite(den) || den == 0) return(NA_real_)
    sum(SA * SB) / den
}

## ---- centre/scale one block, dropping degenerate columns -------------------
.mgcca_scale_block <- function(X, sd_tol = 1e-8) {
    X   <- as.matrix(X)
    ctr <- colMeans(X)
    sdv <- apply(X, 2L, stats::sd)
    keep <- is.finite(ctr) & is.finite(sdv) & sdv > sd_tol
    if (!any(keep))
        stop("a table has no column with non-zero variance in this fold")
    Xs <- scale(X[, keep, drop = FALSE],
                center = ctr[keep], scale = sdv[keep])
    list(X = Xs, center = ctr[keep], scale = sdv[keep], keep = keep,
         n_dropped = sum(!keep))
}

## ---- the in-memory penalized estimator (Gram form) -------------------------
## Port of estimator.R:33-77, which the thesis verified line by line against what
## mgcca computes for method = "penalized". Every object here is n x n or n x L:
## no p x p matrix is ever formed, so it runs unchanged at p >> n.
##
##   X_j  n x p_j, absent individuals padded with zero rows
##   G_j  = X_j X_j'                    R_j = (G_j + lambda_j I)^-1
##   M_j  = G_j R_j                     D   = diag(K_ii^-1/2), K_ii = #tables
##   S    = D (sum_j M_j) D             Y   = sqrt(J) D V_{1:L}
.mgcca_fit_penalized_dense <- function(tables, lambda, nfac, do_scale = TRUE,
                                       sd_tol = 1e-8) {
    J <- length(tables)
    if (length(lambda) == 1L) lambda <- rep(lambda, J)
    ids <- Reduce(union, lapply(tables, rownames))
    n   <- length(ids)
    nfac <- as.integer(nfac)
    if (nfac < 1L || nfac >= n)
        stop(sprintf("nfac (%d) must be in 1..n-1 (n = %d)", nfac, n))

    prep <- lapply(tables, function(m)
        if (isTRUE(do_scale)) .mgcca_scale_block(m, sd_tol)
        else list(X = as.matrix(m), center = NULL, scale = NULL,
                  keep = rep(TRUE, ncol(m)), n_dropped = 0L))

    Xlist <- present <- vector("list", J)
    for (j in seq_len(J)) {
        xi <- prep[[j]]$X
        Xj <- matrix(0, n, ncol(xi), dimnames = list(ids, colnames(xi)))
        Xj[rownames(xi), ] <- xi                      # absent rows stay zero
        Xlist[[j]]   <- Xj
        present[[j]] <- ids %in% rownames(xi)
    }
    Kcount <- Reduce(`+`, lapply(present, as.numeric))
    if (any(Kcount == 0))
        stop("some individual is present in no table")
    D <- 1 / sqrt(Kcount)

    Glist <- lapply(Xlist, tcrossprod)
    Rlist <- lapply(seq_len(J), function(j)
        solve(Glist[[j]] + lambda[j] * diag(n)))
    M <- Reduce(`+`, lapply(seq_len(J), function(j) Glist[[j]] %*% Rlist[[j]]))
    S <- outer(D, D) * M
    S <- (S + t(S)) / 2
    eig <- eigen(S, symmetric = TRUE)
    V   <- eig$vectors
    Y   <- sqrt(J) * (D * V[, seq_len(nfac), drop = FALSE])
    rownames(Y) <- ids

    list(ids = ids, n = n, J = J, nfac = nfac, lambda = lambda,
         Xlist = Xlist, present = present, D = D, Rlist = Rlist,
         V = V, mu = eig$values, Y = Y,
         center = lapply(prep, `[[`, "center"),
         scale  = lapply(prep, `[[`, "scale"),
         keep   = lapply(prep, `[[`, "keep"),
         n_dropped = vapply(prep, `[[`, integer(1), "n_dropped"))
}

## ---- per-block weights, in the dual form -----------------------------------
## B_j = (X_j'X_j + lambda_j I)^-1 X_j'Y  ==  X_j'(X_jX_j' + lambda_j I)^-1 Y.
## The right-hand side reuses the resolvent the fit already has, so a block with
## 400k variables costs one n x n solve that has already happened.
.mgcca_block_weights <- function(fit, j) {
    crossprod(fit$Xlist[[j]], fit$Rlist[[j]] %*% fit$Y)
}

## ---- the scale-matched block scale s_j -------------------------------------
## Port of estimator.R:96-124. s_j = tr(G_j)/rank(G_j) = the mean POSITIVE Gram
## eigenvalue, i.e. the scale of the eigenvalues M_j actually shrinks. A FIXED
## absolute lambda does not hold the amount of regularisation constant as p_j
## grows (tr(G_j) grows with p_j), so the dimensionless gamma_j = lambda_j / s_j
## is the quantity that is comparable across blocks and dimensions.
.mgcca_lambda_scale <- function(tables, do_scale = TRUE, rank_tol = 1e-8,
                                sd_tol = 1e-8) {
    vapply(seq_along(tables), function(j) {
        Xj <- if (isTRUE(do_scale)) .mgcca_scale_block(tables[[j]], sd_tol)$X
              else as.matrix(tables[[j]])
        ev <- eigen(tcrossprod(Xj), symmetric = TRUE, only.values = TRUE)$values
        ev <- ev[is.finite(ev)]
        ev[ev < 0] <- 0                               # G_j is PSD; clamp noise
        mx <- max(ev)
        if (mx <= 0)
            stop(sprintf("block %d has numerical rank 0; its scale is undefined", j))
        sum(ev) / sum(ev > rank_tol * mx)
    }, numeric(1))
}

## ---- one CV pass: mean held-out cross-block agreement ----------------------
## Port of cv_crossblock (estimator.R:445-476). `lambda_fold` is a list of length
## K holding the lambda vector to use in each fold (they differ when lambda is
## scale-matched, because s_j is recomputed from the training rows only).
.mgcca_cv_agreement <- function(tables, lambda_fold, nfac, fold, min_pair,
                                do_scale = TRUE, sd_tol = 1e-8) {
    ids <- names(fold)
    J   <- length(tables)
    K   <- length(lambda_fold)
    rv_by_fold <- rep(NA_real_, K)
    n_cells <- 0L
    for (k in seq_len(K)) {
        lam <- lambda_fold[[k]]
        if (is.null(lam)) next
        te <- ids[fold == k]
        trtab <- lapply(tables, function(m)
            m[!rownames(m) %in% te, , drop = FALSE])
        if (any(vapply(trtab, nrow, integer(1)) <= nfac + 1L)) next
        fit <- tryCatch(.mgcca_fit_penalized_dense(trtab, lam, nfac, do_scale,
                                                   sd_tol),
                        error = function(e) e)
        if (inherits(fit, "error")) next

        S <- vector("list", J)
        for (j in seq_len(J)) {
            raw <- tables[[j]][rownames(tables[[j]]) %in% te, , drop = FALSE]
            if (!nrow(raw)) next
            xs <- raw[, fit$keep[[j]], drop = FALSE]
            if (isTRUE(do_scale))
                xs <- scale(xs, center = fit$center[[j]], scale = fit$scale[[j]])
            sc <- xs %*% .mgcca_block_weights(fit, j)
            rownames(sc) <- rownames(raw)
            if (all(is.finite(sc))) S[[j]] <- sc
        }
        agrees <- numeric(0)
        for (a in seq_len(J - 1L)) for (b in (a + 1L):J) {
            if (is.null(S[[a]]) || is.null(S[[b]])) next
            com <- intersect(rownames(S[[a]]), rownames(S[[b]]))
            if (length(com) < min_pair) next
            v <- .mgcca_subspace_rv(S[[a]][com, , drop = FALSE],
                                    S[[b]][com, , drop = FALSE])
            if (is.finite(v)) { agrees <- c(agrees, v); n_cells <- n_cells + 1L }
        }
        if (length(agrees)) rv_by_fold[k] <- mean(agrees)
    }
    list(rv_by_fold = rv_by_fold, n_cells = n_cells)
}


#' Choose the ridge penalty by held-out cross-block agreement
#'
#' @description Selects the penalty \code{lambda} for
#'   \code{mgcca(method = "penalized")} by \eqn{K}-fold cross-validation over
#'   individuals, scoring each candidate by how well the per-block scores of
#'   \emph{held-out} individuals agree across blocks. It returns the chosen
#'   penalty, the agreement diagnostic and the whole candidate grid, so the
#'   choice lives in the user's record rather than inside a fit.
#'
#'   \code{mgcca()} does not call this function. \code{method = "penalized"}
#'   still requires an explicit \code{lambda}; nothing about the default
#'   behaviour of \code{\link{mgcca}} changes because this function exists.
#'
#' @section Why not choose lambda by stability:
#'   The intuitive rule -- keep the \code{lambda} whose fit moves least under
#'   resampling -- is not merely imperfect, it is actively misleading, and this
#'   function exists to replace it. Predicted instability is \emph{minimised} in
#'   the over-regularised regime, where each block has fallen onto its own
#'   high-variance but non-shared axis: the fit no longer moves because there is
#'   no longer anything shared in it that could move. In the specific-dominant
#'   calibration regime of the reference below, the stability-optimal
#'   \code{lambda} recovers the true shared subspace \strong{0.666 worse} than
#'   the best available \code{lambda}, while the criterion implemented here loses
#'   \strong{0.002}. Instability is a legitimate diagnostic \emph{conditional on}
#'   a chosen \eqn{(\lambda, L)}; it is not a selection rule.
#'
#'   Cross-block agreement rejects that regime by construction. When the blocks
#'   have latched onto their own private axes they disagree about where a
#'   held-out individual sits, so the score collapses exactly where the stability
#'   criterion is happiest.
#'
#' @section The criterion:
#'   For each fold: fit the penalized missing-row estimator on the training
#'   individuals; form each block's weights
#'   \eqn{B_j = X_j'(X_jX_j' + \lambda_j I)^{-1} Y} from the training fit;
#'   project the held-out individuals of block \eqn{j} through \eqn{B_j} using
#'   the \emph{training} centring and scaling; and average the RV coefficient
#'   between the held-out scores of every pair of blocks over the individuals
#'   both blocks hold. RV is invariant to rotation, sign and scale, which is all
#'   a subspace is defined up to. The score never sees any outcome or truth.
#'
#' @param x A named \code{list} of matrices (individuals \eqn{\times} variables)
#'   with \code{rownames} giving the individual IDs. Tables need not share
#'   individuals; they are aligned on the union, as in \code{\link{mgcca}}.
#' @param nfac Number of shared components the penalty is being chosen for.
#'   Default 2. The selected \code{lambda} is conditional on it.
#' @param gamma Candidate grid for the \emph{dimensionless} ridge
#'   (\code{scale_matched = TRUE}), where \eqn{\lambda_j = \gamma s_j} and
#'   \eqn{s_j = \mathrm{tr}(G_j)/\mathrm{rank}(G_j)} is the mean positive Gram
#'   eigenvalue of block \eqn{j}. Default \code{10^seq(-1, 1, by = 0.25)}.
#' @param lambda Candidate grid of absolute penalties
#'   (\code{scale_matched = FALSE}), applied identically to every block. Default
#'   \code{c(0.05, 0.2, 0.5, 1, 2, 5, 10, 30, 100, 300, 1000)} -- the grid of the
#'   reference calibration.
#' @param scale_matched If \code{TRUE} (default) the grid is over \code{gamma}
#'   and each block's penalty is put on the scale of that block's own Gram
#'   spectrum, recomputed inside every training fold. This is what keeps the
#'   amount of regularisation comparable across blocks of very different
#'   dimension; a fixed absolute \code{lambda} does not.
#' @param K Number of folds. Default 5.
#' @param seed Seed for the fold assignment, default 1. The whole procedure is a
#'   deterministic function of it. Pass \code{NULL} to leave the session's RNG
#'   stream untouched and let the folds fall where the current stream puts them.
#' @param rule \code{"argmax"} (default) takes the best mean CV agreement, as in
#'   the reference calibration. \code{"1se"} takes the most regularised candidate
#'   whose mean is within one standard error of the best.
#' @param scale If \code{TRUE} (default) each block is centred and scaled inside
#'   every fold using the training rows only.
#' @param min_pair Minimum number of held-out individuals a pair of blocks must
#'   share for that pair to be scored in a fold. Default \code{nfac + 3}.
#' @param rank_tol Relative tolerance for the numerical rank of \eqn{G_j} when
#'   \code{scale_matched = TRUE}. Default \code{1e-8}.
#' @param sd_tol Columns whose training-fold standard deviation is below this are
#'   dropped inside the fold. Default \code{1e-8}.
#'
#' @return An object of class \code{"mgcca_lambda"}: a list with
#'   \item{lambda}{the selected penalty, one per table, named.}
#'   \item{gamma}{the selected dimensionless ridge (\code{NA} when
#'     \code{scale_matched = FALSE}).}
#'   \item{s}{the per-block scales \eqn{s_j} on the full data (\code{NULL} when
#'     \code{scale_matched = FALSE}).}
#'   \item{agreement}{the mean held-out cross-block RV at the selection -- the
#'     diagnostic that justifies it.}
#'   \item{grid}{the candidate trace: the tuning value, \code{cv_mean},
#'     \code{cv_se}, \code{n_cells} (how many (fold, pair) cells were scorable),
#'     \code{in_1se} and \code{selected}.}
#'   \item{cv_by_fold}{the raw fold-by-candidate agreement matrix.}
#'   \item{lambda_by_candidate}{the penalty each candidate corresponds to.}
#'   \item{status}{\code{"ok"}, or \code{"boundary_optimum"} when the best
#'     candidate is at an edge of the grid, which is also warned about. When no
#'     candidate can be scored at all the function stops rather than returning a
#'     choice nothing supports.}
#'   \item{criterion, rule, nfac, K, seed, fold, call}{the record of how the
#'     choice was made.}
#'
#' @section Reporting it:
#'   Report the selected \code{lambda} together with \code{agreement} and the
#'   grid. A wide near-optimal region (many candidates inside one standard error)
#'   means the penalty is not sharply identified and the region, not the single
#'   number, is the honest summary; the function says so.
#'
#' @references
#'   The criterion, the guardrail and the 0.002 / 0.666 comparison come from the
#'   thesis reliability calibration
#'   (\code{reliability/05_lambda_selection.md}; estimator in
#'   \code{reliability/calibration/estimator.R}). van de Velden, M. and Bijmolt,
#'   T. H. A. (2006) Generalized canonical correlation analysis of matrices with
#'   missing rows. \emph{Psychometrika} 71, 323-331.
#' @seealso \code{\link{mgcca}}
#' @examples
#' ## A small three-block slice of the shipped cohort. The tables overlap only
#' ## partly, which is the situation the criterion is written for.
#' data(cardiovascular, package = "mgcca")
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:150]
#' mk  <- function(d, j = NULL) {
#'     m <- as.matrix(d[rownames(d) %in% ids, , drop = FALSE])
#'     if (!is.null(j)) m <- m[, j, drop = FALSE]
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = mk(X1, 1:20), clinical = mk(X2), cells = mk(X3))
#'
#' ## Three candidates for the dimensionless ridge, three folds. A real
#' ## analysis uses the wider default grid; this one is sized for the example.
#' sel <- mgcca_select_lambda(X, nfac = 2, gamma = c(0.01, 0.1, 1), K = 3)
#' sel
#'
#' ## The whole candidate trace, not only the winner: report the near-optimal
#' ## region whenever several candidates sit within one standard error.
#' sel$grid
#'
#' ## The selected penalty is what mgcca(method = "penalized") then takes.
#' fit <- mgcca(X, nfac = 2, method = "penalized", lambda = sel$lambda)
#' @export
mgcca_select_lambda <- function(x, nfac = 2, gamma = NULL, lambda = NULL,
                                scale_matched = TRUE, K = 5, seed = 1,
                                rule = c("argmax", "1se"), scale = TRUE,
                                min_pair = NULL, rank_tol = 1e-8,
                                sd_tol = 1e-8) {

    cl   <- match.call()
    rule <- match.arg(rule)

    ## ---- what we were given ------------------------------------------------
    if (!is.list(x) || length(x) < 2L)
        stop("'x' must be a list of at least two tables: cross-block agreement ",
             "needs a pair of blocks to compare")
    tables <- lapply(x, as.matrix)
    if (is.null(names(tables)) || any(!nzchar(names(tables))))
        names(tables) <- paste0("table", seq_along(tables))
    if (any(vapply(tables, function(m) is.null(rownames(m)), logical(1))))
        stop("every table needs rownames: they are the individual IDs the ",
             "tables are aligned on")
    nfac <- suppressWarnings(as.integer(nfac))
    ids  <- Reduce(union, lapply(tables, rownames))
    if (is.na(nfac) || length(nfac) != 1L || nfac < 1L || nfac >= length(ids))
        stop("'nfac' must be a single integer in 1..n-1 (n = ", length(ids), ")")
    K <- as.integer(K)
    if (is.na(K) || K < 2L) stop("'K' must be at least 2")
    if (is.null(min_pair)) min_pair <- nfac + 3L       # estimator.R:470
    J <- length(tables)

    if (!is.null(gamma) && !is.null(lambda))
        stop("give a grid over 'gamma' or over 'lambda', not both")
    if (isTRUE(scale_matched) && !is.null(lambda))
        stop("'lambda' is the grid for scale_matched = FALSE; use 'gamma' ",
             "for the scale-matched parameterisation")
    grid_vals <- if (isTRUE(scale_matched)) {
        if (is.null(gamma)) 10^seq(-1, 1, by = 0.25) else gamma
    } else {
        if (is.null(lambda)) c(0.05, 0.2, 0.5, 1, 2, 5, 10, 30, 100, 300, 1000)
        else lambda
    }
    grid_vals <- sort(unique(as.numeric(grid_vals)))
    if (!length(grid_vals) || any(!is.finite(grid_vals)) || any(grid_vals <= 0))
        stop("the candidate grid must be finite and strictly positive")
    G <- length(grid_vals)

    ## ---- folds (estimator.R:449-450) ---------------------------------------
    ## The seed is the caller's: it is what makes the fold assignment, and so the
    ## whole selection, reproducible. `seed = NULL` opts out and leaves the
    ## session's RNG stream alone; the default `seed = 1` is unchanged.
    if (!is.null(seed)) set.seed(seed)
    fold <- sample(rep(seq_len(K), length.out = length(ids)))
    names(fold) <- ids

    ## ---- per-fold block scales, from the TRAINING rows only ----------------
    s_fold <- if (isTRUE(scale_matched))
        lapply(seq_len(K), function(k) {
            trtab <- lapply(tables, function(m)
                m[!rownames(m) %in% ids[fold == k], , drop = FALSE])
            tryCatch(.mgcca_lambda_scale(trtab, scale, rank_tol, sd_tol),
                     error = function(e) NULL)
        }) else vector("list", K)

    ## ---- the CV surface ----------------------------------------------------
    cv_by_fold <- matrix(NA_real_, G, K,
                         dimnames = list(NULL, paste0("fold", seq_len(K))))
    n_cells <- integer(G)
    for (gi in seq_len(G)) {
        lam_fold <- lapply(seq_len(K), function(k)
            if (isTRUE(scale_matched)) {
                if (is.null(s_fold[[k]])) NULL else grid_vals[gi] * s_fold[[k]]
            } else rep(grid_vals[gi], J))
        cv <- .mgcca_cv_agreement(tables, lam_fold, nfac, fold, min_pair,
                                  scale, sd_tol)
        cv_by_fold[gi, ] <- cv$rv_by_fold
        n_cells[gi]      <- cv$n_cells
    }
    fmean <- function(r) { r <- r[is.finite(r)]
                           if (length(r)) mean(r) else NA_real_ }
    fse   <- function(r) { r <- r[is.finite(r)]
                           if (length(r) > 1L) stats::sd(r) / sqrt(length(r))
                           else NA_real_ }
    cv_mean <- apply(cv_by_fold, 1L, fmean)
    cv_se   <- apply(cv_by_fold, 1L, fse)

    s_full <- if (isTRUE(scale_matched))
        .mgcca_lambda_scale(tables, scale, rank_tol, sd_tol) else NULL
    lam_by_cand <- if (isTRUE(scale_matched))
        outer(grid_vals, s_full) else matrix(grid_vals, G, J)
    dimnames(lam_by_cand) <- list(NULL, names(tables))

    ## ---- the pick ----------------------------------------------------------
    if (all(!is.finite(cv_mean)))
        stop("no candidate could be scored: no (fold, pair) had at least ",
             min_pair, " held-out individuals in common. Lower 'K' or ",
             "'min_pair', or accept that the blocks share too few individuals.")
    best <- which.max(cv_mean)
    thr  <- cv_mean[best] - (if (is.finite(cv_se[best])) cv_se[best] else 0)
    in_1se <- is.finite(cv_mean) & cv_mean >= thr
    sel <- if (identical(rule, "1se")) max(which(in_1se)) else best

    selected <- logical(G); selected[sel] <- TRUE
    grid <- data.frame(cv_mean = cv_mean, cv_se = cv_se, n_cells = n_cells,
                       in_1se = in_1se, selected = selected)
    grid <- cbind(stats::setNames(data.frame(grid_vals),
                                  if (isTRUE(scale_matched)) "gamma" else "lambda"),
                  grid)
    row.names(grid) <- NULL

    ## ---- honesty about the shape of the surface ----------------------------
    status <- if (best %in% c(1L, G)) "boundary_optimum" else "ok"
    if (identical(status, "boundary_optimum")) {
        ## The message is assembled first because warningCondition() takes ONE
        ## `message`: anything further in `...` becomes a field of the condition
        ## object, not more text. Same string as before, built outside the call.
        boundary_msg <- paste0(
            "the best candidate is at the ", if (best == 1L) "lower" else "upper",
            " edge of the grid (", format(grid_vals[best], digits = 4),
            "): the criterion has not turned over, so the optimum may lie ",
            "outside it.\n  Widen the grid in that direction, or report the ",
            "near-optimal region rather than a single value.")
        warning(warningCondition(boundary_msg,
            class = "mgcca_lambda_boundary", call = NULL))
    }
    flat <- sum(in_1se, na.rm = TRUE) >= max(3L, floor(G / 2))
    if (flat && identical(status, "ok"))
        message("mgcca_select_lambda: ", sum(in_1se, na.rm = TRUE),
                " of ", G, " candidates are within one standard error of the ",
                "best; the penalty is not sharply identified. Report the region.")

    out <- list(
        lambda    = stats::setNames(as.numeric(lam_by_cand[sel, ]), names(tables)),
        gamma     = if (isTRUE(scale_matched)) grid_vals[sel] else NA_real_,
        s         = if (isTRUE(scale_matched))
                        stats::setNames(s_full, names(tables)) else NULL,
        agreement = unname(cv_mean[sel]),
        grid      = grid,
        cv_by_fold = cv_by_fold,
        lambda_by_candidate = lam_by_cand,
        threshold = unname(thr),
        status    = status,
        flat_1se  = flat,
        criterion = "cv_crossblock_agreement",
        rule      = rule,
        nfac      = nfac,
        K         = K,
        seed      = seed,
        min_pair  = as.integer(min_pair),
        scale_matched = isTRUE(scale_matched),
        fold      = fold,
        call      = cl)
    class(out) <- "mgcca_lambda"
    out
}

#' Print a lambda selection
#'
#' @param x an \code{mgcca_lambda} object, from
#'   \code{\link{mgcca_select_lambda}}.
#' @param ... ignored.
#' @return \code{x}, invisibly. Called for the one-screen summary it writes to
#'   the console: the criterion, the selected penalty, the held-out cross-block
#'   agreement that justifies it, and how many candidates lie within one
#'   standard error of the best.
#' @seealso \code{\link{mgcca_select_lambda}}
#' @examples
#' data(cardiovascular, package = "mgcca")
#' ids <- Reduce(union, list(rownames(X1), rownames(X2), rownames(X3)))[1:150]
#' mk  <- function(d, j = NULL) {
#'     m <- as.matrix(d[rownames(d) %in% ids, , drop = FALSE])
#'     if (!is.null(j)) m <- m[, j, drop = FALSE]
#'     storage.mode(m) <- "double"
#'     m
#' }
#' X <- list(methylation = mk(X1, 1:20), clinical = mk(X2), cells = mk(X3))
#'
#' sel <- mgcca_select_lambda(X, nfac = 2, gamma = c(0.01, 0.1, 1), K = 3)
#' print(sel)
#' @export
#' @method print mgcca_lambda
print.mgcca_lambda <- function(x, ...) {
    cat("mgcca ridge selection by held-out cross-block agreement\n")
    cat(sprintf("  criterion   : %s (%d-fold over individuals, rule = %s)\n",
                x$criterion, x$K, x$rule))
    cat(sprintf("  components  : %d\n", x$nfac))
    if (x$scale_matched)
        cat(sprintf("  gamma       : %.4g   (lambda_j = gamma * s_j)\n", x$gamma))
    cat(sprintf("  lambda      : %s\n",
                paste(sprintf("%s = %.4g", names(x$lambda), x$lambda),
                      collapse = ", ")))
    cat(sprintf("  agreement   : %.4f   (mean held-out cross-block RV)\n",
                x$agreement))
    cat(sprintf("  near-optimal: %d of %d candidates within 1 SE\n",
                sum(x$grid$in_1se, na.rm = TRUE), nrow(x$grid)))
    if (!identical(x$status, "ok"))
        cat(sprintf("  status      : %s\n", x$status))
    cat("  Instability is a diagnostic at this choice, never the way to make it.\n")
    invisible(x)
}
