## In-memory estimator: the mgcca() path taken when filename is NULL.
##
## This is a pure-R fit for tables that fit in RAM. It is the union of two pieces
## that were already validated against the port oracle: the covariance pipeline of
## tests/harness/01_oracle.R (route "cov", p < n) and the Gram/dual form of
## reliability/calibration/estimator.R::fit_pgcca extended to the four methods and
## to scores/corsY/AVE (route "dual", p >= n). Every quantity is produced with the
## SAME dispatch and the SAME formulas the C++/HDF5 core uses, so the returned
## object is numerically equivalent (to the declared tolerances) and structurally
## identical to what mgcca(collect = TRUE) returns from a file. src/ is never
## touched: the file-backed path is left bit-for-bit as it was.
##
## Parity-deciding details reproduced from the C++ engine (src/mgcca_phases.hpp):
##   * union order: first appearance when all tables share the row count,
##     alphabetical otherwise (run_getK).
##   * scaling: base-R scale() (sd with n-1) over the PRESENT rows, before padding.
##   * eigen: symmetric direct solver, descending order, keig = nfac+1 so the
##     external gap mu_L - mu_{L+1} is defined.
##   * corsY over present individuals; pval = 2 * pt(|t|, n_present - 2).
##   * dual weights with tol = lambda_max * 1e-12 (run_svd / run_scores_svd).

## Two-tailed correlation p-value, exactly as the oracle / the C++ core: the df is
## the number of PRESENT individuals minus two (not the union size).
.mgcca_cor_pval_mem <- function(r, n_present) {
    df <- n_present - 2
    t  <- r * sqrt(df) / sqrt(1 - r * r)
    2 * stats::pt(abs(t), df, lower.tail = FALSE)
}

## Per-table Mj (projector), plus whatever the scores stage needs to reuse. Which
## route is taken is decided exactly as the C++ orchestrator decides it.
.mgcca_mem_block <- function(Xj, inv, lam_j, dual) {
    if (dual) {
        # Gram form: G = X X' (m x m), eigen of G gives the left singular vectors U
        # (lambda_i = d_i^2). M_j = U diag(w) U'. Never forms p x p.
        G  <- base::tcrossprod(Xj)
        eg <- eigen(G, symmetric = TRUE)
        li <- eg$values; U <- eg$vectors
        lmax <- max(li); tol <- (if (lmax > 0) lmax else 1) * 1e-12
        w  <- if (inv == 2L) li / (li + lam_j) else as.numeric(li > tol)
        list(Mi = U %*% (w * t(U)),
             cache = list(dual = TRUE, U = U, li = li, tol = tol))
    } else {
        # Covariance form: Mgram = X'X (= X'KX after padding), inverse per method,
        # M_j = X xkx X'. The padded (absent) rows of X are zero, so K drops out.
        xk  <- t(Xj)                                    # p x m
        Mg  <- xk %*% Xj                                # X'X (p x p)
        xkx <- if (inv == 1L) chol2inv(chol(Mg))
               else if (inv == 2L) chol2inv(chol(Mg + base::diag(nrow(Mg)) * lam_j))
               else MASS::ginv(Mg)                      # geninv / ginv
        list(Mi = Xj %*% xkx %*% xk,
             cache = list(dual = FALSE, xkx = xkx, xk = xk))
    }
}

## The in-memory driver. Arguments mirror mgcca(); `inv` is the resolved integer
## method code (1 solve, 2 penalized, 3 geninv/ginv), as in the C++ core.
.mgcca_memory <- function(x, group, datasets, nfac, scale, method, lambda,
                          inv, scores, route, collect, outputs) {

    if (is.character(x) && length(x) == 1L)
        stop("in-memory mode cannot read from an HDF5 file path; supply ",
             "'filename = ' to run the file-backed pipeline on '", x, "'.",
             call. = FALSE)
    if (!isTRUE(collect))
        stop("collect = FALSE has no meaning for an in-memory fit: there is no ",
             "file to leave the results in. Use collect = TRUE, or supply a ",
             "'filename' to run the HDF5-backed pipeline instead.", call. = FALSE)

    tables <- .mgcca_as_tables(x)                       # named list, ind x vars
    nms <- names(tables)
    if (!is.null(datasets)) {
        if (!all(datasets %in% nms))
            stop("'datasets' must be a subset of the table names in 'x'", call. = FALSE)
        tables <- tables[datasets]; nms <- datasets
    }
    if (length(tables) < 2L) stop("mgcca needs at least two tables", call. = FALSE)
    for (nm in nms)
        if (is.null(rownames(tables[[nm]])))
            stop("table '", nm, "' has no rownames (individual IDs required)",
                 call. = FALSE)

    nfac <- as.integer(nfac)
    J    <- length(tables)
    lam  <- if (inv == 2L) { if (length(lambda) == 1L) rep(lambda, J) else lambda }
            else rep(0, J)
    if (inv == 2L && length(lam) != J)
        stop("method 'penalized' requires 'lambda' of length ", J, call. = FALSE)

    ## ---- scaling over present rows, before padding (matches run_normalize) ---
    scaling <- NULL
    Xs <- tables
    if (isTRUE(scale)) {
        scaling <- vector("list", J); names(scaling) <- nms
        for (j in seq_len(J)) {
            sc  <- scale(tables[[j]])                   # base R: sd with (n-1)
            cen <- attr(sc, "scaled:center"); sca <- attr(sc, "scaled:scale")
            attr(sc, "scaled:center") <- NULL; attr(sc, "scaled:scale") <- NULL
            Xs[[j]] <- sc
            csm <- cbind(center = cen, scale = sca)
            rownames(csm) <- colnames(tables[[j]])
            scaling[[j]] <- csm
        }
    }

    ## ---- getK: union of IDs, padding, presence ------------------------------
    ns <- vapply(Xs, nrow, integer(1))
    ids_list <- lapply(Xs, rownames)
    rn <- if (max(ns) == min(ns)) Reduce(union, ids_list)
          else sort(Reduce(union, ids_list))
    m  <- length(rn)
    pad <- function(xi) {
        na.ids <- rn[!rn %in% rownames(xi)]
        xna <- matrix(0, length(na.ids), ncol(xi),
                      dimnames = list(na.ids, colnames(xi)))
        as.matrix(rbind(xi, xna)[rn, , drop = FALSE])
    }
    Xp <- lapply(Xs, pad)                               # padded, m x p_j
    Kd <- lapply(Xs, function(xi) as.integer(rn %in% rownames(xi)))
    p  <- vapply(Xp, ncol, integer(1))

    ## ---- per-table route + Mj -----------------------------------------------
    route_dual <- vapply(seq_len(J), function(j)
        switch(route, cov = FALSE, dual = TRUE, auto = (p[j] >= m)), logical(1))
    names(route_dual) <- nms
    blocks <- lapply(seq_len(J), function(j)
        .mgcca_mem_block(Xp[[j]], inv, lam[j], route_dual[j]))

    ## ---- operator, eigen, Y -------------------------------------------------
    M       <- Reduce(`+`, lapply(blocks, `[[`, "Mi"))
    Ksum    <- Reduce(`+`, Kd)
    Ksum05  <- Ksum^(-0.5)
    MKsum05 <- (Ksum05 * M) * rep(Ksum05, each = m)     # D M D (whitened operator)

    keig <- min(max(nfac + 1L, 2L), m)
    e    <- eigen(MKsum05, symmetric = TRUE)            # descending
    valk <- e$values[seq_len(keig)]
    Vk   <- e$vectors[, seq_len(keig), drop = FALSE]
    Yast <- Vk[, seq_len(nfac), drop = FALSE]
    comp <- paste0("comp", seq_len(nfac))
    Y    <- sqrt(J) * (Ksum05 * Yast)                   # sqrt(J) K^{-1/2} Y*
    dimnames(Y) <- list(rn, comp)

    ## eigen diagnostics: external gap + relative eigen-residual (rel_resid guard)
    SV   <- MKsum05 %*% Yast
    rnum <- sqrt(sum((SV - Yast * rep(valk[seq_len(nfac)], each = m))^2))
    rden <- sqrt(sum(SV^2))
    reig <- if (rden <= 1e-30) (if (rnum <= 1e-30) 0 else NA_real_) else rnum / rden
    have <- length(valk) > nfac
    gabs <- if (have) valk[nfac] - valk[nfac + 1L] else NA_real_
    grel <- if (have) gabs / max(abs(valk[nfac]), 1e-300) else NA_real_
    eig_report <- list(gap = gabs, gap_rel = grel, residual = reig)

    ## ---- corsY / pval / AVE (present individuals only) -----------------------
    corsY <- pvalc <- vector("list", J); names(corsY) <- names(pvalc) <- nms
    AVE_X <- matrix(0, nfac, J)
    for (j in seq_len(J)) {
        pres <- Kd[[j]] == 1L; npr <- sum(pres)
        if (npr < 3L)
            stop("block '", nms[j], "' has fewer than three present individuals; ",
                 "no p-value is defined", call. = FALSE)
        C <- stats::cor(Xp[[j]][pres, , drop = FALSE], Y[pres, , drop = FALSE])
        dimnames(C) <- list(colnames(Xp[[j]]), comp)
        corsY[[j]] <- C
        Pv <- .mgcca_cor_pval_mem(C, npr)
        dimnames(Pv) <- list(colnames(Xp[[j]]), comp)
        pvalc[[j]] <- Pv
        AVE_X[, j] <- colMeans(C^2)
    }
    AVE <- list(AVE_X           = AVE_X,
                AVE_outer_model = as.numeric((AVE_X %*% p) / sum(p)),
                AVE_inner_model = valk[seq_len(nfac)])

    ## ---- scores (optional) --------------------------------------------------
    scores_l <- weights_l <- NULL
    if (isTRUE(scores)) {
        scores_l <- weights_l <- vector("list", J)
        names(scores_l) <- names(weights_l) <- nms
        for (j in seq_len(J)) {
            Xj <- Xp[[j]]; ch <- blocks[[j]]$cache; pres <- Kd[[j]] == 1L
            if (ch$dual) {
                wA <- if (inv == 2L) 1 / (ch$li + lam[j])
                      else ifelse(ch$li > ch$tol, 1 / ch$li, 0)
                inner <- ch$U %*% (wA * (t(ch$U) %*% Y))   # m x nfac
                A <- t(Xj) %*% inner                        # p x nfac
            } else {
                A <- ch$xkx %*% (ch$xk %*% Y)               # p x nfac
            }
            KXA <- Xj %*% A                                 # m x nfac
            vv  <- apply(KXA[pres, , drop = FALSE], 2, stats::sd)  # sd over present
            As  <- A / rep(vv, each = nrow(A))
            sc  <- Xj %*% As
            sc[!pres, ] <- NA_real_
            dimnames(sc) <- list(rn, comp)
            dimnames(As) <- list(colnames(Xj), comp)
            scores_l[[j]]  <- sc
            weights_l[[j]] <- As
        }
    }

    ## ---- overlap report (from the presence masks held in RAM) ---------------
    P   <- matrix(unlist(Kd) != 0, nrow = m, ncol = J, dimnames = list(NULL, nms))
    ovl <- .mgcca_overlap_core(P, nms)

    ## ---- assemble the object, honouring outputs = (subtractive) -------------
    want <- .mgcca_output_plan(outputs, scores = scores, pval = TRUE)
    ans <- list(
        Y        = if (want[["Y"]])       Y          else NULL,
        corsY    = if (want[["corsY"]])   corsY      else NULL,
        scores   = if (want[["scores"]])  scores_l   else NULL,
        pval.cor = if (want[["pval"]])    pvalc      else NULL,
        weights  = if (want[["weights"]]) weights_l  else NULL,
        scaling  = if (want[["scaling"]]) scaling    else NULL,
        AVE      = if (want[["AVE"]])     AVE        else NULL,
        eigen    = if (want[["eigen"]])   eig_report else NULL,
        overlap  = if (want[["overlap"]]) ovl        else NULL)
    class(ans) <- "mgcca"
    if (!is.null(outputs)) attr(ans, "outputs") <- .mgcca_output_canonical(outputs)

    ## ---- descriptor: same shape as the file-backed one, memory backend ------
    desc <- list(
        filename      = NA_character_,
        group         = NA_character_,
        datasets      = nms,
        final_group   = "FINAL_RESULTS",
        nfac          = nfac,
        m             = m,
        eig_values    = valk[seq_len(nfac)],
        route         = route,
        route_dual    = route_dual,
        scores        = isTRUE(scores),
        method        = method,
        lambda        = if (inv == 2L) lam else NULL,
        input_group   = NA_character_,
        scale         = isTRUE(scale),
        backend       = "memory",
        eigen         = eig_report,
        mgcca_version = as.character(utils::packageVersion("mgcca")))
    attr(ans, "desc") <- desc
    ans
}

## Early, explicit refusal for the reliability layer on an in-memory fit. It is
## called at the very top of mgcca_sensitivity() / mgcca_stability() so the frozen
## logic of those functions is left untouched; a file-backed fit passes straight
## through (backend is NULL for those). The reliability kernels read the SOURCE
## blocks from HDF5, which an in-memory fit does not have.
.mgcca_refuse_memory <- function(x) {
    d <- attr(x, "desc")
    if (inherits(x, "mgcca") && !is.null(d) && identical(d$backend, "memory"))
        stop("reliability analyses need the HDF5-backed source blocks; ",
             "re-run mgcca() with a filename", call. = FALSE)
    invisible(TRUE)
}
