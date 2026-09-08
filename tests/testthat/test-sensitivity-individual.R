# --- Bioconductor build time -------------------------------------------------
# This file is one of the four heaviest in the suite. Bioconductor's builders
# cap `R CMD check` at 10 minutes and this suite is the long pole, so the file
# is skipped THERE only (`IS_BIOC_BUILD_MACHINE`); it runs in full everywhere
# else, including on CRAN-style checks and in development.
testthat::skip_on_bioc()

# Per-individual sensitivity decomposition: the `by_individual` element.
#
# The fixture is the one the empirical probe used, reproduced line for line
# (seed, sizes, missingness patterns, ridge, query and grouping). Reproducing it
# is not decoration: the reference magnitudes asserted below were measured on
# THIS fixture, so a different one would leave the tolerances unanchored.
#
# The oracle is plain R over the package's own internal algebra. It shares the
# reference fit with the public path -- comparing against a differently
# conditioned basis would measure the basis, not the decomposition.
#
# ⚠️ TOLERANCE DOCTRINE, frozen. Parity and the decomposition identity are
# checked as  max|delta| <= 1e-12 * max(1, scale)  where `scale` is a BLOCK-level
# magnitude: max(abs(oracle accX_j)) for parity, |S_j| for the identity. A
# per-entry relative criterion is never used, because a participant absent from a
# block contributes ~1e-34 there and a relative error against that measures
# rounding against nothing.

perind_fit <- function(tabs, L = 2L) {
    h5 <- tempfile(fileext = ".h5")
    fit <- suppressMessages(mgcca(tabs, filename = h5, nfac = L, scale = TRUE,
                                  method = "solve", scores = TRUE))
    BigDataStatMeth::hdf5_close_all()
    fit
}

# Built once: the fit is small but every test needs the same one, and refitting
# would only add wall time.
.perind_cache <- new.env(parent = emptyenv())

perind_fixture <- function() {
    if (!is.null(.perind_cache$fx)) return(.perind_cache$fx)
    set.seed(20260903)
    n <- 70L; J <- 3L; L <- 2L
    ids <- sprintf("id%03d", seq_len(n))
    p <- c(blkA = 12L, blkB = 9L, blkC = 15L)
    lambda <- c(0.8, 1.2, 1.0)
    Z <- matrix(stats::rnorm(n * L), n, L)
    Xfull <- lapply(seq_len(J), function(j) {
        X <- Z %*% matrix(stats::rnorm(L * p[j]), L, p[j]) +
            matrix(stats::rnorm(n * p[j]), n, p[j])
        dimnames(X) <- list(ids, sprintf("%s_v%02d", names(p)[j], seq_len(p[j])))
        X
    })
    names(Xfull) <- names(p)
    av <- list(blkA = rep(TRUE, n),
               blkB = !(ids %in% sprintf("id%03d", 61:68)),
               blkC = !(ids %in% sprintf("id%03d", c(5, 17, 23, 42, 55))))
    av$blkC[ids == "id070"] <- FALSE
    tabs <- Map(function(X, m) X[m, , drop = FALSE], Xfull, av)
    set.seed(20260905)
    fx <- list(n = n, L = L, ids = ids, Xfull = Xfull, av = av, tabs = tabs,
               lambda = lambda,
               query = stats::setNames(as.numeric(Z[, 1] + 0.5 * stats::rnorm(n)), ids),
               group = stats::setNames(rep(c("g1", "g2"), each = n / 2), ids))
    fx$fit <- perind_fit(tabs, L)
    .perind_cache$fx <- fx
    fx
}

# The probe's `oraculo_perind`, in plain R over the package's internals. Returns
# the per-participant table and, per block, the scale the parity tolerance uses.
perind_oracle <- function(fitobj, lambda, query, group) {
    ctx <- mgcca:::.mgcca_rel_context(fitobj)
    ids <- ctx$ids
    blocks <- suppressMessages(mgcca:::.mgcca_rel_block_grams(ctx, "memory", 0L, NULL))
    fr <- mgcca:::reliability_reference_fit(blocks$Glist, blocks$present, lambda, ctx$L)
    BigDataStatMeth::hdf5_close_all()
    Q <- mgcca:::.mgcca_rel_query_matrix(query, ids)
    zs <- Q[, 1L] / fr$D
    G <- mgcca:::.mgcca_rel_group(group, ids)
    wrap <- list(V = fr$V, mu = fr$mu, D = fr$D, Rlist = fr$Rlist,
                 L = ctx$L, n = length(ids), lambda = lambda)
    A <- tcrossprod(zs) / sum(zs^2)
    C <- mgcca:::.mgcca_rel_C(wrap, A)
    prm0 <- G$index; grp0 <- G$code
    out <- lapply(seq_along(ctx$datasets), function(j) {
        bi <- mgcca:::.mgcca_rel_block_inputs(wrap, j)
        A0 <- bi$Rr %*% C %*% t(bi$Sm)
        Ho <- bi$alpha * (A0 + t(A0))
        Tm <- blocks$Glist[[j]][prm0 + 1L, , drop = FALSE] %*% Ho
        Tm[, !blocks$present[[j]]] <- 0
        cm <- colMeans(Tm)
        accX <- sweep(Tm, 2, cm)
        gs <- rowsum(Tm, grp0) / as.vector(table(grp0))
        accB <- gs[grp0 + 1L, , drop = FALSE] -
            matrix(cm, nrow(Tm), ncol(Tm), byrow = TRUE)
        list(tab = data.frame(id = ids[prm0 + 1L], block = ctx$datasets[j],
                              S_total = rowSums(accX^2),
                              S_between = rowSums(accB^2),
                              stringsAsFactors = FALSE),
             scale_total = max(abs(accX)), scale_between = max(abs(accB)))
    })
    d <- do.call(rbind, lapply(out, `[[`, "tab")); rownames(d) <- NULL
    list(tab = d,
         scale = data.frame(block = ctx$datasets,
                            scale_total = vapply(out, `[[`, numeric(1), "scale_total"),
                            scale_between = vapply(out, `[[`, numeric(1), "scale_between"),
                            stringsAsFactors = FALSE))
}

perind_key <- function(d) paste(d$query, d$id, d$block, sep = "|")
perind_bkey <- function(d) paste(d$id, d$block, sep = "|")

# Per-participant total, summed over blocks, in a stable id order.
perind_total <- function(s) {
    v <- rowsum(s$by_individual$S_total, s$by_individual$id)
    stats::setNames(v[, 1L], rownames(v))
}

perind_sens <- function(...) suppressMessages(mgcca_sensitivity(...))

test_that("the per-individual table matches a plain-R oracle at block scale", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    bi <- s$by_individual

    o <- perind_oracle(f$fit, f$lambda, f$query, f$group)
    i <- match(perind_bkey(bi), perind_bkey(o$tab))
    expect_false(anyNA(i))
    j <- match(bi$block, o$scale$block)
    tol_t <- 1e-12 * pmax(1, o$scale$scale_total[j])
    tol_b <- 1e-12 * pmax(1, o$scale$scale_between[j])
    d_tot <- abs(bi$S_total - o$tab$S_total[i])
    d_btw <- abs(bi$S_between - o$tab$S_between[i])
    expect_true(all(d_tot <= tol_t))
    expect_true(all(d_btw <= tol_b))
    cat(sprintf("\n[perind] oracle parity: max|dS_total| = %.3e, max|dS_between| = %.3e (tol >= %.1e)\n",
                max(d_tot), max(d_btw), min(tol_t)))

    # Finite and non-negative: these are squared row norms.
    expect_true(all(is.finite(bi$S_total)) && all(bi$S_total >= 0))
    expect_true(all(is.finite(bi$S_between)) && all(bi$S_between >= 0))

    # Absent-from-block entries are a NUMERICAL zero, not an exact one, and are
    # deliberately NOT thresholded away.
    absent <- bi$S_total[bi$block == "blkC" & bi$id == "id070"]
    expect_lt(absent, 1e-20)
})

test_that("A3: rows carry the right participant, and the contribution follows the ID", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    bi <- s$by_individual

    # Exact ID vector and row count: the public rows ARE the ordered participant
    # set of the grouped calculation. Certified here, not promised in the docs.
    ctx <- mgcca:::.mgcca_rel_context(f$fit)
    G <- mgcca:::.mgcca_rel_group(f$group, ctx$ids)
    expect_identical(nrow(bi), 70L * 3L)
    for (b in ctx$datasets)
        expect_identical(bi$id[bi$block == b], ctx$ids[G$index + 1L])

    # (query, id, block) is a key.
    expect_identical(anyDuplicated(perind_key(bi)), 0L)

    # ⚠️ THE GATE THIS RELEASE EXISTS FOR. A row permutation of the kernel's
    # accumulator leaves every block TRACE unchanged, so 1.2.0 could not have
    # noticed one. Attaching those rows to participant IDs makes the same
    # permutation produce perfect block totals with wrong participant labels.
    # So: permute the stored rows of every table and of the query and grouping,
    # refit, and require the contribution to follow the PARTICIPANT.
    set.seed(20260906)
    tabsP <- lapply(f$tabs, function(X) X[sample.int(nrow(X)), , drop = FALSE])
    qP <- f$query[sample.int(f$n)]
    gP <- f$group[sample.int(f$n)]
    fitP <- perind_fit(tabsP, f$L)
    sP <- perind_sens(fitP, query = qP, group = gP, lambda = f$lambda,
                      backend = "memory")
    biP <- sP$by_individual

    # The permutation must actually have reached the kernel, or the test is
    # vacuous. It does: the blocks are STORED in the permuted order, and that is
    # the order the kernel indexes and maps through `mids_to_fit`.
    expect_false(identical(lapply(tabsP, rownames), lapply(f$tabs, rownames)))
    ctxP <- mgcca:::.mgcca_rel_context(fitP)
    bP <- suppressMessages(mgcca:::.mgcca_rel_block_grams(ctxP, "memory", 0L, NULL))
    b0 <- suppressMessages(mgcca:::.mgcca_rel_block_grams(ctx, "memory", 0L, NULL))
    BigDataStatMeth::hdf5_close_all()
    expect_false(identical(bP$block_ids, b0$block_ids))
    # The EXPOSED row order, by contrast, is canonical -- which is the point:
    # consumers key by (query, id, block), and the contribution follows the
    # participant wherever that participant's data happened to be stored.
    expect_setequal(perind_bkey(biP), perind_bkey(bi))
    i <- match(perind_bkey(biP), perind_bkey(bi))
    d <- max(abs(biP$S_total - bi$S_total[i]))
    expect_lt(d, 1e-12)
    cat(sprintf("[perind] permutation : max|dS_total| after matching by ID = %.3e\n", d))
})

test_that("A4: the decomposition identity holds within tolerance and is reported", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    Q <- cbind(one = f$query, two = stats::setNames(rev(as.numeric(f$query)), f$ids))
    s <- perind_sens(f$fit, query = Q, group = f$group,
                     lambda = f$lambda, backend = "memory")
    bi <- s$by_individual
    expect_identical(anyDuplicated(perind_key(bi)), 0L)

    gate_t <- 0; gate_b <- 0
    for (q in unique(s$by_block$query)) for (b in unique(s$by_block$block)) {
        m <- bi$query == q & bi$block == b
        r <- s$by_block$query == q & s$by_block$block == b
        St <- s$by_block$S_total[r]; Sb <- s$by_block$S_between[r]
        gate_t <- max(gate_t, abs(sum(bi$S_total[m]) - St) / max(1, abs(St)))
        gate_b <- max(gate_b, abs(sum(bi$S_between[m]) - Sb) / max(1, abs(Sb)))
    }
    expect_lt(gate_t, 1e-12)
    expect_lt(gate_b, 1e-12)

    # Reported in its OWN element -- `validity` must stay what it always was, so
    # that a regression bridge can compare it with identical().
    expect_type(s$by_individual_check, "list")
    expect_named(s$by_individual_check, c("sum_rel_total", "sum_rel_between"))
    expect_lt(s$by_individual_check$sum_rel_total, 1e-12)
    expect_lt(s$by_individual_check$sum_rel_between, 1e-12)
    expect_false("by_individual_sum_rel" %in% names(s$validity))
    expect_named(s$validity, c("gap", "rel_gap", "valid", "reason"))
    cat(sprintf("[perind] sum identity: gate_total = %.3e, gate_between = %.3e, reported = %.3e / %.3e\n",
                gate_t, gate_b, s$by_individual_check$sum_rel_total,
                s$by_individual_check$sum_rel_between))

    # Ungrouped results report only the total arm.
    u <- perind_sens(f$fit, query = f$query, lambda = f$lambda, backend = "memory")
    expect_named(u$by_individual_check, "sum_rel_total")
})

test_that("A1: shares are per block, sum to one there, and are NA when undefined", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    bi <- s$by_individual

    d_sh <- 0
    for (b in unique(bi$block)) {
        m <- bi$block == b
        r <- s$by_block$block == b
        # share_ij = s_ij / S_j, with S_j the value published in by_block.
        expect_equal(bi$share_total[m], bi$S_total[m] / s$by_block$S_total[r],
                     tolerance = 1e-14)
        expect_equal(bi$share_between[m], bi$S_between[m] / s$by_block$S_between[r],
                     tolerance = 1e-14)
        d_sh <- max(d_sh, abs(sum(bi$share_total[m]) - 1),
                    abs(sum(bi$share_between[m]) - 1))
    }
    expect_lt(d_sh, 1e-12)
    cat(sprintf("[perind] share sums  : max|sum_i share_ij - 1| over query x block = %.3e\n",
                d_sh))

    expect_true(all(is.finite(bi$share_total)) && all(is.finite(bi$share_between)))
    expect_true(all(bi$share_total >= 0 & bi$share_total <= 1))
    expect_true(all(bi$share_between >= 0 & bi$share_between <= 1))

    # NO THRESHOLDING of the raw contributions. A participant absent from a block
    # contributes ~1e-34 there; that survives into a share as a real, positive,
    # tiny number -- it is neither rounded to zero nor turned into NA, because it
    # is a defined quantity that merely happens to be negligible.
    absent <- bi$share_total[bi$block == "blkC" & bi$id == "id070"]
    expect_true(is.finite(absent) && absent > 0 && absent < 1e-20)

    # The zero/non-finite denominator branch returns NA rather than 0. It is not
    # reachable on this fixture -- every block carries a strictly positive
    # sensitivity -- so what is asserted here is that nothing on the reachable
    # path silently produces the 0 that the branch exists to avoid.
    expect_false(any(bi$share_total == 0))
    expect_false(any(bi$share_between == 0))
})

test_that("the per-individual vector inherits the invariances, and the ridge dependence", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    base <- perind_sens(f$fit, query = f$query, group = f$group,
                        lambda = f$lambda, backend = "memory")
    bv <- perind_total(base)
    bsh <- base$by_individual$share_total

    for (cc in c(0.5, 2, 10)) {
        s <- perind_sens(f$fit, query = cc * f$query, group = f$group,
                         lambda = f$lambda, backend = "memory")
        v <- perind_total(s)
        i <- match(perind_bkey(s$by_individual), perind_bkey(base$by_individual))
        d_ratio <- max(abs(v[names(bv)] / bv - 1))
        d_share <- max(abs(s$by_individual$share_total - bsh[i]))
        expect_lt(d_ratio, 1e-13)
        expect_lt(d_share, 1e-15)
        cat(sprintf("[perind] q -> %sq : max|s ratio - 1| = %.3e, max|dshare| = %.3e\n",
                    format(cc), d_ratio, d_share))
    }

    # Rescaling a whole block changes nothing either: the projector depends only
    # on the block's column space.
    X10 <- f$Xfull; X10$blkA <- 10 * X10$blkA
    tabs10 <- Map(function(X, m) X[m, , drop = FALSE], X10, f$av)
    fit10 <- perind_fit(tabs10, f$L)
    s10 <- perind_sens(fit10, query = f$query, group = f$group,
                       lambda = f$lambda, backend = "memory")
    v10 <- perind_total(s10)
    i <- match(perind_bkey(s10$by_individual), perind_bkey(base$by_individual))
    d_ratio <- max(abs(v10[names(bv)] / bv - 1))
    d_share <- max(abs(s10$by_individual$share_total - bsh[i]))
    expect_lt(d_ratio, 1e-12)
    expect_lt(d_share, 1e-12)
    cat(sprintf("[perind] block x10   : max|s ratio - 1| = %.3e, max|dshare| = %.3e\n",
                d_ratio, d_share))

    # NEGATIVE CONTROL: the layer's ridge is NOT a nuisance parameter here. If a
    # future "helpful" normalisation ever makes this invariant, the quantity
    # being computed is no longer the one documented.
    s2l <- perind_sens(f$fit, query = f$query, group = f$group,
                       lambda = 2 * f$lambda, backend = "memory")
    v2l <- perind_total(s2l)
    d_ridge <- max(abs(v2l[names(bv)] / bv - 1))
    expect_gt(d_ridge, 0.1)
    cat(sprintf("[perind] ridge x2    : max|s ratio - 1| = %.3e (NOT invariant, as documented)\n",
                d_ridge))
})

test_that("the per-individual vector is conditioned on the grouped participant set", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    full <- perind_sens(f$fit, query = f$query, group = f$group,
                        lambda = f$lambda, backend = "memory")
    keep <- f$ids[seq_len(60)]
    part <- perind_sens(f$fit, query = f$query, group = f$group[keep],
                        lambda = f$lambda, backend = "memory")
    expect_identical(nrow(part$by_individual), 60L * 3L)
    ctx <- mgcca:::.mgcca_rel_context(f$fit)
    G60 <- mgcca:::.mgcca_rel_group(f$group[keep], ctx$ids)
    expect_identical(part$by_individual$id[part$by_individual$block == "blkA"],
                     ctx$ids[G60$index + 1L])

    a <- part$by_individual; b <- full$by_individual
    i <- match(perind_bkey(a), perind_bkey(b))
    expect_false(anyNA(i))
    # It CHANGES, and that is the documented conditioning, not a defect.
    d <- max(abs(a$S_total - b$S_total[i]))
    expect_gt(d, 1e-12)
    cat(sprintf("[perind] group 60/70 : max|dS_total| vs full group = %.3e (conditioning)\n", d))
})

test_that("the schema is exact, with and without a grouping", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    g <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    u <- perind_sens(f$fit, query = f$query, lambda = f$lambda, backend = "memory")
    expect_identical(names(g$by_individual),
                     c("query", "id", "block", "S_total", "S_between",
                       "share_total", "share_between"))
    # ABSENT MEANS ABSENT: no NA placeholder columns.
    expect_identical(names(u$by_individual),
                     c("query", "id", "block", "S_total", "share_total"))
    i <- match(perind_bkey(u$by_individual), perind_bkey(g$by_individual))
    d <- max(abs(u$by_individual$S_total - g$by_individual$S_total[i]))
    expect_lt(d, 1e-15)
    cat(sprintf("[perind] no group    : max|dS_total| vs grouped = %.3e\n", d))
})

test_that("adding the per-individual table leaves the existing outputs alone", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    s <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    expect_identical(names(s$overall),
                     c("query", "n_used", "T", "S_total", "S_between", "R_between"))
    expect_identical(names(s$by_block),
                     c("query", "block", "S_total", "S_between",
                       "share_total", "share_between"))

    # overall and by_block still equal the oracle's own block totals.
    o <- perind_oracle(f$fit, f$lambda, f$query, f$group)$tab
    ot <- rowsum(o$S_total, o$block); ob <- rowsum(o$S_between, o$block)
    j <- match(s$by_block$block, rownames(ot))
    expect_equal(s$by_block$S_total, as.numeric(ot[j, 1L]), tolerance = 1e-12)
    expect_equal(s$by_block$S_between, as.numeric(ob[j, 1L]), tolerance = 1e-12)
    expect_equal(s$overall$S_total, sum(ot[, 1L]), tolerance = 1e-12)
    expect_equal(s$overall$S_between, sum(ob[, 1L]), tolerance = 1e-12)
    # by_block's own shares are still across blocks and still sum to one.
    expect_equal(sum(s$by_block$share_total), 1, tolerance = 1e-12)

    # Backends agree on the new table as they already did on the old ones.
    h <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "hdf5")
    i <- match(perind_bkey(h$by_individual), perind_bkey(s$by_individual))
    d_ind <- max(abs(h$by_individual$S_total - s$by_individual$S_total[i]))
    d_blk <- max(abs(h$by_block$S_total - s$by_block$S_total))
    expect_lt(d_ind, 1e-12)
    expect_lt(d_blk, 1e-12)
    cat(sprintf("[perind] memory/hdf5 : max|dS_total| by_individual = %.3e, by_block = %.3e\n",
                d_ind, d_blk))
})

test_that("the per-individual table is persisted with the others", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    # A fresh fit: storing refuses to overwrite, so the shared one is left alone.
    fit <- perind_fit(f$tabs, f$L)
    s <- perind_sens(fit, query = f$query, group = f$group, lambda = f$lambda,
                     backend = "memory", store = TRUE)
    expect_true("by_individual" %in% s$storage$tables)
})

test_that("the individuals plot accounts for the query's WHOLE sensitivity", {
    skip_if_not_installed("BigDataStatMeth")
    skip_if_not_installed("ggplot2")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    Q <- cbind(one = f$query, two = stats::setNames(rev(as.numeric(f$query)), f$ids))
    s <- perind_sens(f$fit, query = Q, group = f$group,
                     lambda = f$lambda, backend = "memory")
    p <- plot(s, type = "individuals", top_n = 10)
    expect_s3_class(p, "ggplot")

    # ⚠️ THE REASON THE "all others" BAR EXISTS. Every bar drawn, top N and
    # remainder together, must add back up to the query's total sensitivity --
    # otherwise the panel shows a top N as though it were the whole thing.
    bl <- ggplot2::ggplot_build(p)
    d <- bl$data[[1]]
    lay <- bl$layout$layout
    tot <- tapply(d$xmax - d$xmin, d$PANEL, sum)
    # Panels are matched by their query LABEL, never by position: facet order and
    # `overall` row order are set independently and do not have to agree.
    qn <- as.character(lay$query[match(names(tot), as.character(lay$PANEL))])
    ref <- s$overall$S_total[match(qn, s$overall$query)]
    d_bar <- max(abs(as.numeric(tot) / ref - 1))
    expect_lt(d_bar, 1e-12)
    # top_n named participants + one remainder bar, stacked over every block.
    expect_identical(nrow(d), 2L * (10L + 1L) * 3L)
    cat(sprintf("[perind] plot bars   : max rel dev of panel totals vs overall = %.3e\n",
                d_bar))

    # With top_n at or above the participant count there is no remainder to draw.
    p2 <- plot(s, type = "individuals", top_n = 70)
    expect_identical(nrow(ggplot2::ggplot_build(p2)$data[[1]]), 2L * 70L * 3L)
})

test_that("the individuals plot works ungrouped, and refuses what it cannot draw", {
    skip_if_not_installed("BigDataStatMeth")
    skip_if_not_installed("ggplot2")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    g <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    u <- perind_sens(f$fit, query = f$query, lambda = f$lambda, backend = "memory")
    expect_s3_class(plot(u, type = "individuals"), "ggplot")

    # Grouped results carry the group in the axis label; ungrouped ones cannot.
    lab <- function(p) as.character(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels())
    expect_true(any(grepl("^id\\d+ \\(g[12]\\)$", lab(plot(g, type = "individuals")))))
    expect_true(any(grepl("^id\\d+$", lab(plot(u, type = "individuals")))))
    # The remainder bar names no group: it is not one participant.
    expect_true(any(grepl("^all others \\(55\\)$", lab(plot(g, type = "individuals")))))

    # An object from a version without the decomposition is refused, not drawn.
    old <- g; old$by_individual <- NULL
    expect_error(plot(old, type = "individuals"), "older version")

    for (bad in list(0, -1, c(5, 6), "five", NA_real_, Inf))
        expect_error(plot(g, type = "individuals", top_n = bad), "top_n")
})

test_that("print and summary report the table without claiming more than it is", {
    skip_if_not_installed("BigDataStatMeth")
    f <- perind_fixture()
    on.exit(BigDataStatMeth::hdf5_close_all(), add = TRUE)
    g <- perind_sens(f$fit, query = f$query, group = f$group,
                     lambda = f$lambda, backend = "memory")
    u <- perind_sens(f$fit, query = f$query, lambda = f$lambda, backend = "memory")
    expect_output(print(g), "individuals: 70 per query")
    expect_output(print(u), "individuals: 70 per query")
    expect_output(summary(g), "Per-individual sensitivity contributions \\(exploratory\\)")
    expect_output(summary(u), "Per-individual sensitivity contributions \\(exploratory\\)")
    expect_output(summary(g), "share_of_query")
    expect_output(summary(g), "MIXTO calibration evidence")
    expect_output(summary(g), "no individual-level calibration has been established")

    # The displayed share is against the QUERY total, so the five shown add to
    # less than one; the table's own share_total is a within-block share.
    out <- utils::capture.output(summary(g))
    expect_false(any(grepl("patient", out, ignore.case = TRUE)))

    # ⚠️ An object from a version that had no per-individual table prints as it
    # did. Only `by_individual` is dropped here, NOT `by_individual_check`:
    # `$` partial-matches, so a lookup of `x$by_individual` on such an object
    # silently returns the check element instead. Dropping both would hide that.
    old <- g; old$by_individual <- NULL
    expect_output(print(old), "mgcca sensitivity")
    expect_output(summary(old), "Per-block decomposition")
    out <- utils::capture.output(summary(old))
    expect_false(any(grepl("Per-individual", out, fixed = TRUE)))
})
