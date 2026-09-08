## Pre-fit availability and design audit.
##
## This module reads an availability design -- who is present in which block --
## and nothing else. It never calls the estimator, it never reads a data value,
## and no quantity of a fit can move because of it. It runs BEFORE a fit, which
## is the point: the questions it answers ("is this pair of blocks supported by
## anybody?", "does this stratum have enough weight behind it?") are settled by
## the design, not by the numbers that come out of it.
##
## What it does NOT do, deliberately:
##
##   * It does not conclude that anything is or is not identified. Whether the
##     working model is identified depends on structural conditions this audit
##     does not evaluate -- hence `$structural_checks`, always "not_evaluated"
##     here, present so that a consumer never has to guess whether the check was
##     run and passed or never run at all.
##   * It does not refuse. Every well-formed design gets an object back, with
##     the values and the advisory codes side by side, and the reader decides. A
##     diagnostic that refuses has made the decision for them, and it makes it
##     on a line somebody drew by hand.
##   * It does not weight anything, correct anything, or feed anything back into
##     an estimator. Declared weights are read to describe their own
##     concentration; they are never applied to data here.
##
## The design lines (`ess_min`, `availability_min`, `max_weight`) are DECLARED
## lines carrying margin. They are not estimated, they are not thresholds at
## which something fails, and the value that raised a code is always reported
## next to the code so the reader can see how close to the line it really is.
##
## Input polymorphism is a cohesion requirement, not a convenience: auditing the
## very object the fit will consume is what makes the audit a statement about
## the actual analysis. The availability table is derived through the package's
## own readers -- the same union-of-IDs rule the estimator uses in both of its
## modes -- and only presence/absence is read; no block data is materialised.


## ---- frozen vocabulary ------------------------------------------------------
## Pair codes partition the J*(J-1)/2 pairs; advisory codes are additive and a
## pair or a cell can carry none of them. No IDENTIFIED_* code exists.
.mgcca_av_pair_codes <- c("DIRECT_SUPPORT", "MODEL_TRANSFER_ELIGIBLE",
                          "UNSUPPORTED", "STRUCTURAL_ZERO")
.mgcca_av_advisory_codes <- c("LOW_SUPPORT", "SINGLE_BRIDGE",
                              "SEVERE_REPRESENTATION_WARNING")

.mgcca_av_default_lines <- list(ess_min = 30, availability_min = 0.05,
                                max_weight = 20)

.mgcca_av_store_group <- "availability_audit"


## ---- small graph primitives -------------------------------------------------
## The support graph has J nodes, and J is the number of omics blocks: two to a
## couple of dozen. Everything below is deliberately the textbook algorithm in
## base R rather than a dependency -- at this size the constant factors are
## irrelevant and the code is checkable by hand against the tests.

## Connected components, iterative depth-first. `adj` is J x J logical with a
## FALSE diagonal.
.mgcca_av_components <- function(adj) {
    J <- nrow(adj)
    comp <- integer(J)
    if (!J) return(comp)
    cur <- 0L
    for (v in seq_len(J)) {
        if (comp[v] > 0L) next
        cur <- cur + 1L
        stack <- v
        while (length(stack)) {
            u <- stack[1L]
            stack <- stack[-1L]
            if (comp[u] > 0L) next
            comp[u] <- cur
            nb <- which(adj[u, ])
            stack <- c(stack, nb[comp[nb] == 0L])
        }
    }
    comp
}

## Max-flow between two nodes with every undirected edge at unit capacity
## (Edmonds-Karp: BFS augmenting paths). By Menger's theorem this is both the
## number of edge-disjoint j--k routes and the size of the smallest set of edges
## whose removal separates them, which is what kappa reports.
.mgcca_av_maxflow <- function(adj, s, t) {
    J <- nrow(adj)
    cap <- matrix(0L, J, J)
    cap[adj] <- 1L
    flow <- 0L
    repeat {
        prev <- integer(J)
        prev[s] <- s
        queue <- s
        while (length(queue) && prev[t] == 0L) {
            u <- queue[1L]
            queue <- queue[-1L]
            nxt <- which(cap[u, ] > 0L & prev == 0L)
            if (length(nxt)) {
                prev[nxt] <- u
                queue <- c(queue, nxt)
            }
        }
        if (prev[t] == 0L) break
        v <- t
        while (v != s) {
            u <- prev[v]
            cap[u, v] <- cap[u, v] - 1L
            cap[v, u] <- cap[v, u] + 1L
            v <- u
        }
        flow <- flow + 1L
    }
    flow
}

## kappa_jk for every pair. The diagonal is NA: the edge connectivity of a node
## with itself is not a quantity, and a 0 there would be read as "disconnected".
.mgcca_av_edge_connectivity <- function(adj) {
    J <- nrow(adj)
    K <- matrix(NA_integer_, J, J, dimnames = dimnames(adj))
    if (J < 2L) return(K)
    for (j in seq_len(J - 1L)) for (k in (j + 1L):J) {
        v <- .mgcca_av_maxflow(adj, j, k)
        K[j, k] <- K[k, j] <- v
    }
    K
}

## A node is an articulation block when deleting it leaves the remaining blocks
## in more pieces than they were already in.
.mgcca_av_articulation <- function(adj) {
    J <- nrow(adj)
    base_comp <- .mgcca_av_components(adj)
    out <- vapply(seq_len(J), function(v) {
        keep <- setdiff(seq_len(J), v)
        if (length(keep) < 2L) return(FALSE)
        sub <- adj[keep, keep, drop = FALSE]
        max(.mgcca_av_components(sub)) > length(unique(base_comp[keep]))
    }, logical(1))
    stats::setNames(out, rownames(adj))
}

## Widest route from one block: a max-min Dijkstra. `w` is the edge-support
## matrix (0 where there is no edge). best[v] is the largest B such that some
## j -> v route has every edge supported by at least B; prev[] reconstructs it.
## Ties are broken towards the lowest block index, so the result is a
## deterministic function of the column order of the availability table.
.mgcca_av_widest_from <- function(w, s) {
    J <- nrow(w)
    best <- rep(-Inf, J)
    prev <- rep(NA_integer_, J)
    done <- rep(FALSE, J)
    best[s] <- Inf
    repeat {
        cand <- which(!done & best > -Inf)
        if (!length(cand)) break
        u <- cand[which.max(best[cand])]
        done[u] <- TRUE
        for (v in seq_len(J)) {
            if (done[v] || w[u, v] <= 0) next
            val <- min(best[u], w[u, v])
            if (val > best[v]) {
                best[v] <- val
                prev[v] <- u
            }
        }
    }
    list(best = best, prev = prev)
}

## Kish's effective sample size, (sum w)^2 / sum(w^2). An empty cell, or one
## whose declared weights are all zero, has ESS 0 -- not NA: the cell really
## does carry no effective observation, and an NA would make every consumer
## special-case it.
.mgcca_av_kish <- function(w) {
    if (!length(w)) return(0)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) return(0)
    sw^2 / sum(w^2)
}

## Paste the names of whichever flags fired, "" when none did.
.mgcca_av_paste_codes <- function(...) {
    fl <- list(...)
    nm <- names(fl)
    if (!length(fl[[1L]])) return(character(0))
    vapply(seq_along(fl[[1L]]), function(i) {
        hit <- nm[vapply(fl, function(f) isTRUE(f[i]), logical(1))]
        if (length(hit)) paste(hit, collapse = ";") else ""
    }, character(1))
}


## ---- input: the availability table -----------------------------------------
## An availability table handed over directly (design audits before any data
## exists, and the form the tests are written against).
.mgcca_av_matrix <- function(availability) {
    if (is.data.frame(availability)) {
        ok <- vapply(availability, function(z)
            is.logical(z) || is.numeric(z), logical(1))
        if (!all(ok))
            stop("the availability table must hold logical or 0/1 columns; ",
                 "column(s) ", paste(names(availability)[!ok], collapse = ", "),
                 " are neither.", call. = FALSE)
        rn <- rownames(availability)
        availability <- as.matrix(availability)
        rownames(availability) <- rn
    }
    if (!is.matrix(availability))
        stop("the availability table must be an n x J logical (or 0/1) matrix ",
             "or data frame, one row per individual and one column per block.",
             call. = FALSE)
    if (!nrow(availability) || ncol(availability) < 2L)
        stop("the availability table needs at least one individual and at ",
             "least two blocks; with one block there is no pair to report on.",
             call. = FALSE)
    if (anyNA(availability))
        stop("the availability table must not contain NA: a missing entry is ",
             "not the same statement as `this individual is absent from this ",
             "block`, and the audit will not guess which was meant.",
             call. = FALSE)
    if (is.numeric(availability) && !all(availability %in% c(0, 1)))
        stop("the availability table must be logical or 0/1; it holds other ",
             "values.", call. = FALSE)
    blocks <- colnames(availability)
    if (is.null(blocks)) blocks <- paste0("block", seq_len(ncol(availability)))
    if (anyNA(blocks) || any(!nzchar(blocks)))
        stop("every column of the availability table needs a block name.",
             call. = FALSE)
    if (anyDuplicated(blocks))
        stop("block names must be unique; duplicated: ",
             paste(unique(blocks[duplicated(blocks)]), collapse = ", "),
             call. = FALSE)
    matrix(availability != 0, nrow(availability), ncol(availability),
           dimnames = list(rownames(availability), blocks))
}

## Assemble the availability table from per-block individual IDs, using the
## SAME union rule the estimator uses (mgcca_memory.R / mgcca_getK): first
## appearance when every table has the same number of rows, sorted otherwise.
## If the two ever diverged, the audit would be describing a different design
## from the one that gets fitted.
.mgcca_av_assemble <- function(ids) {
    if (length(ids) < 2L)
        stop("mgcca needs at least two tables; the audit needs at least two ",
             "blocks to have a pair to report on.", call. = FALSE)
    for (nm in names(ids))
        if (is.null(ids[[nm]]) || !length(ids[[nm]]))
            stop("table '", nm, "' has no rownames (individual IDs required)",
                 call. = FALSE)
    ns <- vapply(ids, length, integer(1))
    rn <- if (max(ns) == min(ns)) Reduce(union, ids) else sort(Reduce(union, ids))
    A <- matrix(FALSE, length(rn), length(ids),
                dimnames = list(rn, names(ids)))
    for (j in seq_along(ids)) A[, j] <- rn %in% ids[[j]]
    A
}

## In-memory blocks: exactly the input types mgcca() takes with filename = NULL.
## Only the individual IDs are touched -- for a list of matrices nothing is
## copied at all, and for the assay-based classes only the sample names are
## read. Deliberately NOT .mgcca_as_tables(), which materialises every block.
.mgcca_av_ids_memory <- function(x) {
    if (is.list(x) && !methods::is(x, "MultiAssayExperiment")) {
        nms <- names(x)
        if (is.null(nms) || any(!nzchar(nms)))
            nms <- paste0("table", seq_along(x))
        ids <- lapply(x, rownames)
        names(ids) <- nms
        return(ids)
    }
    if (methods::is(x, "MultiAssayExperiment")) {
        if (!requireNamespace("MultiAssayExperiment", quietly = TRUE))
            stop("package 'MultiAssayExperiment' is required for this input",
                 call. = FALSE)
        ex <- MultiAssayExperiment::experiments(x)
        ## the importer stores assays transposed (individuals x variables), so
        ## the individual IDs are the assay's COLUMN names.
        ids <- lapply(seq_along(ex), function(j)
            colnames(SummarizedExperiment::assay(ex[[j]])))
        names(ids) <- names(ex)
        return(ids)
    }
    if (methods::is(x, "SummarizedExperiment"))
        return(list(assay = colnames(SummarizedExperiment::assay(x))))
    if (methods::is(x, "ExpressionSet"))
        return(list(exprs = colnames(Biobase::exprs(x))))
    stop("unsupported input type for mgcca_audit: ", class(x)[1L],
         ". Give it an availability matrix, the block object mgcca() would ",
         "fit, or the path/descriptor of an HDF5 file holding the blocks.",
         call. = FALSE)
}

## HDF5-backed blocks: the dimnames dataset only. rownames() on the R6 reader
## goes through BigDataStatMeth's own dimnames reader, so not one value of any
## block is read.
.mgcca_av_ids_hdf5 <- function(filename, group, datasets) {
    on.exit(suppressMessages(try(BigDataStatMeth::hdf5_close_all(),
                                 silent = TRUE)), add = TRUE)
    out <- vector("list", length(datasets))
    names(out) <- datasets
    for (k in seq_along(datasets)) {
        hm <- BigDataStatMeth::hdf5_matrix(filename,
                                           paste0(group, "/", datasets[k]))
        rn <- tryCatch(rownames(hm), error = function(e) NULL)
        if (!is.null(hm) && is.function(hm$close)) try(hm$close(), silent = TRUE)
        if (is.null(rn) || !length(rn))
            stop("dataset '", datasets[k], "' in ", filename,
                 " has no rownames (individual IDs required)", call. = FALSE)
        out[[k]] <- as.character(rn)
    }
    out
}

## The dispatcher. Returns the availability table plus where it came from, and
## the HDF5 descriptor when there is one (save_hdf5 needs it).
.mgcca_av_input <- function(data) {
    if (is.matrix(data) || is.data.frame(data))
        return(list(availability = .mgcca_av_matrix(data),
                    source = "availability table", hdf5 = NULL))

    ## HDF5-backed: a path, or the descriptor mgcca_import_hdf5() and
    ## mgcca(collect = FALSE) return (which is how a non-default group and a
    ## dataset subset get expressed without a second argument here).
    desc <- NULL
    if (is.character(data) && length(data) == 1L)
        desc <- list(filename = data, group = "MGCCA_IN", datasets = NULL)
    else if (is.list(data) && is.character(data[["filename"]]) &&
             length(data[["filename"]]) == 1L)
        ## `[[` on purpose: `$` partially matches, so a block named
        ## `filename_something` would be mistaken for a descriptor.
        desc <- list(filename = data[["filename"]],
                     group = if (is.null(data[["group"]])) "MGCCA_IN"
                             else data[["group"]],
                     datasets = data[["datasets"]])
    if (!is.null(desc)) {
        if (!file.exists(desc$filename))
            stop("the HDF5 file does not exist: ", desc$filename, call. = FALSE)
        ## read-only branch of the package's own importer: for a path it lists
        ## the datasets and returns, writing nothing.
        d <- mgcca_import_hdf5(desc$filename, filename = desc$filename,
                               group = desc$group, datasets = desc$datasets)
        suppressMessages(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE))
        return(list(availability = .mgcca_av_matrix(
                        .mgcca_av_assemble(.mgcca_av_ids_hdf5(
                            d$filename, d$group, d$datasets))),
                    source = "HDF5", hdf5 = d))
    }

    list(availability = .mgcca_av_matrix(
             .mgcca_av_assemble(.mgcca_av_ids_memory(data))),
         source = "in-memory blocks", hdf5 = NULL)
}

.mgcca_av_lines <- function(design_lines) {
    if (is.null(design_lines)) design_lines <- list()
    if (!is.list(design_lines))
        stop("`design_lines` must be a named list.", call. = FALSE)
    unknown <- setdiff(names(design_lines), names(.mgcca_av_default_lines))
    if (length(unknown))
        stop("unknown design line(s): ", paste(unknown, collapse = ", "),
             ". Available: ",
             paste(names(.mgcca_av_default_lines), collapse = ", "),
             call. = FALSE)
    out <- utils::modifyList(.mgcca_av_default_lines, design_lines)
    for (nm in names(out)) {
        v <- out[[nm]]
        if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0)
            stop("design line `", nm,
                 "` must be a single finite non-negative number.", call. = FALSE)
        ## always double, so a line declared as 30L and one declared as 30 give
        ## the same object -- and survive a round trip through HDF5 unchanged.
        out[[nm]] <- as.numeric(v)
    }
    out
}


## ---- persistence -------------------------------------------------------------
## The audit travels with the data: `save_hdf5 = TRUE` writes the COMPLETE audit
## into an `availability_audit` group of the same file, and mgcca_audit_load()
## brings it back whole.
##
## HOW IT IS ENCODED, and why it is not a blob. Everything is either a double
## dataset or a character dimnames vector, so any HDF5 browser can read the
## certificate without R:
##
##   * numeric and logical columns -> double datasets, column names carried as
##     the dataset's colnames;
##   * block names, individual IDs and stratum labels -> character DIMNAMES of
##     the datasets they name, written and read by BigDataStatMeth's own
##     dimnames machinery (hdf5_create_matrix / hdf5_matrix), which is the
##     package's character support;
##   * status codes -> integers, with the frozen vocabulary stored beside them
##     as the `pair_codes` and `advisory_codes` legend datasets (the code is the
##     row index, the label is the row name), so an integer is never orphaned
##     from its meaning. Advisory codes are additive, so they are stored as a
##     bit sum against the `advisory_codes` legend;
##   * transfer routes -> the `paths` dataset, one row per pair, the route as a
##     sequence of integer block indices, zero-padded;
##   * the manifest (version, call, source, the rank-L sentence) -> group
##     attributes, through the same attribute writer FINAL_RESULTS and
##     RELIABILITY/ already use.
##
## rhdf5 is not used anywhere in this module, deliberately and by rule.
##
## HOW IT IS RESTORED. The audit is a deterministic function of five stored
## inputs -- the availability table, L, the declared strata, the declared
## weights and the design lines -- so the loader rebuilds the object from those
## and then CHECKS the rebuild against the stored derived tables, warning if the
## two disagree. That way the file is a full, inspectable certificate AND the
## restored object cannot silently drift from it if the package changes.
##
## Rewriting is idempotent: the dataset set and the attribute set are fixed, so
## every one of them is overwritten on every save and no piece of an older audit
## can survive beside a newer one.

.mgcca_av_advisory_bits <- function(codes) {
    bits <- stats::setNames(2^(seq_along(.mgcca_av_advisory_codes) - 1L),
                            .mgcca_av_advisory_codes)
    vapply(codes, function(s) {
        if (!nzchar(s)) return(0)
        sum(bits[strsplit(s, ";", fixed = TRUE)[[1L]]])
    }, numeric(1), USE.NAMES = FALSE)
}

## Write one double matrix, closing the handle the R6 reader leaves open.
.mgcca_av_write <- function(filename, path, m) {
    storage.mode(m) <- "double"
    hm <- BigDataStatMeth::hdf5_create_matrix(filename, path, data = m,
                                              dtype = "double", overwrite = TRUE)
    if (!is.null(hm) && is.function(hm$close)) try(hm$close(), silent = TRUE)
    suppressMessages(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE))
    invisible(path)
}

.mgcca_av_read <- function(filename, path) {
    hm <- BigDataStatMeth::hdf5_matrix(filename, path)
    m <- as.matrix(hm)
    dn <- try(dimnames(hm), silent = TRUE)
    if (is.function(hm$close)) try(hm$close(), silent = TRUE)
    suppressMessages(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE))
    if (!inherits(dn, "try-error") && length(dn) == 2L) dimnames(m) <- dn
    ## The audit never produces a NaN -- every stored number is a count, a
    ## proportion, an ESS or a missing value -- so any NaN that comes back is a
    ## missing value that lost its payload on the way through, and treating it
    ## as one keeps NA_real_ round-tripping as NA_real_ rather than NaN.
    m[is.nan(m)] <- NA_real_
    m
}

## The encoded form of the derived tables: what goes in the file, and what the
## loader re-derives and compares against.
.mgcca_av_encode <- function(x) {
    b <- x$blocks
    slev <- x$strata_levels
    p <- x$provenance$pairs
    P <- nrow(p)

    bfrom <- bto <- rep(0, P)
    for (e in seq_len(P)) {
        pe <- x$provenance$paths[[e]]
        if (is.null(pe)) next
        wk <- which.min(pe$edges$s_e)
        bfrom[e] <- match(pe$edges$from[wk], b)
        bto[e] <- match(pe$edges$to[wk], b)
    }

    pairs <- cbind(
        block1_index = match(p$block1, b), block2_index = match(p$block2, b),
        code_index = match(p$code, .mgcca_av_pair_codes),
        n_jk = p$n_jk, pi_hat_jk = p$pi_hat_jk,
        n_jkc_min = p$n_jkc_min, pi_jkc_min = p$pi_jkc_min,
        direct = as.numeric(p$direct),
        same_component = as.numeric(p$same_component),
        s_e = p$s_e, s_e_ess = p$s_e_ess, kappa = p$kappa,
        path_length = p$path_length, bottleneck = p$bottleneck,
        bottleneck_from_index = bfrom, bottleneck_to_index = bto,
        single_bridge = as.numeric(p$single_bridge),
        alternative_paths = as.numeric(p$alternative_paths),
        L = p$L, n_jk_ge_L = as.numeric(p$n_jk_ge_L),
        bottleneck_ge_L = as.numeric(p$bottleneck_ge_L),
        advisory_bits = .mgcca_av_advisory_bits(x$status$pairs$advisory))

    ## routes as integer block-index sequences, zero-padded
    steps <- max(1L, max(c(0L, p$path_length + 1L), na.rm = TRUE))
    paths <- matrix(0, P, steps,
                    dimnames = list(NULL, paste0("step", seq_len(steps))))
    for (e in seq_len(P)) {
        pe <- x$provenance$paths[[e]]
        if (is.null(pe)) next
        idx <- match(pe$path, b)
        paths[e, seq_along(idx)] <- idx
    }

    r <- x$representation
    representation <- cbind(
        stratum_index = match(r$stratum, slev), block_index = match(r$block, b),
        n_stratum = r$n_stratum, n_obs = r$n_obs, proportion = r$proportion,
        ess = r$ess, max_weight = r$max_weight,
        weight_concentration = r$weight_concentration,
        below_availability_min = as.numeric(r$below_availability_min),
        below_ess_min = as.numeric(r$below_ess_min),
        above_max_weight = as.numeric(r$above_max_weight),
        advisory_bits = .mgcca_av_advisory_bits(r$advisory))

    ab <- x$availability$by_block
    availability_by_block <- matrix(
        c(ab$n, ab$n_j, ab$pi_hat_j), nrow(ab), 3L,
        dimnames = list(b, c("n", "n_j", "pi_hat_j")))

    as_ <- x$availability$by_stratum
    availability_by_stratum <- cbind(
        stratum_index = match(as_$stratum, slev),
        block_index = match(as_$block, b),
        n_stratum = as_$n_stratum, n_jc = as_$n_jc, pi_hat_jc = as_$pi_hat_jc)

    ps <- x$availability$pairs_by_stratum
    if (is.null(ps)) {
        ## no strata declared: the single-stratum joint counts are the pair
        ## counts themselves, written so the dataset set never varies.
        pairs_by_stratum <- cbind(
            stratum_index = rep(1, P),
            block1_index = match(p$block1, b), block2_index = match(p$block2, b),
            n_stratum = rep(x$n, P), n_jkc = p$n_jk, pi_hat_jkc = p$pi_hat_jk)
    } else {
        pairs_by_stratum <- cbind(
            stratum_index = match(ps$stratum, slev),
            block1_index = match(ps$block1, b),
            block2_index = match(ps$block2, b),
            n_stratum = ps$n_stratum, n_jkc = ps$n_jkc,
            pi_hat_jkc = ps$pi_hat_jkc)
    }

    graph <- matrix(c(as.numeric(x$graph$component),
                      as.numeric(x$graph$articulation), ab$n_j),
                    length(b), 3L,
                    dimnames = list(b, c("component", "articulation", "n_j")))

    design <- matrix(c(x$n, x$J, x$L, x$design_lines$ess_min,
                       x$design_lines$availability_min,
                       x$design_lines$max_weight,
                       as.numeric(x$strata_declared),
                       as.numeric(x$weights_declared),
                       x$graph$n_components, as.numeric(x$graph$connected),
                       as.numeric(x$graph$complete)), 1L, 11L,
                     dimnames = list(NULL,
                                     c("n", "J", "L", "ess_min",
                                       "availability_min", "max_weight",
                                       "strata_declared", "weights_declared",
                                       "n_components", "connected", "complete")))

    list(pairs = pairs, paths = paths, representation = representation,
         availability_by_block = availability_by_block,
         availability_by_stratum = availability_by_stratum,
         pairs_by_stratum = pairs_by_stratum, graph = graph, design = design)
}

.mgcca_av_store <- function(x, filename, strata_index, w,
                            group = .mgcca_av_store_group) {
    on.exit(suppressMessages(try(BigDataStatMeth::hdf5_close_all(),
                                 silent = TRUE)), add = TRUE)
    enc <- .mgcca_av_encode(x)
    b <- x$blocks
    slev <- x$strata_levels

    tab <- matrix(as.numeric(x$table), x$n, x$J,
                  dimnames = list(x$individuals, b))
    kap <- x$graph$kappa
    storage.mode(kap) <- "double"

    fixed <- list(
        table = tab,
        strata = matrix(as.numeric(strata_index), x$n, 1L,
                        dimnames = list(NULL, "stratum_index")),
        weights = matrix(as.numeric(w), x$n, 1L,
                         dimnames = list(NULL, "weight")),
        strata_levels = matrix(seq_along(slev), length(slev), 1L,
                               dimnames = list(slev, "stratum_index")),
        pair_codes = matrix(seq_along(.mgcca_av_pair_codes),
                            length(.mgcca_av_pair_codes), 1L,
                            dimnames = list(.mgcca_av_pair_codes, "code_index")),
        advisory_codes = matrix(2^(seq_along(.mgcca_av_advisory_codes) - 1L),
                                length(.mgcca_av_advisory_codes), 1L,
                                dimnames = list(.mgcca_av_advisory_codes,
                                                "bit_value")),
        adjacency = matrix(as.numeric(x$graph$adjacency), x$J, x$J,
                           dimnames = list(b, b)),
        kappa = kap)
    all_tables <- c(fixed, enc)
    for (nm in names(all_tables))
        .mgcca_av_write(filename, paste0(group, "/", nm), all_tables[[nm]])

    att <- list(
        mgcca_version = as.character(utils::packageVersion("mgcca")),
        mgcca_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        analysis = "mgcca_audit",
        structural_checks = x$structural_checks,
        source = x$source,
        call = paste(deparse(x$call), collapse = " "),
        rank_assumption = x$provenance$rank_assumption,
        hdf5_group = as.character(x$hdf5$group),
        hdf5_datasets = as.character(x$hdf5$datasets),
        blocks = b,
        has_individuals = as.integer(!is.null(x$individuals)),
        tables = names(all_tables),
        n = as.integer(x$n), J = as.integer(x$J), L = as.integer(x$L),
        strata_declared = as.integer(x$strata_declared),
        weights_declared = as.integer(x$weights_declared),
        ess_min = as.numeric(x$design_lines$ess_min),
        availability_min = as.numeric(x$design_lines$availability_min),
        max_weight = as.numeric(x$design_lines$max_weight))
    mgcca_write_attrs_rcpp(filename, group, "", att)
    suppressMessages(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE))
    list(filename = filename, group = group, tables = names(all_tables))
}


#' Read a stored availability audit back from HDF5
#'
#' @description Restores, whole, the audit that \code{\link{mgcca_audit}} wrote
#'   into a file with \code{save_hdf5 = TRUE}. The result is the same
#'   \code{"mgcca_audit"} object, with the same tables, codes and routes.
#'
#' @details The group holds every table the audit produced, as double datasets
#'   with their names carried in the dimnames and the frozen vocabulary stored
#'   beside the integer codes as a legend, plus the audit's five inputs -- the
#'   availability table, the declared rank \code{L}, the declared strata, the
#'   declared weights and the design lines.
#'
#'   The object is reassembled \strong{component by component from what is
#'   stored}: nothing is deserialised and nothing is recomputed, so what comes
#'   back is what went in, and a fault in the encoding shows up as a difference
#'   rather than being papered over by a recomputation. Because the inputs are
#'   stored too, the audit can also be re-derived from them and compared against
#'   what was stored; \code{check_inputs} does exactly that and reports a
#'   disagreement -- which would mean the file was written by a different
#'   version of the audit -- instead of leaving it silent.
#'
#'   Nothing here reads a block of data, and no estimator is called. This is a
#'   provenance record: the certificate travelling with the file it describes.
#'
#' @param file Path to the HDF5 file holding the audit.
#' @param group Group the audit was written to. Default
#'   \code{"availability_audit"}, which is where \code{mgcca_audit()} puts it.
#' @param check_inputs If \code{TRUE} (default), also re-derive the audit from
#'   the inputs stored beside it and compare, so a file written by a different
#'   version of the audit is reported rather than trusted in silence. The object
#'   returned is the stored one either way.
#' @return An object of class \code{"mgcca_audit"}; see \code{\link{mgcca_audit}}.
#' @seealso \code{\link{mgcca_audit}}
#' @examples
#' # a small HDF5-backed design: three blocks over a union of 60 individuals
#' ids <- sprintf("ind%03d", 1:60)
#' blocks <- list(
#'     methylation = matrix(rnorm(30 * 4), 30, 4,
#'                          dimnames = list(ids[1:30], paste0("m", 1:4))),
#'     expression  = matrix(rnorm(60 * 3), 60, 3,
#'                          dimnames = list(ids, paste0("e", 1:3))),
#'     proteins    = matrix(rnorm(30 * 2), 30, 2,
#'                          dimnames = list(ids[31:60], paste0("p", 1:2))))
#' h5 <- tempfile(fileext = ".h5")
#' desc <- mgcca_import_hdf5(blocks, filename = h5, overwriteFile = TRUE)
#'
#' a <- mgcca_audit(desc, L = 2, save_hdf5 = TRUE)
#' b <- mgcca_audit_load(h5)
#' identical(b$provenance$pairs, a$provenance$pairs)
#' unlink(h5)
#' @export
mgcca_audit_load <- function(file, group = .mgcca_av_store_group,
                             check_inputs = TRUE) {
    if (!is.character(file) || length(file) != 1L)
        stop("`file` must be the path to a single HDF5 file.", call. = FALSE)
    if (!file.exists(file))
        stop("the HDF5 file does not exist: ", file, call. = FALSE)
    on.exit(suppressMessages(try(BigDataStatMeth::hdf5_close_all(),
                                 silent = TRUE)), add = TRUE)

    att <- tryCatch(mgcca_read_attrs_rcpp(file, group, ""),
                    error = function(e) NULL)
    suppressMessages(try(BigDataStatMeth::hdf5_close_all(), silent = TRUE))
    if (is.null(att) || !identical(as.character(att$analysis), "mgcca_audit"))
        stop("no mgcca audit is stored in '", group, "' of ", file,
             ". Write one with mgcca_audit(..., save_hdf5 = TRUE).",
             call. = FALSE)

    ds <- .mgcca_av_expected_tables()
    missing_ds <- setdiff(ds, as.character(att$tables))
    if (length(missing_ds))
        stop("the audit stored in '", group, "' of ", file,
             " is incomplete: ", paste(missing_ds, collapse = ", "),
             " missing. Write it again with ",
             "mgcca_audit(..., save_hdf5 = TRUE).", call. = FALSE)
    stored <- lapply(stats::setNames(ds, ds), function(nm)
        .mgcca_av_read(file, paste0(group, "/", nm)))

    out <- .mgcca_av_decode(stored, att)
    out$hdf5 <- list(filename = file, group = as.character(att$hdf5_group),
                     datasets = as.character(att$hdf5_datasets))
    out$saved <- list(filename = file, group = group,
                      tables = as.character(att$tables))

    if (isTRUE(check_inputs)) {
        enc <- .mgcca_av_encode(.mgcca_av_from_stored_inputs(stored, att))
        drift <- character(0)
        for (nm in names(enc))
            if (!isTRUE(all.equal(unname(stored[[nm]]), unname(enc[[nm]]),
                                  tolerance = 1e-12)))
                drift <- c(drift, nm)
        if (length(drift))
            warning("the audit stored in '", group, "' of ", file,
                    " does not match what its own inputs imply today (",
                    paste(drift, collapse = ", "),
                    "). The file was most likely written by a different ",
                    "version of mgcca_audit(). What is returned is what was ",
                    "stored, not what those inputs imply now.", call. = FALSE)
    }
    out
}

## The fixed dataset set. Fixed on purpose: every save writes all of them, so a
## rewrite cannot leave a piece of an older audit behind, and a load can say
## exactly what is missing.
.mgcca_av_expected_tables <- function()
    c("table", "strata", "weights", "strata_levels", "pair_codes",
      "advisory_codes", "adjacency", "kappa", "pairs", "paths",
      "representation", "availability_by_block", "availability_by_stratum",
      "pairs_by_stratum", "graph", "design")

## Re-derive the audit from the inputs stored beside it. Only ever used to
## CHECK what was stored; it is never what mgcca_audit_load() returns.
.mgcca_av_from_stored_inputs <- function(stored, att) {
    blocks <- as.character(att$blocks)
    tab <- stored$table
    ind <- if (isTRUE(as.integer(att$has_individuals) == 1L))
        rownames(tab) else NULL
    A <- matrix(tab != 0, nrow(tab), ncol(tab), dimnames = list(ind, blocks))
    strata <- NULL
    if (isTRUE(as.integer(att$strata_declared) == 1L)) {
        slev <- rownames(stored$strata_levels)
        strata <- factor(slev[as.integer(stored$strata[, 1L])], levels = slev)
    }
    weights <- if (isTRUE(as.integer(att$weights_declared) == 1L))
        as.numeric(stored$weights[, 1L]) else NULL
    mgcca_audit(A, L = as.integer(att$L), strata = strata, weights = weights,
                design_lines = list(ess_min = as.numeric(att$ess_min),
                                    availability_min = as.numeric(att$availability_min),
                                    max_weight = as.numeric(att$max_weight)),
                save_hdf5 = FALSE)
}

## Split an advisory bit sum back into its frozen code names.
.mgcca_av_bits_to_codes <- function(bits) {
    vapply(as.numeric(bits), function(b) {
        hit <- .mgcca_av_advisory_codes[
            bitwAnd(as.integer(b), 2L^(seq_along(.mgcca_av_advisory_codes) - 1L)) > 0L]
        if (length(hit)) paste(hit, collapse = ";") else ""
    }, character(1), USE.NAMES = FALSE)
}

## Assemble the object from the stored datasets and the stored manifest, one
## component at a time. Nothing is deserialised and nothing is recomputed.
.mgcca_av_decode <- function(stored, att) {
    blocks <- as.character(att$blocks)
    slev <- rownames(stored$strata_levels)
    J <- length(blocks)
    n <- as.integer(att$n)
    L <- as.integer(att$L)
    strata_declared <- isTRUE(as.integer(att$strata_declared) == 1L)
    weights_declared <- isTRUE(as.integer(att$weights_declared) == 1L)
    design <- stored$design
    dcol <- function(nm) unname(design[1L, nm])

    ind <- if (isTRUE(as.integer(att$has_individuals) == 1L))
        rownames(stored$table) else NULL
    A <- matrix(stored$table != 0, n, J, dimnames = list(ind, blocks))

    lines <- list(ess_min = as.numeric(att$ess_min),
                  availability_min = as.numeric(att$availability_min),
                  max_weight = as.numeric(att$max_weight))

    ## --- availability ----------------------------------------------------
    ab <- stored$availability_by_block
    av_block <- data.frame(block = blocks,
                           n = as.integer(ab[, "n"]),
                           n_j = as.integer(ab[, "n_j"]),
                           pi_hat_j = as.numeric(ab[, "pi_hat_j"]),
                           stringsAsFactors = FALSE)

    bs <- stored$availability_by_stratum
    av_stratum <- data.frame(
        stratum = slev[as.integer(bs[, "stratum_index"])],
        block = blocks[as.integer(bs[, "block_index"])],
        n_stratum = as.integer(bs[, "n_stratum"]),
        n_jc = as.integer(bs[, "n_jc"]),
        pi_hat_jc = as.numeric(bs[, "pi_hat_jc"]),
        stringsAsFactors = FALSE)

    p <- stored$pairs
    b1 <- blocks[as.integer(p[, "block1_index"])]
    b2 <- blocks[as.integer(p[, "block2_index"])]
    av_pairs <- data.frame(block1 = b1, block2 = b2,
                           n_jk = as.integer(p[, "n_jk"]),
                           pi_hat_jk = as.numeric(p[, "pi_hat_jk"]),
                           n_jkc_min = as.integer(p[, "n_jkc_min"]),
                           pi_jkc_min = as.numeric(p[, "pi_jkc_min"]),
                           stringsAsFactors = FALSE)

    av_pairs_stratum <- NULL
    if (strata_declared) {
        ps <- stored$pairs_by_stratum
        av_pairs_stratum <- data.frame(
            stratum = slev[as.integer(ps[, "stratum_index"])],
            block1 = blocks[as.integer(ps[, "block1_index"])],
            block2 = blocks[as.integer(ps[, "block2_index"])],
            n_stratum = as.integer(ps[, "n_stratum"]),
            n_jkc = as.integer(ps[, "n_jkc"]),
            pi_hat_jkc = as.numeric(ps[, "pi_hat_jkc"]),
            stringsAsFactors = FALSE)
    }

    ## --- graph ------------------------------------------------------------
    adj <- matrix(stored$adjacency != 0, J, J, dimnames = list(blocks, blocks))
    kappa <- matrix(as.integer(stored$kappa), J, J,
                    dimnames = list(blocks, blocks))
    gnode <- stored$graph
    comp <- stats::setNames(as.integer(gnode[, "component"]), blocks)
    artic <- stats::setNames(as.logical(gnode[, "articulation"]), blocks)

    direct <- as.logical(p[, "direct"])
    s_e <- as.numeric(p[, "s_e"])
    s_e_ess <- as.numeric(p[, "s_e_ess"])
    kap_pair <- as.integer(p[, "kappa"])
    edges <- data.frame(block1 = b1[direct], block2 = b2[direct],
                        n_jk = av_pairs$n_jk[direct], s_e = s_e[direct],
                        s_e_ess = s_e_ess[direct], kappa = kap_pair[direct],
                        is_bridge = kap_pair[direct] == 1L,
                        stringsAsFactors = FALSE)
    rownames(edges) <- NULL

    ## --- provenance -------------------------------------------------------
    path_len <- as.integer(p[, "path_length"])
    pmat <- stored$paths
    route_of <- function(e) {
        if (is.na(path_len[e])) return(NULL)
        as.integer(pmat[e, seq_len(path_len[e] + 1L)])
    }
    path_chr <- vapply(seq_along(path_len), function(e) {
        r <- route_of(e)
        if (is.null(r)) NA_character_ else paste(blocks[r], collapse = " -> ")
    }, character(1))
    bfrom <- as.integer(p[, "bottleneck_from_index"])
    bto <- as.integer(p[, "bottleneck_to_index"])
    bneck_edge <- ifelse(bfrom == 0L, NA_character_,
                         paste(blocks[pmax(bfrom, 1L)],
                               blocks[pmax(bto, 1L)], sep = "--"))

    prov_pairs <- data.frame(
        block1 = b1, block2 = b2,
        code = .mgcca_av_pair_codes[as.integer(p[, "code_index"])],
        n_jk = av_pairs$n_jk, pi_hat_jk = av_pairs$pi_hat_jk,
        n_jkc_min = av_pairs$n_jkc_min, pi_jkc_min = av_pairs$pi_jkc_min,
        direct = direct, same_component = as.logical(p[, "same_component"]),
        s_e = s_e, s_e_ess = s_e_ess, kappa = kap_pair,
        path = path_chr, path_length = path_len,
        bottleneck = as.numeric(p[, "bottleneck"]),
        bottleneck_edge = bneck_edge,
        single_bridge = as.logical(p[, "single_bridge"]),
        alternative_paths = as.logical(p[, "alternative_paths"]),
        L = as.integer(p[, "L"]),
        n_jk_ge_L = as.logical(p[, "n_jk_ge_L"]),
        bottleneck_ge_L = as.logical(p[, "bottleneck_ge_L"]),
        stringsAsFactors = FALSE)

    edge_row <- function(from, to)
        which((b1 == from & b2 == to) | (b1 == to & b2 == from))
    paths <- vector("list", nrow(prov_pairs))
    names(paths) <- paste(b1, b2, sep = "--")
    for (e in seq_along(paths)) {
        r <- route_of(e)
        if (is.null(r)) next
        from <- blocks[r[-length(r)]]
        to <- blocks[r[-1L]]
        m <- vapply(seq_along(from), function(i) edge_row(from[i], to[i]),
                    integer(1))
        paths[[e]] <- list(
            path = blocks[r],
            edges = data.frame(from = from, to = to, s_e = s_e[m],
                               s_e_ess = s_e_ess[m], kappa = kap_pair[m],
                               is_bridge = kap_pair[m] == 1L,
                               stringsAsFactors = FALSE),
            bottleneck = prov_pairs$bottleneck[e], kappa = kap_pair[e],
            alternative_paths = prov_pairs$alternative_paths[e])
    }

    ## --- representation and status ---------------------------------------
    r <- stored$representation
    rep_adv <- .mgcca_av_bits_to_codes(r[, "advisory_bits"])
    rep_tab <- data.frame(
        stratum = slev[as.integer(r[, "stratum_index"])],
        block = blocks[as.integer(r[, "block_index"])],
        n_stratum = as.integer(r[, "n_stratum"]),
        n_obs = as.integer(r[, "n_obs"]),
        proportion = as.numeric(r[, "proportion"]),
        ess = as.numeric(r[, "ess"]),
        max_weight = as.numeric(r[, "max_weight"]),
        weight_concentration = as.numeric(r[, "weight_concentration"]),
        below_availability_min = as.logical(r[, "below_availability_min"]),
        below_ess_min = as.logical(r[, "below_ess_min"]),
        above_max_weight = as.logical(r[, "above_max_weight"]),
        advisory = rep_adv, stringsAsFactors = FALSE)

    pair_bits <- as.integer(p[, "advisory_bits"])
    status_pairs <- data.frame(
        block1 = b1, block2 = b2, code = prov_pairs$code,
        advisory = .mgcca_av_bits_to_codes(pair_bits),
        low_support = bitwAnd(pair_bits, 1L) > 0L,
        single_bridge = bitwAnd(pair_bits, 2L) > 0L,
        stringsAsFactors = FALSE)
    rep_bits <- as.integer(r[, "advisory_bits"])
    status_strata <- data.frame(
        stratum = rep_tab$stratum, block = rep_tab$block,
        advisory = rep_adv,
        low_support = bitwAnd(rep_bits, 1L) > 0L,
        severe_representation_warning = bitwAnd(rep_bits, 4L) > 0L,
        stringsAsFactors = FALSE)

    out <- list(
        call = tryCatch(str2lang(as.character(att$call)),
                        error = function(e) NULL),
        source = as.character(att$source), hdf5 = NULL,
        table = A, individuals = ind,
        n = n, J = J, blocks = blocks, L = L,
        strata_declared = strata_declared, strata_levels = slev,
        weights_declared = weights_declared,
        design_lines = lines,
        availability = list(by_block = av_block, by_stratum = av_stratum,
                            pairs = av_pairs,
                            pairs_by_stratum = av_pairs_stratum),
        graph = list(adjacency = adj, edges = edges, component = comp,
                     n_components = as.integer(dcol("n_components")),
                     connected = dcol("connected") == 1,
                     complete = dcol("complete") == 1, kappa = kappa,
                     bridges = edges[edges$is_bridge, , drop = FALSE],
                     articulation = artic),
        provenance = list(pairs = prov_pairs, paths = paths,
                          rank_assumption = as.character(att$rank_assumption)),
        representation = rep_tab,
        status = list(pairs = status_pairs, strata = status_strata),
        structural_checks = as.character(att$structural_checks),
        saved = NULL)
    class(out) <- "mgcca_audit"
    out
}



## ---- the audit --------------------------------------------------------------

#' Pre-fit availability and design audit
#'
#' @description Audits the availability design of a multi-block study
#'   \emph{before} anything is fitted: how many individuals each block and each
#'   declared stratum actually holds, which pairs of blocks share enough
#'   individuals to be supported directly, which pairs are reachable only
#'   through other blocks, how fragile those routes are, and how concentrated
#'   any declared availability weights are. It reads presence and absence, the
#'   declared strata and the declared weights, and nothing else: it never calls
#'   \code{\link{mgcca}}, never reads a data value, and cannot change any
#'   quantity a fit produces.
#'
#' @details
#' \strong{It audits the object the fit will consume.} Besides a ready-made
#' availability table, \code{data} accepts the same inputs \code{\link{mgcca}}
#' takes in either mode -- an in-memory list of blocks (or a
#' \code{MultiAssayExperiment}, \code{SummarizedExperiment} or
#' \code{ExpressionSet}), or an HDF5 file holding them -- and derives the
#' availability table itself through the package's own readers, with the same
#' union-of-IDs rule the estimator uses. Only the individual IDs are read; no
#' block data is materialised, on either route. The union order is the order the
#' rows of the audit are in, and \code{$individuals} reports it, so declared
#' \code{strata} and \code{weights} can be lined up with it.
#'
#' \strong{Support graph.} Blocks \eqn{j} and \eqn{k} are joined by an edge when
#' the realised joint availability \eqn{\hat\pi_{jk,c}} reaches
#' \code{design_lines$availability_min} in \emph{every} declared stratum
#' \eqn{c}. With no strata declared there is one stratum, the whole table. A
#' pair that clears the line in the sample as a whole but has an empty joint
#' cell in one declared stratum gets no edge, and the per-stratum table shows
#' which stratum emptied it.
#'
#' \strong{Pair codes} (frozen vocabulary, and a partition of the pairs):
#' \describe{
#'   \item{\code{DIRECT_SUPPORT}}{the pair has its own edge: individuals
#'     observed in both blocks carry it.}
#'   \item{\code{MODEL_TRANSFER_ELIGIBLE}}{no edge, but a route of edges runs
#'     from \eqn{j} to \eqn{k}. The pair is then identified conditional on the
#'     declared rank-\code{L} model and its required structural conditions,
#'     which are not empirically checked by this audit. \code{$provenance}
#'     names the route, the support of every edge on it, and its bottleneck.}
#'   \item{\code{UNSUPPORTED}}{no edge and no route, but the pair does have
#'     jointly observed individuals -- too few, or absent from a declared
#'     stratum.}
#'   \item{\code{STRUCTURAL_ZERO}}{no edge, no route, and not one individual
#'     anywhere is observed in both blocks. The joint cell is empty by design.}
#' }
#'
#' \strong{Advisory codes} are additive and carry no verdict:
#' \code{LOW_SUPPORT} (the support the pair leans on sits below the declared
#' \code{ess_min} line), \code{SINGLE_BRIDGE} (\eqn{\kappa_{jk}=1}: one edge
#' carries the whole connection) and \code{SEVERE_REPRESENTATION_WARNING} (a
#' stratum-block cell crosses one of the declared representation lines). The
#' value that raised a code is always reported beside it.
#'
#' \strong{Design-support fragility.} \eqn{\kappa_{jk}} is the number of
#' edge-disjoint routes between two blocks, equivalently the smallest number of
#' edges whose loss would separate them. It measures the \emph{topological
#' redundancy of the design's support} -- how many independent ways the design
#' has of connecting two blocks -- and it is silent about how well anything is
#' estimated. A bridge is an edge with \eqn{\kappa = 1}: lose those individuals
#' and the two sides stop being connected at all.
#'
#' \strong{Support counts are descriptive.} The default edge-support metric is
#' the raw joint count \eqn{n_{jk}}. When \code{weights} are declared, an
#' ESS-based support is reported \emph{in addition}, in its own column; it never
#' replaces the raw count. Comparing either against \code{L} is descriptive
#' support and nothing more: \eqn{n_{jk} \ge L} does not prove that the
#' corresponding cross-block quantity has numerical rank \code{L}.
#'
#' \strong{Positivity.} A finite entry in the availability table \emph{evidences}
#' positivity for that cell -- somebody was observed there. It does not prove
#' it: an empty cell may be empty by chance or impossible by design, and the
#' table cannot tell the two apart.
#'
#' \strong{Design lines are declared, with margin.} \code{ess_min},
#' \code{availability_min} and \code{max_weight} are conventional lines set by
#' the analyst, not thresholds at which identification fails and not quantities
#' estimated from the data. Everything the audit reports is continuous: the
#' object carries the value (an ESS of 23.7, say) together with any code it
#' raised, and no design makes the audit refuse to return one.
#'
#' \strong{Structural checks are not run.} \code{$structural_checks} is always
#' \code{"not_evaluated"}. Empirical rank and conditioning checks on the blocks
#' and on the transfer routes are a different, later question; the field is
#' present so that a consumer can never confuse "checked and fine" with "not
#' checked".
#'
#' Among the pre-fit diagnostics we surveyed for multi-block methods with
#' missing individuals, we found none reporting the transfer route and its
#' bottleneck alongside the representation of the declared strata.
#'
#' @param data The design to audit. Either an \eqn{n \times J} logical (or 0/1)
#'   availability matrix or data frame (columns named for the blocks,
#'   \code{TRUE} where the individual is present), or the very object
#'   \code{\link{mgcca}} would fit: a named list of matrices, a
#'   \code{MultiAssayExperiment}, a \code{SummarizedExperiment}, an
#'   \code{ExpressionSet}, a path to an HDF5 file holding the blocks, or the
#'   descriptor list (\code{filename}, \code{group}, \code{datasets}) returned by
#'   \code{\link{mgcca_import_hdf5}}. A bare path carries no declared order or
#'   subset of the blocks, so the file's own listing decides both -- as it does
#'   for \code{\link{mgcca}}; pass the descriptor to fix them.
#' @param L Declared rank of the working model, an integer \eqn{\ge 1}. It is
#'   recorded and compared against the support counts descriptively; nothing is
#'   fitted at rank \code{L} here.
#' @param strata Optional declared strata, a factor or vector of length \eqn{n},
#'   in the row order of the audit (\code{$individuals}). Availability,
#'   representation and the edge rule are then reported within every stratum.
#'   Default \code{NULL}: one stratum, the whole table.
#' @param weights Optional user-declared availability weights, a non-negative
#'   numeric vector of length \eqn{n}, in the same row order. They are described
#'   (ESS, maximum, concentration) and used for the additional ESS-based support
#'   column. They are never applied to data, and never replace the raw counts.
#' @param design_lines Declared lines, with margin: \code{ess_min} (minimum
#'   effective count, default 30, applied to stratum ESS and to pair support),
#'   \code{availability_min} (minimum realised availability, default 0.05) and
#'   \code{max_weight} (largest declared weight tolerated in a cell, default
#'   20). A partial list is merged with the defaults.
#' @param save_hdf5 If \code{TRUE} and \code{data} is HDF5-backed, write the
#'   complete audit into an \code{availability_audit} group of that same file,
#'   so the certificate travels with the data, and read it back whole with
#'   \code{\link{mgcca_audit_load}}. The encoding is structured and inspectable
#'   rather than a blob: numeric tables as datasets, block names, individual IDs
#'   and stratum labels as their dimnames, status codes as integers with the
#'   frozen vocabulary stored beside them as a legend, transfer routes as
#'   sequences of integer block indices. Rewriting is idempotent. The returned
#'   object is the same either way, and is \code{saveRDS}-able.
#'
#' @return An object of class \code{"mgcca_audit"}, a list with
#'   \describe{
#'     \item{\code{$availability}}{\code{by_block} and \code{by_stratum}
#'       (\eqn{n_{jc}} and \eqn{\hat\pi_{jc}}), \code{pairs} (\eqn{n_{jk}}) and
#'       \code{pairs_by_stratum} (\eqn{n_{jk,c}}; \code{NULL} when no strata are
#'       declared).}
#'     \item{\code{$graph}}{\code{adjacency}, \code{edges}, \code{component},
#'       \code{n_components}, \code{connected}, \code{complete}, \code{kappa},
#'       \code{bridges} and \code{articulation} -- design-support fragility and
#'       topological redundancy, not identification robustness.}
#'     \item{\code{$provenance}}{one row per pair with its code, supports,
#'       \eqn{\kappa}, best route, bottleneck and fragility; \code{$paths} holds
#'       the route of every connected pair edge by edge;
#'       \code{$rank_assumption} states what transfer is conditional on.}
#'     \item{\code{$representation}}{per declared stratum and block: observed
#'       proportion, Kish ESS (only when \code{weights} are declared), maximum
#'       weight, weight concentration and status.}
#'     \item{\code{$status}}{machine-readable codes, \code{$pairs} and
#'       \code{$strata}.}
#'     \item{\code{$structural_checks}}{always \code{"not_evaluated"}.}
#'   }
#'   plus \code{$table} (the availability table it worked from),
#'   \code{$individuals}, \code{$blocks}, \code{$L}, \code{$design_lines},
#'   \code{$source} and \code{$saved}.
#'
#' @seealso \code{\link{mgcca}} for the estimator this audit deliberately does
#'   not touch, and \code{\link{mgcca_audit_load}} for reading a stored audit
#'   back out of the file it was written into.
#'
#' @examples
#' # A design in which one block is only reachable through another: 40
#' # individuals carry methylation and expression, 60 carry expression and
#' # proteins, and nobody carries methylation and proteins together.
#' av <- cbind(methylation = rep(c(TRUE, FALSE), c(40, 60)),
#'             expression  = rep(TRUE, 100),
#'             proteins    = rep(c(FALSE, TRUE), c(40, 60)))
#' a <- mgcca_audit(av, L = 2)
#' a
#'
#' # methylation--proteins has no shared individual, but a route through
#' # expression, and the audit reports the route and its bottleneck.
#' a$provenance$pairs[, c("block1", "block2", "code", "path", "bottleneck")]
#'
#' # Declared strata and declared weights: the audit adds the per-stratum
#' # representation of each block and the concentration of the weights.
#' st <- rep(c("cohort A", "cohort B"), each = 50)
#' w  <- rep(c(1, 3), each = 50)
#' summary(mgcca_audit(av, L = 2, strata = st, weights = w))
#'
#' # The same audit taken from the blocks themselves: the availability table is
#' # derived from the individual IDs, with mgcca's own union rule.
#' ids <- sprintf("ind%03d", 1:100)
#' blocks <- list(
#'     methylation = matrix(0, 40, 3, dimnames = list(ids[1:40], paste0("m", 1:3))),
#'     expression  = matrix(0, 100, 2, dimnames = list(ids, paste0("e", 1:2))),
#'     proteins    = matrix(0, 60, 2, dimnames = list(ids[41:100], paste0("p", 1:2))))
#' mgcca_audit(blocks, L = 2)$graph$component
#' @export
mgcca_audit <- function(data, L, strata = NULL, weights = NULL,
                        design_lines = list(ess_min = 30,
                                            availability_min = 0.05,
                                            max_weight = 20),
                        save_hdf5 = FALSE) {
    cl <- match.call()
    src <- .mgcca_av_input(data)
    A <- src$availability
    n <- nrow(A)
    J <- ncol(A)
    blocks <- colnames(A)
    individuals <- rownames(A)
    lines <- .mgcca_av_lines(design_lines)

    if (!is.logical(save_hdf5) || length(save_hdf5) != 1L || is.na(save_hdf5))
        stop("`save_hdf5` must be TRUE or FALSE.", call. = FALSE)
    if (isTRUE(save_hdf5) && is.null(src$hdf5))
        stop("`save_hdf5 = TRUE` needs an HDF5-backed design to write into: ",
             "the audit was built from ", src$source,
             ", and there is no file to put it beside. Audit the HDF5 file (or ",
             "its descriptor) instead, or keep the object with saveRDS().",
             call. = FALSE)

    if (missing(L) || !is.numeric(L) || length(L) != 1L || !is.finite(L) ||
        L < 1 || L != round(L))
        stop("`L` must be the declared rank of the working model: a single ",
             "integer >= 1.", call. = FALSE)
    L <- as.integer(L)

    strata_declared <- !is.null(strata)
    if (strata_declared) {
        if (length(strata) != n)
            stop("`strata` must have one entry per individual, in the row ",
                 "order of the audit: length ", length(strata), " for ", n,
                 " individuals.", call. = FALSE)
        if (anyNA(strata))
            stop("`strata` must not contain NA; a declared stratum is a ",
                 "statement about the design, and there is no default for it.",
                 call. = FALSE)
        sf <- droplevels(factor(strata))
    } else {
        sf <- factor(rep("(all)", n))
    }
    slev <- levels(sf)

    weights_declared <- !is.null(weights)
    if (weights_declared) {
        if (!is.numeric(weights) || length(weights) != n)
            stop("`weights` must be a numeric vector with one entry per ",
                 "individual, in the row order of the audit: length ",
                 length(weights), " for ", n, " individuals.", call. = FALSE)
        if (anyNA(weights) || any(!is.finite(weights)))
            stop("`weights` must be finite.", call. = FALSE)
        if (any(weights < 0))
            stop("`weights` must be non-negative.", call. = FALSE)
        w <- as.numeric(weights)
    } else {
        w <- rep(1, n)
    }

    ## ---- (A) availability summaries -------------------------------------
    n_j <- vapply(seq_len(J), function(j) sum(A[, j]), integer(1))
    av_block <- data.frame(block = blocks, n = n, n_j = n_j,
                           pi_hat_j = n_j / n, stringsAsFactors = FALSE)

    n_c <- as.integer(table(sf))
    gr <- expand.grid(block = blocks, stratum = slev, stringsAsFactors = FALSE)
    gr <- gr[, c("stratum", "block")]
    gr$n_stratum <- n_c[match(gr$stratum, slev)]
    gr$n_jc <- vapply(seq_len(nrow(gr)), function(r)
        sum(A[sf == gr$stratum[r], gr$block[r]]), integer(1))
    gr$pi_hat_jc <- ifelse(gr$n_stratum > 0L, gr$n_jc / gr$n_stratum, 0)
    av_stratum <- gr[order(match(gr$stratum, slev),
                           match(gr$block, blocks)), , drop = FALSE]
    rownames(av_stratum) <- NULL

    ij <- utils::combn(J, 2L)
    n_jk <- vapply(seq_len(ncol(ij)), function(e)
        sum(A[, ij[1L, e]] & A[, ij[2L, e]]), integer(1))
    av_pairs <- data.frame(block1 = blocks[ij[1L, ]], block2 = blocks[ij[2L, ]],
                           n_jk = n_jk, pi_hat_jk = n_jk / n,
                           stringsAsFactors = FALSE)

    ## Joint counts per declared stratum. The edge rule reads the minimum of
    ## pi_hat_jkc over the declared strata, so these are computed either way and
    ## only EXPOSED as a table when strata were actually declared.
    pjkc <- matrix(0, ncol(ij), length(slev), dimnames = list(NULL, slev))
    njkc <- matrix(0L, ncol(ij), length(slev), dimnames = list(NULL, slev))
    for (ci in seq_along(slev)) {
        keep <- sf == slev[ci]
        nc <- sum(keep)
        for (e in seq_len(ncol(ij))) {
            cnt <- sum(A[keep, ij[1L, e]] & A[keep, ij[2L, e]])
            njkc[e, ci] <- cnt
            pjkc[e, ci] <- if (nc > 0L) cnt / nc else 0
        }
    }
    av_pairs$n_jkc_min <- as.integer(apply(njkc, 1L, min))
    av_pairs$pi_jkc_min <- apply(pjkc, 1L, min)

    av_pairs_stratum <- NULL
    if (strata_declared) {
        av_pairs_stratum <- do.call(rbind, lapply(seq_along(slev), function(ci)
            data.frame(stratum = slev[ci],
                       block1 = blocks[ij[1L, ]], block2 = blocks[ij[2L, ]],
                       n_stratum = n_c[ci], n_jkc = njkc[, ci],
                       pi_hat_jkc = pjkc[, ci], stringsAsFactors = FALSE)))
        rownames(av_pairs_stratum) <- NULL
    }

    ## ---- (B) structural support graph -----------------------------------
    is_edge <- av_pairs$n_jk > 0L &
        av_pairs$pi_jkc_min >= lines$availability_min
    adj <- matrix(FALSE, J, J, dimnames = list(blocks, blocks))
    adj[cbind(ij[1L, ][is_edge], ij[2L, ][is_edge])] <- TRUE
    adj[cbind(ij[2L, ][is_edge], ij[1L, ][is_edge])] <- TRUE

    comp <- stats::setNames(.mgcca_av_components(adj), blocks)
    kappa <- .mgcca_av_edge_connectivity(adj)
    artic <- .mgcca_av_articulation(adj)

    ## Edge support: the raw joint count is the default metric. The ESS-based
    ## support is an ADDITIONAL column when weights are declared; it never takes
    ## the raw count's place, here or anywhere downstream.
    s_e <- as.numeric(av_pairs$n_jk)
    s_e_ess <- rep(NA_real_, length(s_e))
    if (weights_declared)
        s_e_ess <- vapply(seq_len(ncol(ij)), function(e)
            .mgcca_av_kish(w[A[, ij[1L, e]] & A[, ij[2L, e]]]), numeric(1))

    kap_pair <- kappa[cbind(ij[1L, ], ij[2L, ])]
    edges <- data.frame(block1 = av_pairs$block1[is_edge],
                        block2 = av_pairs$block2[is_edge],
                        n_jk = av_pairs$n_jk[is_edge],
                        s_e = s_e[is_edge], s_e_ess = s_e_ess[is_edge],
                        kappa = kap_pair[is_edge],
                        is_bridge = kap_pair[is_edge] == 1L,
                        stringsAsFactors = FALSE)
    rownames(edges) <- NULL

    ## ---- (C) provenance --------------------------------------------------
    W <- matrix(0, J, J, dimnames = list(blocks, blocks))
    W[cbind(ij[1L, ][is_edge], ij[2L, ][is_edge])] <- s_e[is_edge]
    W[cbind(ij[2L, ][is_edge], ij[1L, ][is_edge])] <- s_e[is_edge]
    wide <- lapply(seq_len(J), function(j) .mgcca_av_widest_from(W, j))

    ## unname() on purpose: `comp` is named by block, and a named logical
    ## column silently becomes the data frame's ROW NAMES when it is the only
    ## named column -- which happens at J = 2, where the pair's names do not
    ## collide. The pair tables would then be shaped differently at J = 2 than
    ## at J = 3.
    same_comp <- unname(comp[ij[1L, ]] == comp[ij[2L, ]])
    direct <- is_edge
    code <- ifelse(direct, "DIRECT_SUPPORT",
            ifelse(same_comp, "MODEL_TRANSFER_ELIGIBLE",
            ifelse(av_pairs$n_jk > 0L, "UNSUPPORTED", "STRUCTURAL_ZERO")))

    bottleneck <- rep(NA_real_, ncol(ij))
    path_chr <- rep(NA_character_, ncol(ij))
    path_len <- rep(NA_integer_, ncol(ij))
    bneck_edge <- rep(NA_character_, ncol(ij))
    paths <- vector("list", ncol(ij))
    names(paths) <- paste(av_pairs$block1, av_pairs$block2, sep = "--")

    for (e in seq_len(ncol(ij))) {
        j <- ij[1L, e]
        k <- ij[2L, e]
        if (!same_comp[e]) next
        wf <- wide[[j]]
        route <- k
        v <- k
        while (v != j) {
            v <- wf$prev[v]
            if (is.na(v)) { route <- NULL; break }
            route <- c(v, route)
        }
        if (is.null(route)) next
        bottleneck[e] <- wf$best[k]
        path_chr[e] <- paste(blocks[route], collapse = " -> ")
        path_len[e] <- length(route) - 1L
        ed <- data.frame(from = blocks[route[-length(route)]],
                         to = blocks[route[-1L]],
                         s_e = vapply(seq_len(length(route) - 1L), function(i)
                             W[route[i], route[i + 1L]], numeric(1)),
                         stringsAsFactors = FALSE)
        ed$s_e_ess <- NA_real_
        ed$kappa <- NA_integer_
        ed$is_bridge <- NA
        for (i in seq_len(nrow(ed))) {
            m <- which((av_pairs$block1 == ed$from[i] &
                            av_pairs$block2 == ed$to[i]) |
                           (av_pairs$block1 == ed$to[i] &
                                av_pairs$block2 == ed$from[i]))
            ed$s_e_ess[i] <- s_e_ess[m]
            ed$kappa[i] <- kap_pair[m]
            ed$is_bridge[i] <- kap_pair[m] == 1L
        }
        wk <- which.min(ed$s_e)
        bneck_edge[e] <- paste(ed$from[wk], ed$to[wk], sep = "--")
        paths[[e]] <- list(path = blocks[route], edges = ed,
                           bottleneck = wf$best[k], kappa = kap_pair[e],
                           alternative_paths = kap_pair[e] >= 2L)
    }

    support_used <- ifelse(same_comp, bottleneck, as.numeric(av_pairs$n_jk))
    support_used_ess <- ifelse(direct, s_e_ess, NA_real_)
    low_support <- support_used < lines$ess_min |
        (weights_declared & !is.na(support_used_ess) &
             support_used_ess < lines$ess_min)
    low_support[is.na(low_support)] <- FALSE
    single_bridge <- same_comp & kap_pair == 1L

    prov_pairs <- data.frame(
        block1 = av_pairs$block1, block2 = av_pairs$block2,
        code = code,
        n_jk = av_pairs$n_jk, pi_hat_jk = av_pairs$pi_hat_jk,
        n_jkc_min = av_pairs$n_jkc_min, pi_jkc_min = av_pairs$pi_jkc_min,
        direct = direct, same_component = same_comp,
        s_e = s_e, s_e_ess = s_e_ess, kappa = kap_pair,
        path = path_chr, path_length = path_len,
        bottleneck = bottleneck, bottleneck_edge = bneck_edge,
        single_bridge = single_bridge,
        alternative_paths = same_comp & kap_pair >= 2L,
        L = L,
        n_jk_ge_L = av_pairs$n_jk >= L,
        bottleneck_ge_L = bottleneck >= L,
        stringsAsFactors = FALSE)

    rank_assumption <- sprintf(
        paste("Pairs coded MODEL_TRANSFER_ELIGIBLE are identified conditional",
              "on the declared rank-%d model and its required structural",
              "conditions, which are not empirically checked by this audit.",
              "n_jk >= L and bottleneck >= L are descriptive support only:",
              "they do not prove numerical rank %d."), L, L)

    ## ---- (D) representation ---------------------------------------------
    rep_tab <- av_stratum
    rep_tab$n_obs <- rep_tab$n_jc
    rep_tab$proportion <- rep_tab$pi_hat_jc
    rep_tab$n_jc <- NULL
    rep_tab$pi_hat_jc <- NULL
    rep_tab$ess <- NA_real_
    rep_tab$max_weight <- NA_real_
    rep_tab$weight_concentration <- NA_real_
    if (weights_declared) {
        for (r in seq_len(nrow(rep_tab))) {
            ww <- w[sf == rep_tab$stratum[r] & A[, rep_tab$block[r]]]
            ess <- .mgcca_av_kish(ww)
            rep_tab$ess[r] <- ess
            rep_tab$max_weight[r] <- if (length(ww)) max(ww) else NA_real_
            ## Concentration: the share of the cell's nominal size lost to the
            ## spread of the declared weights. 0 = every declared weight equal.
            rep_tab$weight_concentration[r] <-
                if (length(ww)) 1 - ess / length(ww) else NA_real_
        }
    }

    fire_avail <- rep_tab$proportion < lines$availability_min
    fire_ess <- weights_declared & rep_tab$ess < lines$ess_min
    fire_maxw <- weights_declared & !is.na(rep_tab$max_weight) &
        rep_tab$max_weight > lines$max_weight
    fire_ess[is.na(fire_ess)] <- FALSE
    fire_maxw[is.na(fire_maxw)] <- FALSE
    severe <- fire_avail | fire_ess | fire_maxw
    low_cell <- !severe & rep_tab$n_obs < lines$ess_min

    rep_tab$below_availability_min <- fire_avail
    rep_tab$below_ess_min <- fire_ess
    rep_tab$above_max_weight <- fire_maxw
    rep_tab$advisory <- .mgcca_av_paste_codes(
        LOW_SUPPORT = low_cell, SEVERE_REPRESENTATION_WARNING = severe)
    rownames(rep_tab) <- NULL

    ## ---- (E) machine-readable status ------------------------------------
    status_pairs <- data.frame(
        block1 = prov_pairs$block1, block2 = prov_pairs$block2,
        code = prov_pairs$code,
        advisory = .mgcca_av_paste_codes(LOW_SUPPORT = low_support,
                                         SINGLE_BRIDGE = single_bridge),
        low_support = low_support, single_bridge = single_bridge,
        stringsAsFactors = FALSE)
    status_strata <- data.frame(
        stratum = rep_tab$stratum, block = rep_tab$block,
        advisory = rep_tab$advisory,
        low_support = low_cell, severe_representation_warning = severe,
        stringsAsFactors = FALSE)

    out <- list(
        call = cl, source = src$source, hdf5 = src$hdf5,
        table = A, individuals = individuals,
        n = n, J = J, blocks = blocks, L = L,
        strata_declared = strata_declared, strata_levels = slev,
        weights_declared = weights_declared,
        design_lines = lines,
        availability = list(by_block = av_block, by_stratum = av_stratum,
                            pairs = av_pairs,
                            pairs_by_stratum = av_pairs_stratum),
        graph = list(adjacency = adj, edges = edges, component = comp,
                     n_components = max(comp), connected = max(comp) == 1L,
                     complete = all(is_edge), kappa = kappa,
                     bridges = edges[edges$is_bridge, , drop = FALSE],
                     articulation = artic),
        provenance = list(pairs = prov_pairs, paths = paths,
                          rank_assumption = rank_assumption),
        representation = rep_tab,
        status = list(pairs = status_pairs, strata = status_strata),
        structural_checks = "not_evaluated",
        saved = NULL)
    class(out) <- "mgcca_audit"

    if (isTRUE(save_hdf5))
        out$saved <- .mgcca_av_store(out, src$hdf5$filename,
                                     strata_index = as.integer(sf), w = w)
    out
}


## ---- methods ----------------------------------------------------------------

#' @export
print.mgcca_audit <- function(x, ...) {
    cat("mgcca availability audit (pre-fit; the estimator is not called)\n")
    cat(sprintf("  source      : %s\n", x$source))
    cat(sprintf("  individuals : %d ; blocks : %d ; declared rank L = %d\n",
                x$n, x$J, x$L))
    cat(sprintf("  strata      : %s\n",
                if (x$strata_declared)
                    sprintf("%d declared (%s)", length(x$strata_levels),
                            paste(utils::head(x$strata_levels, 4),
                                  collapse = ", "))
                else "none declared -- one stratum, the whole table"))
    cat(sprintf("  weights     : %s\n",
                if (x$weights_declared)
                    "declared (ESS reported in addition to the raw counts)"
                else "none declared -- ESS is not reported"))
    cat(sprintf("  design lines: ess_min = %s ; availability_min = %s ; max_weight = %s\n",
                format(x$design_lines$ess_min),
                format(x$design_lines$availability_min),
                format(x$design_lines$max_weight)))
    if (!is.null(x$saved))
        cat(sprintf("  stored in   : %s [%s]\n", x$saved$filename,
                    x$saved$group))
    cat("\nSupport graph\n")
    cat(sprintf("  components %d%s ; edges %d of %d pairs%s\n",
                x$graph$n_components,
                if (x$graph$connected) " (connected)" else "",
                nrow(x$graph$edges), nrow(x$provenance$pairs),
                if (x$graph$complete) " (complete)" else ""))
    nb <- nrow(x$graph$bridges)
    na <- sum(x$graph$articulation)
    cat(sprintf("  design-support fragility: %d bridge%s, %d articulation block%s\n",
                nb, if (nb == 1L) "" else "s", na, if (na == 1L) "" else "s"))
    cat("\nPair codes\n")
    tb <- table(factor(x$provenance$pairs$code, levels = .mgcca_av_pair_codes))
    for (nm in names(tb)) cat(sprintf("  %-24s %d\n", nm, tb[[nm]]))
    adv <- x$status$pairs[nzchar(x$status$pairs$advisory), , drop = FALSE]
    if (nrow(adv)) {
        cat("\nPairs carrying an advisory\n")
        print(format(adv[, c("block1", "block2", "code", "advisory")]),
              row.names = FALSE)
    }
    sadv <- x$status$strata[nzchar(x$status$strata$advisory), , drop = FALSE]
    if (nrow(sadv)) {
        cat("\nStratum-block cells carrying an advisory\n")
        keep <- c("stratum", "block", "n_obs", "proportion", "ess",
                  "max_weight", "advisory")
        d <- merge(sadv[, c("stratum", "block")], x$representation,
                   by = c("stratum", "block"), sort = FALSE)
        print(format(d[, keep], digits = 4), row.names = FALSE)
    }
    cat("\n  structural_checks: not_evaluated. Empirical rank and conditioning\n")
    cat("  are not checked here, so a MODEL_TRANSFER_ELIGIBLE pair is identified\n")
    cat("  conditional on the declared rank-L model and its required structural\n")
    cat("  conditions, which are not empirically checked by this audit.\n")
    cat("  The design lines above are declared lines with margin. Every value\n")
    cat("  that raised a code is reported beside it; nothing here is a verdict,\n")
    cat("  and no design makes this audit refuse to return one.\n")
    invisible(x)
}

#' @export
summary.mgcca_audit <- function(object, ...) {
    print(object)
    cat("\nAvailability by block\n")
    print(format(object$availability$by_block, digits = 4), row.names = FALSE)
    if (object$strata_declared) {
        cat("\nAvailability by declared stratum\n")
        print(format(object$availability$by_stratum, digits = 4),
              row.names = FALSE)
    }
    cat("\nPairs\n")
    keep <- c("block1", "block2", "code", "n_jk", "pi_jkc_min", "s_e",
              "s_e_ess", "kappa", "path", "bottleneck")
    if (!object$weights_declared) keep <- setdiff(keep, "s_e_ess")
    print(format(object$provenance$pairs[, keep], digits = 4), row.names = FALSE)
    cat("\n  s_e is the raw joint count n_jk. ")
    if (object$weights_declared)
        cat("s_e_ess is the additional ESS-based\n  support; it does not replace s_e.\n")
    else cat("No weights were declared, so no\n  ESS-based support is reported.\n")
    cat("  bottleneck is the widest-route value: the largest support a single\n")
    cat("  route can guarantee on its weakest edge. kappa counts edge-disjoint\n")
    cat("  routes -- design-support fragility, not identification robustness.\n")
    cat("  n_jk >= L is descriptive support; it does not prove numerical rank L.\n")
    cat("\nRepresentation\n")
    rk <- c("stratum", "block", "n_stratum", "n_obs", "proportion")
    if (object$weights_declared)
        rk <- c(rk, "ess", "max_weight", "weight_concentration")
    rk <- c(rk, "advisory")
    print(format(object$representation[, rk], digits = 4), row.names = FALSE)
    cat("\n  A finite entry in this table evidences positivity for that cell.\n")
    cat("  It does not prove it: an empty cell may be empty by chance or\n")
    cat("  impossible by design, and the table cannot tell the two apart.\n")
    cat("\n  ", object$provenance$rank_assumption, "\n", sep = "")
    invisible(object)
}


## ---- the visual layer -------------------------------------------------------
##
## Frozen visual rules (V1-V7), which bind every legend, label and caption here:
##
##   V1  No traffic-light semantics. These are design properties, not verdicts,
##       so no red/green anywhere -- not for the classification, not for the
##       fragility levels, not for the declared lines.
##   V2  DIRECT_SUPPORT must be the visually strongest class: everything else
##       leans on something extra.
##   V3  Edge width encodes co-observation count and nothing else, and says so
##       in its own legend and in a frozen caption.
##   V4  Topology is not conditioning: "design-support redundancy", never
##       "robust identification".
##   V5  Support lines are DECLARED lines, drawn in a neutral style.
##   V6  STRUCTURAL_ZERO and UNSUPPORTED must be distinguishable at a glance.
##   V7  Any provenance view carries the declared rank-L condition.
##
## Colour encodes exactly one semantic variable per panel. In the matrix that is
## the frozen classification -- never support, which is why fill saturation is
## not used for counts: a thinly supported DIRECT cell must not come out looking
## like a MODEL_TRANSFER one.

## Class fills, from the house Okabe-Ito palette. Blue is the strongest and goes
## to DIRECT_SUPPORT (V2); amber to MODEL_TRANSFER_ELIGIBLE; the two unsupported
## classes are neutral greys, separated by lightness and, for the structural
## zero, by a redundant corner glyph (V6). No red, no green (V1), and blue/amber
## is the canonical colour-blind-safe pair, with the greys separated only by
## lightness, which every form of colour vision deficiency preserves.
.mgcca_av_class_fill <- c(DIRECT_SUPPORT          = "#0072B2",
                          MODEL_TRANSFER_ELIGIBLE = "#E69F00",
                          UNSUPPORTED             = "#EAEAEA",
                          STRUCTURAL_ZERO         = "#999999")
.mgcca_av_class_text <- c(DIRECT_SUPPORT          = "white",
                          MODEL_TRANSFER_ELIGIBLE = "grey10",
                          UNSUPPORTED             = "grey35",
                          STRUCTURAL_ZERO         = "grey15")

## Advisory glyph shapes. Point shapes, not Unicode, so they render in every
## device and appear in a legend of their own.
.mgcca_av_flag_shape <- c(LOW_SUPPORT = 17L, SINGLE_BRIDGE = 15L)

## Redundancy levels for the graph. Purple and sky blue -- neither red nor green
## (V1) -- and reinforced by line type, so the distinction survives in greyscale.
.mgcca_av_redundancy <- c("single bridge (kappa = 1)" = "#CC79A7",
                          "redundant (kappa >= 2)"    = "#56B4E9")

.mgcca_av_edge_caption <- "Edge width denotes raw co-observation count only."

## Categorical palettes for the panels that carry a grouping of their own. Every
## one of these encodes a VARIABLE, never decoration: node colour is the
## connected component, bar colour is the declared stratum. All entries are
## Okabe-Ito, and the two traffic-light hues are deliberately absent from the
## pool (V1), so nothing here can be read as a grade.
##
## The graph's own edges already use purple and sky blue for redundancy, so the
## node pool avoids those two: within one panel, one colour must not carry two
## meanings.
.mgcca_av_component_pal <- c("#0072B2", "#E69F00", "#F0E442", "#000000",
                             "#999999")
.mgcca_av_stratum_pal <- c("#0072B2", "#E69F00", "#CC79A7", "#56B4E9",
                           "#F0E442", "#999999")

## Sequential ramp for realised availability. Deliberately NOT the class blue:
## the availability heatmap and the classification matrix are read side by side
## in the overview, and a shared hue would invite the reader to connect a pale
## tile with a DIRECT_SUPPORT cell, which means nothing. Teal is far enough from
## both the class blue and the class amber to be unmistakable, and far enough
## from 120 degrees not to read as a green "good" (V1).
.mgcca_av_ramp_high <- "#0E7C84"

## Accent for the pattern-frequency bars: their own panel, their own colour, so
## the set does not read as one long blue document.
.mgcca_av_accent <- "#CC79A7"

## Recycle a categorical palette to the number of levels needed. With more
## groups than colours the labels carry the identity and the colour repeats --
## honest, and the matrix is the panel that carries the analysis anyway.
.mgcca_av_pal <- function(pal, levels)
    stats::setNames(rep_len(pal, length(levels)), levels)

## Captions carry the frozen wording, so they are left-aligned and allowed the
## width of the plot rather than being right-aligned into the margin.
.mgcca_av_caption_theme <- function()
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0,
                                                        colour = "grey35"))

.mgcca_av_no_grid <- function()
    ## element_blank() on the parent does not win against the children the house
    ## theme sets explicitly, so both children are blanked here.
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

## Every size here is a MULTIPLE of the panel's own base size, never a fixed
## point value: the composite has to grow and shrink with `base_size` as one
## piece, and a hard-coded 5.5pt legend beside a scaled panel would drift out of
## proportion the moment the base size moves. The multipliers are the ratios the
## composite was tuned at (its panels are drawn at base_size 8, so 5.5pt became
## 0.6875, 6.5pt became 0.8125, and so on) -- the proportions are unchanged, only
## their anchor is.
.mgcca_av_compact <- function(p, compact, base_size = 12)
    if (!compact) p else
        p + ggplot2::theme(legend.key.size = grid::unit(0.55, "lines"),
                           legend.spacing.y = grid::unit(1, "pt"),
                           legend.margin = ggplot2::margin(1, 1, 1, 1),
                           legend.text = ggplot2::element_text(
                               size = base_size * 0.6875),
                           legend.title = ggplot2::element_text(
                               size = base_size * 0.8125),
                           plot.title = ggplot2::element_text(
                               size = base_size * 1.125),
                           plot.subtitle = ggplot2::element_text(
                               size = base_size * 0.8125),
                           plot.caption = ggplot2::element_text(
                               size = base_size * 0.75))


## ---- matrix -----------------------------------------------------------------
.mgcca_av_plot_matrix <- function(x, base_size = 12, compact = FALSE) {
    b <- x$blocks
    J <- x$J
    p <- x$provenance$pairs
    st <- x$status$pairs

    ## lower triangle: column = block1, row = block2, rows running downwards
    cell <- data.frame(
        xi = match(p$block1, b),
        yi = J + 1L - match(p$block2, b),
        code = factor(p$code, levels = .mgcca_av_pair_codes),
        stringsAsFactors = FALSE)
    ## Cell text is the support quantity, and which quantity depends on the
    ## class: a direct pair is carried by its own participants, a transfer pair
    ## by the weakest edge of its route, and an unsupported pair by nothing worth
    ## printing.
    cell$label <- ifelse(
        p$code == "DIRECT_SUPPORT", as.character(p$n_jk),
        ifelse(p$code == "MODEL_TRANSFER_ELIGIBLE",
               paste0("B=", format(p$bottleneck, trim = TRUE, digits = 4)),
               ifelse(p$code == "STRUCTURAL_ZERO", "0", "")))
    cell$txt <- .mgcca_av_class_text[as.character(cell$code)]

    diag_df <- data.frame(xi = seq_len(J), yi = J + 1L - seq_len(J),
                          label = as.character(x$availability$by_block$n_j),
                          stringsAsFactors = FALSE)

    flags <- do.call(rbind, lapply(
        names(.mgcca_av_flag_shape), function(fl) {
            hit <- if (identical(fl, "LOW_SUPPORT")) st$low_support
                   else st$single_bridge
            if (!any(hit)) return(NULL)
            data.frame(xi = cell$xi[hit] +
                           if (identical(fl, "LOW_SUPPORT")) -0.28 else 0.28,
                       yi = cell$yi[hit] + 0.28,
                       flag = factor(fl, levels = names(.mgcca_av_flag_shape)),
                       stringsAsFactors = FALSE)
        }))
    zero <- cell[cell$code == "STRUCTURAL_ZERO", , drop = FALSE]

    ## Every one of the four classes needs a visible key, including the ones no
    ## pair in THIS design happens to be in -- a reader has to be able to see
    ## what the absent categories would have looked like.
    ##
    ## `drop = FALSE` alone is not enough: it keeps the level's slot and its
    ## label, but the key GLYPH is drawn from a layer's data, so a level with no
    ## rows behind it gets a label and no swatch. This layer gives every level a
    ## row. It is zero-sized, so it draws nothing in the panel, and it sits on an
    ## existing cell so it cannot move the limits either.
    seed <- data.frame(
        xi = 1L, yi = 1L,
        code = factor(.mgcca_av_pair_codes, levels = .mgcca_av_pair_codes),
        stringsAsFactors = FALSE)

    p1 <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = cell,
                           ggplot2::aes(x = .data$xi, y = .data$yi,
                                        fill = .data$code),
                           colour = "white", linewidth = 0.7) +
        ggplot2::geom_tile(data = seed,
                           ggplot2::aes(x = .data$xi, y = .data$yi,
                                        fill = .data$code),
                           width = 0, height = 0) +
        ggplot2::geom_tile(data = diag_df,
                           ggplot2::aes(x = .data$xi, y = .data$yi),
                           fill = "white", colour = "grey75", linewidth = 0.5) +
        ggplot2::geom_text(data = diag_df,
                           ggplot2::aes(x = .data$xi, y = .data$yi,
                                        label = .data$label),
                           size = base_size / 4, colour = "grey30",
                           fontface = "italic")
    if (nrow(zero))
        ## redundant encoding of the class, in its own corner so it cannot be
        ## confused with an advisory glyph (V6)
        p1 <- p1 + ggplot2::geom_point(
            data = zero,
            ggplot2::aes(x = .data$xi - 0.28, y = .data$yi - 0.28),
            shape = 4L, size = base_size / 6, stroke = 0.8, colour = "grey20",
            show.legend = FALSE)
    p1 <- p1 +
        ggplot2::geom_text(data = cell,
                           ggplot2::aes(x = .data$xi, y = .data$yi,
                                        label = .data$label),
                           colour = cell$txt, size = base_size / 4)
    if (!is.null(flags)) {
        ## The same trap as the classification legend, one aesthetic over: a
        ## flag no pair in THIS design carries would get a label and no glyph,
        ## because the key is drawn from a layer's data. This layer gives both
        ## flags a row. A zero SIZE would hide the key as well as the point, so
        ## it is drawn fully transparent instead and the guide puts the opacity
        ## back for the key alone. (No flags at all still means no legend: an
        ## empty advisory legend would be noise, not information.)
        flag_seed <- data.frame(
            xi = 1L, yi = 1L,
            flag = factor(names(.mgcca_av_flag_shape),
                          levels = names(.mgcca_av_flag_shape)),
            stringsAsFactors = FALSE)
        p1 <- p1 + ggplot2::geom_point(
            data = flags,
            ggplot2::aes(x = .data$xi, y = .data$yi, shape = .data$flag),
            size = base_size / 7, colour = "grey15") +
            ggplot2::geom_point(
                data = flag_seed,
                ggplot2::aes(x = .data$xi, y = .data$yi, shape = .data$flag),
                size = base_size / 7, colour = "grey15", alpha = 0) +
            ggplot2::scale_shape_manual(values = .mgcca_av_flag_shape,
                                        name = "Advisory flags", drop = FALSE)
    }
    p1 +
        ggplot2::scale_fill_manual(values = .mgcca_av_class_fill,
                                   name = "Pair classification", drop = FALSE) +
        ## The two neutral classes are nearly white, so their legend keys need a
        ## border or they vanish against the legend background; and the seeded
        ## advisory glyphs need their opacity back for the key. Both overrides
        ## live in one guides() call so neither can quietly drop the other.
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                override.aes = list(colour = "grey60")),
            shape = ggplot2::guide_legend(
                override.aes = list(alpha = 1))) +
        ggplot2::scale_x_continuous(breaks = seq_len(J), labels = b,
                                    limits = c(0.5, J + 0.5), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = seq_len(J), labels = rev(b),
                                    limits = c(0.5, J + 0.5), expand = c(0, 0)) +
        ggplot2::coord_equal() +
        ggplot2::labs(
            x = NULL, y = NULL, title = "Pair classification",
            subtitle = paste("lower triangle; the diagonal carries each block's",
                             "marginal n_j"),
            caption = paste("Colour encodes the frozen classification only, not",
                            "how much support a pair has.\nCell text: n_jk",
                            "(direct), B = route bottleneck (model transfer),",
                            "0 (structural zero).")) +
        .mgcca_theme(base_size) +
        .mgcca_av_no_grid() +
        .mgcca_av_caption_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                           hjust = 1))
}


## ---- graph ------------------------------------------------------------------
## Deterministic circular layout: connected component, then degree (descending,
## so a hub anchors its component), then block name. No force-directed layout --
## a layout that moves between runs cannot be read as evidence of anything.
.mgcca_av_node_order <- function(x) {
    deg <- as.integer(rowSums(x$graph$adjacency))
    order(x$graph$component, -deg, x$blocks, method = "radix")
}

.mgcca_av_plot_graph <- function(x, pair = NULL, base_size = 12,
                                 compact = FALSE) {
    J <- x$J
    ord <- .mgcca_av_node_order(x)
    blocks <- x$blocks[ord]
    ang <- seq(0, 2 * pi, length.out = J + 1L)[seq_len(J)] - pi / 2
    ## Node colour is the connected component -- the audit's central object, and
    ## the thing the reader most needs to see at a glance: blocks of one colour
    ## can reach each other, blocks of different colours cannot.
    comp_lev <- sort(unique(x$graph$component))
    nodes <- data.frame(block = blocks, x = cos(ang), y = sin(ang),
                        component = factor(x$graph$component[blocks],
                                           levels = comp_lev),
                        stringsAsFactors = FALSE)
    ## Labels sit OUTSIDE their node: a label centred on a node on the left or
    ## right of the circle lands on top of it. Anchor by which side it is on.
    nodes$hj <- ifelse(nodes$x > 0.15, 0, ifelse(nodes$x < -0.15, 1, 0.5))
    nodes$lx <- nodes$x * 1.10
    nodes$ly <- nodes$y * 1.10 +
        ifelse(abs(nodes$x) < 0.15, sign(nodes$y) * 0.12, 0)

    route <- NULL
    if (!is.null(pair)) route <- .mgcca_av_route(x, pair)

    e <- x$graph$edges
    seg <- NULL
    if (nrow(e)) {
        i1 <- match(e$block1, nodes$block)
        i2 <- match(e$block2, nodes$block)
        seg <- data.frame(
            x = nodes$x[i1], y = nodes$y[i1],
            xend = nodes$x[i2], yend = nodes$y[i2],
            s_e = e$s_e,
            redundancy = ifelse(e$is_bridge, names(.mgcca_av_redundancy)[1L],
                                names(.mgcca_av_redundancy)[2L]),
            stringsAsFactors = FALSE)
        seg$on_route <- if (is.null(route)) TRUE else
            vapply(seq_len(nrow(e)), function(i)
                any((route$edges$from == e$block1[i] &
                         route$edges$to == e$block2[i]) |
                        (route$edges$from == e$block2[i] &
                             route$edges$to == e$block1[i])), logical(1))
    }

    p <- ggplot2::ggplot()
    if (!is.null(seg)) {
        if (!is.null(route)) {
            off <- seg[!seg$on_route, , drop = FALSE]
            if (nrow(off))
                p <- p + ggplot2::geom_segment(
                    data = off,
                    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                 yend = .data$yend),
                    colour = "grey85", linewidth = 0.5)
            seg <- seg[seg$on_route, , drop = FALSE]
        }
        p <- p + ggplot2::geom_segment(
            data = seg,
            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                         yend = .data$yend, linewidth = .data$s_e,
                         colour = .data$redundancy,
                         linetype = .data$redundancy)) +
            ## a handful of breaks, not one per distinct count: a tall key
            ## column crowds the panel it is meant to explain
            ggplot2::scale_linewidth(range = c(0.4, 2.4),
                                     name = "jointly observed n",
                                     breaks = function(lims) {
                                         k <- if (compact) 3L else 4L
                                         b <- pretty(lims, n = k)
                                         ## a narrow integer range makes
                                         ## pretty() return one break per unit
                                         if (length(b) > k + 1L)
                                             b <- b[round(seq(1, length(b),
                                                              length.out = k))]
                                         b
                                     }) +
            ggplot2::scale_colour_manual(values = .mgcca_av_redundancy,
                                         name = "Design-support redundancy",
                                         drop = FALSE) +
            ggplot2::scale_linetype_manual(
                values = stats::setNames(c("22", "solid"),
                                         names(.mgcca_av_redundancy)),
                name = "Design-support redundancy", drop = FALSE)
    }

    node_hi <- rep(TRUE, nrow(nodes))
    if (!is.null(route)) node_hi <- nodes$block %in% route$path
    ## Off-route nodes fade rather than lose their colour: the component they
    ## belong to is still true of them, so the honest de-emphasis is opacity,
    ## not a recolour to grey.
    p <- p +
        ggplot2::geom_point(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$y, fill = .data$component),
            shape = 21L, colour = "grey25", stroke = 0.6,
            size = ifelse(node_hi, base_size / 3, base_size / 5),
            alpha = ifelse(node_hi, 1, 0.3)) +
        ggplot2::scale_fill_manual(
            values = .mgcca_av_pal(.mgcca_av_component_pal, comp_lev),
            name = "Connected component") +
        ggplot2::guides(fill = ggplot2::guide_legend(
            override.aes = list(alpha = 1, size = base_size / 3.5))) +
        ggplot2::geom_text(
            data = nodes,
            ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$block),
            hjust = nodes$hj, size = base_size / 3.6,
            colour = ifelse(node_hi, "grey10", "grey60"))

    ## In the overview the graph is a small panel beside a legend column, and
    ## the full subtitle runs straight into the legend titles. The short form is
    ## used there; the standalone type keeps the whole sentence, and the frozen
    ## caption (V3) is never shortened in either.
    sub <- if (compact) "an edge is a pair supported directly"
           else paste("an edge is a pair supported directly;",
                      "node order is component, then degree, then name")
    cap <- .mgcca_av_edge_caption
    if (!is.null(route)) {
        wk <- which.min(route$edges$s_e)
        i1 <- match(route$edges$from[wk], nodes$block)
        i2 <- match(route$edges$to[wk], nodes$block)
        bn <- data.frame(x = (nodes$x[i1] + nodes$x[i2]) / 2,
                         y = (nodes$y[i1] + nodes$y[i2]) / 2,
                         label = paste0("bottleneck  B = ",
                                        format(route$bottleneck, trim = TRUE,
                                               digits = 4)),
                         stringsAsFactors = FALSE)
        p <- p + ggplot2::geom_label(
            data = bn,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
            size = base_size / 4, linewidth = 0.2, colour = "grey10",
            fill = "white", alpha = 0.85)
        ## V7: the transfer condition travels with the picture, always.
        sub <- "Model-transfer-eligible under declared rank-L assumption"
        if (isTRUE(route$alternative_paths))
            sub <- paste0(sub, "; alternative routes exist")
        cap <- paste0(cap, "\nRoute shown: ",
                      paste(route$path, collapse = " -> "),
                      ". L = ", x$L, " as declared.")
    }

    .mgcca_av_compact(
        p +
            ggplot2::coord_equal(xlim = c(-1.95, 1.95), ylim = c(-1.4, 1.4)) +
            ggplot2::labs(x = NULL, y = NULL, title = "Support graph",
                          subtitle = sub, caption = cap) +
            .mgcca_theme(base_size) +
            .mgcca_av_no_grid() +
            .mgcca_av_caption_theme() +
            ggplot2::theme(axis.text = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank()),
        compact, base_size)
}

## The stored route for a pair. Nothing is recomputed here: the path, its edges
## and its bottleneck are read back exactly as the audit recorded them.
.mgcca_av_route <- function(x, pair) {
    if (!is.character(pair) || length(pair) != 2L)
        stop("`pair` must be two block names, e.g. ",
             "pair = c(\"", x$blocks[1L], "\", \"", x$blocks[x$J], "\").",
             call. = FALSE)
    unknown <- setdiff(pair, x$blocks)
    if (length(unknown))
        stop("unknown block name(s): ", paste(unknown, collapse = ", "),
             ". Available: ", paste(x$blocks, collapse = ", "), call. = FALSE)
    if (pair[1L] == pair[2L])
        stop("`pair` must name two different blocks.", call. = FALSE)
    pp <- x$provenance$pairs
    i <- which((pp$block1 == pair[1L] & pp$block2 == pair[2L]) |
                   (pp$block1 == pair[2L] & pp$block2 == pair[1L]))
    code <- pp$code[i]
    if (!identical(code, "MODEL_TRANSFER_ELIGIBLE")) {
        ## the class-specific half of the message is assembled into a local
        ## first: pasting inside the condition signal itself is what the
        ## Bioconductor coding-practice check asks packages not to do.
        why <- switch(code,
                      DIRECT_SUPPORT = paste0(
                          "It is carried directly by the ", pp$n_jk[i],
                          " participants observed in both blocks, so there is ",
                          "no transfer route to draw."),
                      UNSUPPORTED = paste0(
                          "It has no route at all: ", pp$n_jk[i],
                          " jointly observed participants, and the declared ",
                          "availability line is not met in every declared ",
                          "stratum."),
                      STRUCTURAL_ZERO = paste0(
                          "It has no route at all, and no participant anywhere ",
                          "is observed in both blocks."))
        stop("the provenance overlay draws the route a MODEL_TRANSFER_ELIGIBLE ",
             "pair depends on, and ", pair[1L], "--", pair[2L], " is coded ",
             code, ". ", why,
             " Use plot(x, type = \"matrix\") for the whole classification.",
             call. = FALSE)
    }
    x$provenance$paths[[i]]
}


## ---- availability -----------------------------------------------------------
## The fill scale is fixed to [0, 1] and never auto-scaled to the observed
## range: a design whose availabilities all sit between 0.90 and 0.94 must not
## be painted to look like the full spread.
.mgcca_av_plot_availability <- function(x, base_size = 12, compact = FALSE) {
    d <- x$availability$by_stratum
    ## the block order is the object's, not the alphabet's, so this panel and
    ## the classification matrix name their columns in the same order
    d$block <- factor(d$block, levels = x$blocks)
    d$stratum <- factor(d$stratum, levels = rev(x$strata_levels))
    .mgcca_av_compact(
        ggplot2::ggplot(d, ggplot2::aes(x = .data$block, y = .data$stratum,
                                        fill = .data$pi_hat_jc)) +
            ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
            ggplot2::geom_text(ggplot2::aes(label = .data$n_jc),
                               size = base_size / 4) +
            ggplot2::scale_fill_gradient(
                low = "white", high = .mgcca_av_ramp_high,
                limits = c(0, 1), name = expression(hat(pi)[jc])) +
            ggplot2::labs(
                x = NULL, y = NULL, title = "Realised availability",
                subtitle = "cell labels are n_jc; the scale is fixed to [0, 1]",
                caption = sprintf(
                    paste("The declared availability line is %s -- a declared",
                          "line with margin, not a\nclass boundary and not a",
                          "threshold at which anything fails."),
                    format(x$design_lines$availability_min))) +
            .mgcca_theme(base_size) +
            .mgcca_av_no_grid() +
            .mgcca_av_caption_theme() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                               hjust = 1)),
        compact, base_size)
}


## ---- representation ---------------------------------------------------------
## Two aligned facets, never a dual axis: an effective count and a proportion do
## not share a scale, and drawing them on one would invite exactly the
## comparison that is not available.
.mgcca_av_plot_representation <- function(x, base_size = 12, compact = FALSE) {
    r <- x$representation
    lab <- if (x$strata_declared) paste(r$stratum, r$block, sep = " / ")
           else r$block
    ordr <- order(match(r$stratum, x$strata_levels), match(r$block, x$blocks),
                  method = "radix")
    lev <- rev(lab[ordr])

    count_name <- if (x$weights_declared) "Effective sample size (Kish)"
                  else "Observed count (no weights declared)"
    count_val <- if (x$weights_declared) r$ess else as.numeric(r$n_obs)
    ## the stratum size travels with every effective count, because comparing
    ## raw ESS across strata of different sizes is the mistake this annotation
    ## exists to head off
    count_lab <- if (x$weights_declared)
        sprintf("ESS %.1f (n=%d)", r$ess, r$n_stratum)
        else sprintf("%d (n=%d)", r$n_obs, r$n_stratum)

    ## Bar colour is the declared stratum, the same colour in both panels, so a
    ## stratum can be followed from its effective count to its retention without
    ## re-reading the labels. One variable, one scale.
    strat <- factor(r$stratum, levels = x$strata_levels)
    d <- rbind(
        data.frame(cell = factor(lab, levels = lev), measure = count_name,
                   value = count_val, label = count_lab, stratum = strat,
                   stringsAsFactors = FALSE),
        data.frame(cell = factor(lab, levels = lev),
                   measure = "Retention (proportion of stratum)",
                   value = r$proportion,
                   label = sprintf("%.2f", r$proportion), stratum = strat,
                   stringsAsFactors = FALSE))
    d$measure <- factor(d$measure,
                        levels = c(count_name,
                                   "Retention (proportion of stratum)"))
    lines <- data.frame(
        measure = factor(c(count_name, "Retention (proportion of stratum)"),
                         levels = levels(d$measure)),
        at = c(x$design_lines$ess_min, x$design_lines$availability_min),
        stringsAsFactors = FALSE)
    lines$label <- "declared design line"
    ## Per-facet headroom: the count panel needs room for its annotation, the
    ## retention panel must still stop at 1 -- a proportion axis running past 1
    ## would be nonsense. facet expansion is shared, an invisible point is not.
    room <- rbind(
        data.frame(measure = factor(count_name, levels = levels(d$measure)),
                   value = max(c(count_val, x$design_lines$ess_min),
                               na.rm = TRUE) * 1.68,
                   cell = factor(lev[1L], levels = lev),
                   stringsAsFactors = FALSE),
        data.frame(measure = factor("Retention (proportion of stratum)",
                                    levels = levels(d$measure)),
                   value = 1.15, cell = factor(lev[1L], levels = lev),
                   stringsAsFactors = FALSE))

    .mgcca_av_compact(
        ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$cell)) +
            ggplot2::geom_col(ggplot2::aes(fill = .data$stratum),
                              width = 0.65, alpha = 0.9) +
            ggplot2::scale_fill_manual(
                values = .mgcca_av_pal(.mgcca_av_stratum_pal,
                                       x$strata_levels),
                name = "Declared stratum") +
            ## with nothing declared there is one stratum and the legend would
            ## be a colour key for the word "everybody"
            ggplot2::guides(fill = if (x$strata_declared) "legend" else "none") +
            ## V5: a declared line, drawn neutrally -- dashed grey, never a red
            ## failure bar.
            ggplot2::geom_vline(data = lines,
                                ggplot2::aes(xintercept = .data$at),
                                linetype = "dashed", colour = "grey45",
                                linewidth = 0.5) +
            ggplot2::geom_text(data = if (compact) lines[0L, ] else lines,
                               ggplot2::aes(x = .data$at, y = length(lev) + 0.4,
                                            label = .data$label),
                               inherit.aes = FALSE, hjust = -0.05,
                               vjust = 0, size = base_size / 5,
                               colour = "grey45") +
            ggplot2::geom_text(ggplot2::aes(label = .data$label),
                               hjust = -0.06, size = base_size / 4.5,
                               colour = "grey20") +
            ggplot2::geom_blank(data = room,
                                ggplot2::aes(x = .data$value, y = .data$cell)) +
            ggplot2::facet_wrap(~ .data$measure, nrow = 1L, scales = "free_x") +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(
                mult = c(0, 0.02))) +
            ggplot2::labs(x = NULL, y = NULL, title = "Representation",
                          subtitle = paste("same ordering in both panels;",
                                           "the two scales are not comparable"),
                          caption = paste("n is the declared stratum's size.",
                                          "Dashed lines are declared design",
                                          "lines with margin,\nnot thresholds",
                                          "at which anything fails.")) +
            .mgcca_theme(base_size) +
            .mgcca_av_caption_theme() +
            ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()),
        compact, base_size)
}


## ---- patterns ---------------------------------------------------------------
## An UpSet-like view of which combinations of blocks the participants actually
## carry, aggregated at plot time from the availability table the object already
## holds. Nothing is added to the object for it.
##
## Both panels share one continuous x: block positions are negative and counts
## are positive, so a single scale can hand each panel its own breaks and labels
## without a second scale package.
.mgcca_av_plot_patterns <- function(x, base_size = 12, compact = FALSE) {
    A <- x$table
    b <- x$blocks
    J <- x$J
    key <- apply(A, 1L, function(r) paste(as.integer(r), collapse = ""))
    tb <- table(key)
    pat <- names(tb)
    cnt <- as.integer(tb)
    ## deterministic: frequency descending, then the pattern key lexicographic
    ord <- order(-cnt, pat, method = "radix")
    pat <- pat[ord]
    cnt <- cnt[ord]
    lev <- rev(pat)

    memb <- do.call(rbind, lapply(seq_along(pat), function(i) {
        present <- strsplit(pat[i], "", fixed = TRUE)[[1L]] == "1"
        data.frame(pattern = factor(pat[i], levels = lev),
                   xi = seq_len(J) - J - 1L, present = present,
                   panel = "Blocks present", stringsAsFactors = FALSE)
    }))
    link <- do.call(rbind, lapply(seq_along(pat), function(i) {
        present <- which(strsplit(pat[i], "", fixed = TRUE)[[1L]] == "1")
        if (length(present) < 2L) return(NULL)
        data.frame(pattern = factor(pat[i], levels = lev),
                   x = min(present) - J - 1L, xend = max(present) - J - 1L,
                   panel = "Blocks present", stringsAsFactors = FALSE)
    }))
    bars <- data.frame(pattern = factor(pat, levels = lev), value = cnt,
                       panel = "Individuals", stringsAsFactors = FALSE)
    pan <- c("Blocks present", "Individuals")
    memb$panel <- factor(memb$panel, levels = pan)
    bars$panel <- factor(bars$panel, levels = pan)
    if (!is.null(link)) link$panel <- factor(link$panel, levels = pan)

    ## The two facets share one scale, so this has to tell them apart from the
    ## limits alone -- and the sign of the upper limit will not do it. The block
    ## panel's data stop at -1, but the 20% right-hand expansion the count
    ## labels need is a fraction of the SPAN, so it grows with the number of
    ## blocks: at J = 6 it reaches exactly 0 and from there on the panel tests
    ## as positive, falls through to the count branch, and the block names
    ## silently disappear. What does separate the two at any J is which side of
    ## zero the panel sits on: the blocks are almost entirely to the left of it,
    ## the counts almost entirely to the right.
    brk <- function(lims) {
        if (!length(lims) || is.na(lims[2L])) return(numeric(0))
        if (lims[2L] < (lims[2L] - lims[1L]) / 2) seq_len(J) - J - 1L
        else pretty(c(0, lims[2L]))
    }
    lab <- function(v) {
        out <- ifelse(is.na(v), "",
                      ifelse(v < 0, b[pmin(pmax(v + J + 1L, 1L), J)],
                             format(v, trim = TRUE)))
        as.character(out)
    }

    p <- ggplot2::ggplot()
    if (!is.null(link))
        p <- p + ggplot2::geom_segment(
            data = link,
            ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$pattern,
                         yend = .data$pattern),
            colour = "grey60", linewidth = 0.5)
    .mgcca_av_compact(
        p +
            ggplot2::geom_point(
                data = memb,
                ggplot2::aes(x = .data$xi, y = .data$pattern,
                             colour = .data$present),
                size = base_size / 4.5, show.legend = FALSE) +
            ggplot2::geom_col(
                data = bars,
                ggplot2::aes(x = .data$value, y = .data$pattern),
                fill = .mgcca_av_accent, width = 0.6, alpha = 0.9) +
            ggplot2::geom_text(
                data = bars,
                ggplot2::aes(x = .data$value, y = .data$pattern,
                             label = .data$value),
                hjust = -0.25, size = base_size / 4.5, colour = "grey20") +
            ggplot2::scale_colour_manual(
                values = c(`TRUE` = "grey15", `FALSE` = "grey85")) +
            ## Pin the row order explicitly. The segment layer comes first but
            ## carries no row for a single-block pattern (nothing to connect),
            ## so leaving the discrete scale to collect levels layer by layer
            ## appends exactly those patterns at the end -- i.e. at the top of
            ## the plot, out of frequency order. Seen on a real 5-block design.
            ggplot2::scale_y_discrete(limits = lev) +
            ggplot2::scale_x_continuous(breaks = brk, labels = lab,
                                        expand = ggplot2::expansion(
                                            mult = c(0.14, 0.2))) +
            ggplot2::facet_grid(cols = ggplot2::vars(.data$panel),
                                scales = "free_x") +
            ggplot2::labs(
                x = NULL, y = NULL, title = "Availability patterns",
                subtitle = paste("which combinations of blocks the",
                                 "participants carry, most frequent first"),
                caption = paste("Counts are participants. A pattern that no",
                                "participant carries is not shown: an absent\n",
                                "combination may be absent by chance or",
                                "impossible by design, and this cannot tell",
                                "them apart.")) +
            .mgcca_theme(base_size) +
            .mgcca_av_caption_theme() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.ticks.y = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_blank(),
                           ## Angled so long block names stay apart at J >= 6;
                           ## the numeric count panel tolerates the same angle.
                           axis.text.x = ggplot2::element_text(size =
                                                                   base_size * 0.7,
                                                               angle = 30,
                                                               hjust = 1)),
        compact, base_size)
}


## ---- overview ---------------------------------------------------------------
## Composition with base grid viewports. `grid` ships with R, so the composite
## costs no dependency; ggplot2's own print method takes a viewport, which is
## all the layout needs.
##
## Curated and stable, not "everything": the classification matrix anchors it
## because it is the precise view, and the graph is the small orientation panel
## rather than the front door. The provenance overlay is deliberately absent --
## it needs a pair chosen by the reader.
## Four panels share one device, so the composite draws them at a FRACTION of
## the caller's base size rather than at a size of its own: kept as a ratio, the
## composite grows and shrinks with `base_size` instead of ignoring it, and the
## four panels stay in proportion with each other at any size.
.mgcca_av_overview_scale <- 2 / 3

.mgcca_av_overview <- function(x, base_size = 12 * .mgcca_av_overview_scale) {
    panels <- list(
        matrix = .mgcca_av_plot_matrix(x, base_size, compact = TRUE),
        graph = .mgcca_av_plot_graph(x, NULL, base_size, compact = TRUE),
        availability = .mgcca_av_plot_availability(x, base_size, compact = TRUE),
        representation = .mgcca_av_plot_representation(x, base_size,
                                                       compact = TRUE))
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
        4L, 2L, widths = grid::unit(c(1.15, 1), "null"),
        heights = grid::unit(rep(1, 4L), "null"))))
    vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
    print(panels$matrix, vp = vp(seq_len(3L), 1L))
    print(panels$graph, vp = vp(4L, 1L))
    print(panels$availability, vp = vp(seq_len(2L), 2L))
    print(panels$representation, vp = vp(3:4, 2L))
    grid::popViewport()
    invisible(panels)
}


#' Plot an availability audit
#'
#' @description Six views of one audit. \code{"overview"} (the default, and what
#'   a bare \code{plot()} gives you) is a curated four-panel composite: the
#'   classification matrix as the anchor, realised availability, representation,
#'   and a small support graph. The other five are single panels and come back
#'   as \code{ggplot} objects you can style, save or assemble yourself.
#'
#' @details
#' \describe{
#'   \item{\code{"overview"}}{Four panels drawn straight to the device through
#'     base \code{grid} viewports, returning the panel list invisibly. Curated
#'     and deliberately stable: it will not silently grow a panel when a new
#'     type is added. Give it room -- roughly 11 by 9 inches.}
#'   \item{\code{"matrix"}}{One triangle of the pair classification, the
#'     diagonal carrying each block's marginal \eqn{n_j}. Fill encodes the
#'     frozen classification and nothing else; the cell text carries the support
#'     quantity that belongs to that class (\eqn{n_{jk}} for a direct pair,
#'     \code{B =} the route bottleneck for a transfer pair, \code{0} for a
#'     structural zero, blank for unsupported); advisory flags are corner
#'     glyphs with their own legend.}
#'   \item{\code{"graph"}}{The support graph on a circle, node order fixed by
#'     connected component, then degree, then block name. With \code{pair} it
#'     becomes a provenance overlay: the stored route between those two blocks
#'     is emphasised, everything else recedes, and the bottleneck is marked on
#'     its weakest edge.}
#'   \item{\code{"availability"}}{Realised availability per block and declared
#'     stratum, on a fill scale fixed to \eqn{[0, 1]}.}
#'   \item{\code{"representation"}}{Two aligned panels in the same order: an
#'     effective count (Kish ESS where weights are declared, the observed count
#'     otherwise, each annotated with its stratum size) and the retention
#'     proportion. Never one dual axis -- a count and a proportion do not share
#'     a scale.}
#'   \item{\code{"patterns"}}{Which combinations of blocks the participants
#'     actually carry, as an UpSet-like dot matrix with frequency bars, ordered
#'     by frequency and then lexicographically. Aggregated at plot time from the
#'     availability table the object already holds. Not part of
#'     \code{"overview"}: it is as tall as the design is varied.}
#' }
#'
#' \strong{What the colours mean, and what they do not.} Colour encodes exactly
#' one variable per panel, and in the matrix that variable is the frozen
#' classification -- never how much support a pair has, which is why a thinly
#' supported \code{DIRECT_SUPPORT} cell is the same blue as a richly supported
#' one, with the count printed in it. There is no traffic-light scale here and
#' no red or green anywhere: these are properties of a design, not verdicts on
#' it. \code{DIRECT_SUPPORT} is the visually strongest class because every other
#' one leans on something extra.
#'
#' Every other panel colours by a variable too, never for decoration: graph
#' nodes by connected component (blocks of one colour can reach each other,
#' blocks of different colours cannot), representation bars by declared stratum
#' (the same colour in both panels, so a stratum can be followed across them),
#' availability by realised proportion on a teal ramp deliberately unlike the
#' classification blue, and the pattern bars in their own accent. Every palette
#' is colour-blind-safe, and the two traffic-light hues appear in none of them.
#'
#' \strong{Edge width, and only edge width.} In the graph, the width of an edge
#' is the raw count of participants observed in both blocks and nothing else.
#' What the graph shows is \emph{design-support redundancy} -- how many
#' independent routes the design gives you between two blocks -- which is a
#' statement about topology, not about conditioning or how well anything is
#' estimated.
#'
#' \strong{Declared lines.} The dashed lines in the representation panel are the
#' declared design lines, drawn neutrally and labelled as such. They are set by
#' the analyst, they carry margin, and they are not thresholds at which anything
#' fails.
#'
#' \strong{Provenance carries its condition.} A \code{pair} overlay is only
#' available for a \code{MODEL_TRANSFER_ELIGIBLE} pair, and it always says so:
#' such a pair is identified conditional on the declared rank-\code{L} model and
#' its required structural conditions, which are not empirically checked by this
#' audit. If the stored provenance records that alternative routes exist, the
#' subtitle says so rather than letting one route look privileged.
#'
#' @param x An object of class \code{"mgcca_audit"}.
#' @param type One of \code{"overview"} (default), \code{"matrix"},
#'   \code{"graph"}, \code{"availability"}, \code{"representation"} or
#'   \code{"patterns"}.
#' @param pair Optional two block names, for \code{type = "graph"} only: draw
#'   the stored transfer route between them. Only the stored provenance is
#'   read; no route is computed here.
#' @param base_size Base font size passed to the theme. Every label, legend and
#'   caption in every panel is a multiple of it, so one number resizes the whole
#'   view. \code{"overview"} draws its four panels at a fixed fraction of it,
#'   since they share one device.
#' @param ... Ignored.
#' @return For every type except \code{"overview"}, a \code{ggplot} object. For
#'   \code{"overview"}, the four panels as a named list, returned invisibly
#'   after the composite is drawn.
#' @seealso \code{\link{mgcca_audit}}
#' @examples
#' # a design where one block is reachable only through a single edge
#' ids <- sprintf("ind%02d", 1:85)
#' avail <- cbind(methylation = seq_along(ids) %in% 1:60,
#'                expression  = seq_along(ids) %in% 1:60,
#'                proteins    = rep(TRUE, length(ids)),
#'                metabolites = seq_along(ids) %in% 61:85)
#' rownames(avail) <- ids
#' a <- mgcca_audit(avail, L = 2,
#'                  strata = rep(c("cohort A", "cohort B"), length.out = 85),
#'                  weights = ifelse(seq_along(ids) %in% 1:60, 1, 3))
#'
#' # the curated composite (give it room: ~11 x 9 inches)
#' panels <- plot(a)
#' names(panels)
#'
#' # the precise view, on its own
#' plot(a, type = "matrix")
#'
#' # which combinations of blocks people actually carry
#' plot(a, type = "patterns")
#'
#' # the route a transfer-eligible pair depends on
#' plot(a, type = "graph", pair = c("methylation", "metabolites"))
#' @export
plot.mgcca_audit <- function(x, type = c("overview", "matrix", "graph",
                                         "availability", "representation",
                                         "patterns"),
                             pair = NULL, base_size = 13, ...) {
    type <- match.arg(type)
    if (!is.null(pair) && !identical(type, "graph"))
        stop("`pair` draws a transfer route on the support graph, so it needs ",
             "type = \"graph\"; it has no meaning for type = \"", type, "\".",
             call. = FALSE)
    switch(type,
           overview = .mgcca_av_overview(x, base_size *
                                             .mgcca_av_overview_scale),
           matrix = .mgcca_av_plot_matrix(x, base_size),
           graph = .mgcca_av_plot_graph(x, pair, base_size),
           availability = .mgcca_av_plot_availability(x, base_size),
           representation = .mgcca_av_plot_representation(x, base_size),
           patterns = .mgcca_av_plot_patterns(x, base_size))
}
