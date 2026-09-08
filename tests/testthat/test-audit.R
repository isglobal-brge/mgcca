# Pre-fit availability and design audit.
#
# Every graph case below is hand-built and hand-computed: the components, the
# edges, the edge connectivities kappa, the widest-route bottlenecks B, the Kish
# ESS values. Nothing is compared against a previous run of the same code, so a
# change of behaviour cannot be absorbed by a regenerated expectation.
#
# The audit reads presence and absence only. It never calls the estimator, so
# none of these tests fit anything, and no number a fit produces can move
# because of them.


# ---- hand-built designs -----------------------------------------------------

# CHAIN A1 -- B -- C. 40 individuals carry (A1, B), 60 carry (B, C), nobody
# carries (A1, C).
#   n_A1B = 40, n_BC = 60, n_A1C = 0
#   one component; A1--C has no edge but the route A1 -> B -> C, whose weakest
#   edge is 40, so B_A1C = 40. Every edge is a bridge (kappa = 1 throughout) and
#   B is the only articulation block.
ad_chain <- function() {
    cbind(A1 = rep(c(TRUE, FALSE), c(40, 60)),
          B  = rep(TRUE, 100),
          C  = rep(c(FALSE, TRUE), c(40, 60)))
}

# CYCLE A -- B -- C -- D -- A, with deliberately unequal supports so that every
# widest route is unique.
#   n_AB = 20, n_BC = 30, n_CD = 40, n_DA = 10, n_AC = n_BD = 0
#   kappa = 2 for every pair (a cycle), no bridge, no articulation block.
#   B_AC = min(20, 30) = 20 by A -> B -> C (A -> D -> C gives only 10)
#   B_BD = min(30, 40) = 30 by B -> C -> D (B -> A -> D gives only 10)
#   B_AD = 20 by A -> B -> C -> D, which BEATS its own direct edge of 10.
ad_cycle <- function() {
    m <- matrix(FALSE, 100, 4, dimnames = list(NULL, c("A", "B", "C", "D")))
    m[1:20, c(1, 2)] <- TRUE
    m[21:50, c(2, 3)] <- TRUE
    m[51:90, c(3, 4)] <- TRUE
    m[91:100, c(4, 1)] <- TRUE
    m
}

# DISCONNECTED: 50 individuals carry (A, B), 50 carry C alone.
#   two components; n_AC = n_BC = 0, so both are structural zeros.
ad_disconnected <- function() {
    m <- matrix(FALSE, 100, 3, dimnames = list(NULL, c("A", "B", "C")))
    m[1:50, c(1, 2)] <- TRUE
    m[51:100, 3] <- TRUE
    m
}

# SINGLE BRIDGE: a triangle A, B, C carried by 30 individuals, plus D attached
# to C alone by 20 more.
#   n_AB = n_AC = n_BC = 30, n_CD = 20, n_AD = n_BD = 0
#   kappa: 2 inside the triangle, 1 for anything involving D. C--D is the only
#   bridge, C the only articulation block. B_AD = B_BD = min(30, 20) = 20.
ad_bridge <- function() {
    m <- matrix(FALSE, 50, 4, dimnames = list(NULL, c("A", "B", "C", "D")))
    m[1:30, c(1, 2, 3)] <- TRUE
    m[31:50, c(3, 4)] <- TRUE
    m
}

# UNSUPPORTED rather than STRUCTURAL_ZERO: A and B do share 4 individuals out of
# 200, which is below the declared availability line of 0.05, so there is no
# edge and no route -- but the joint cell is not empty.
ad_thin <- function() {
    m <- matrix(FALSE, 200, 2, dimnames = list(NULL, c("A", "B")))
    m[1:4, c(1, 2)] <- TRUE
    m[5:102, 1] <- TRUE
    m[103:200, 2] <- TRUE
    m
}

pair_of <- function(a, j, k, col) {
    p <- a$provenance$pairs
    p[[col]][(p$block1 == j & p$block2 == k) | (p$block1 == k & p$block2 == j)]
}


# ---- (1) graph classification ----------------------------------------------

test_that("the chain has one component, three bridges and one articulation block", {
    a <- mgcca_audit(ad_chain(), L = 2)

    expect_equal(unname(a$graph$component), c(1L, 1L, 1L))
    expect_equal(a$graph$n_components, 1L)
    expect_true(a$graph$connected)
    expect_false(a$graph$complete)

    expect_equal(nrow(a$graph$edges), 2L)
    expect_true(all(a$graph$edges$is_bridge))
    expect_equal(a$graph$kappa,
                 matrix(c(NA, 1L, 1L, 1L, NA, 1L, 1L, 1L, NA), 3, 3,
                        dimnames = list(c("A1", "B", "C"), c("A1", "B", "C"))))
    expect_equal(unname(a$graph$articulation), c(FALSE, TRUE, FALSE))

    expect_equal(pair_of(a, "A1", "B", "code"), "DIRECT_SUPPORT")
    expect_equal(pair_of(a, "B", "C", "code"), "DIRECT_SUPPORT")
    expect_equal(pair_of(a, "A1", "C", "code"), "MODEL_TRANSFER_ELIGIBLE")
    expect_equal(pair_of(a, "A1", "B", "n_jk"), 40L)
    expect_equal(pair_of(a, "B", "C", "n_jk"), 60L)
    expect_equal(pair_of(a, "A1", "C", "n_jk"), 0L)
})

test_that("the cycle has kappa 2 everywhere, no bridge and no articulation block", {
    a <- mgcca_audit(ad_cycle(), L = 2)

    expect_equal(a$graph$n_components, 1L)
    expect_equal(nrow(a$graph$edges), 4L)
    expect_equal(nrow(a$graph$bridges), 0L)
    expect_false(any(a$graph$articulation))
    k <- a$graph$kappa
    expect_true(all(k[upper.tri(k)] == 2L))
    expect_true(all(is.na(diag(k))))

    expect_equal(pair_of(a, "A", "C", "code"), "MODEL_TRANSFER_ELIGIBLE")
    expect_equal(pair_of(a, "B", "D", "code"), "MODEL_TRANSFER_ELIGIBLE")
    expect_equal(sum(a$provenance$pairs$code == "DIRECT_SUPPORT"), 4L)
    expect_true(all(a$provenance$pairs$alternative_paths))
    expect_false(any(a$provenance$pairs$single_bridge))
})

test_that("the disconnected design splits into two components", {
    a <- mgcca_audit(ad_disconnected(), L = 2)

    expect_equal(unname(a$graph$component), c(1L, 1L, 2L))
    expect_equal(a$graph$n_components, 2L)
    expect_false(a$graph$connected)
    expect_equal(pair_of(a, "A", "B", "code"), "DIRECT_SUPPORT")
    expect_equal(pair_of(a, "A", "C", "code"), "STRUCTURAL_ZERO")
    expect_equal(pair_of(a, "B", "C", "code"), "STRUCTURAL_ZERO")
    expect_equal(pair_of(a, "A", "C", "kappa"), 0L)
    expect_equal(pair_of(a, "B", "C", "kappa"), 0L)
    expect_true(is.na(pair_of(a, "A", "C", "bottleneck")))
})

test_that("the single-bridge design isolates C--D as the only bridge", {
    a <- mgcca_audit(ad_bridge(), L = 2)

    expect_equal(a$graph$n_components, 1L)
    expect_equal(nrow(a$graph$edges), 4L)
    expect_equal(nrow(a$graph$bridges), 1L)
    expect_equal(a$graph$bridges$block1, "C")
    expect_equal(a$graph$bridges$block2, "D")
    expect_equal(unname(a$graph$articulation), c(FALSE, FALSE, TRUE, FALSE))

    expect_equal(pair_of(a, "A", "B", "kappa"), 2L)
    expect_equal(pair_of(a, "A", "C", "kappa"), 2L)
    expect_equal(pair_of(a, "B", "C", "kappa"), 2L)
    expect_equal(pair_of(a, "C", "D", "kappa"), 1L)
    expect_equal(pair_of(a, "A", "D", "kappa"), 1L)
    expect_equal(pair_of(a, "B", "D", "kappa"), 1L)

    expect_equal(pair_of(a, "A", "D", "code"), "MODEL_TRANSFER_ELIGIBLE")
    expect_true(pair_of(a, "A", "D", "single_bridge"))
    expect_false(pair_of(a, "A", "D", "alternative_paths"))
    expect_true(pair_of(a, "A", "B", "alternative_paths"))
})

test_that("a thin joint cell with no route is UNSUPPORTED, not a structural zero", {
    a <- mgcca_audit(ad_thin(), L = 2)

    expect_equal(pair_of(a, "A", "B", "n_jk"), 4L)
    expect_equal(pair_of(a, "A", "B", "pi_hat_jk"), 0.02)
    expect_equal(pair_of(a, "A", "B", "code"), "UNSUPPORTED")
    expect_equal(a$graph$n_components, 2L)
    expect_equal(nrow(a$graph$edges), 0L)

    # and it becomes an edge the moment the declared line is moved below it --
    # the line is declared, it is not a property of the data
    b <- mgcca_audit(ad_thin(), L = 2, design_lines = list(availability_min = 0.01))
    expect_equal(pair_of(b, "A", "B", "code"), "DIRECT_SUPPORT")
})


# ---- (2) provenance: widest routes and bottlenecks --------------------------

test_that("the chain's transfer route and bottleneck are the hand-computed ones", {
    a <- mgcca_audit(ad_chain(), L = 2)

    expect_equal(pair_of(a, "A1", "C", "path"), "A1 -> B -> C")
    expect_equal(pair_of(a, "A1", "C", "path_length"), 2L)
    expect_equal(pair_of(a, "A1", "C", "bottleneck"), 40)     # min(40, 60)
    expect_equal(pair_of(a, "A1", "C", "bottleneck_edge"), "A1--B")

    ed <- a$provenance$paths[["A1--C"]]$edges
    expect_equal(ed$from, c("A1", "B"))
    expect_equal(ed$to, c("B", "C"))
    expect_equal(ed$s_e, c(40, 60))
    expect_true(all(ed$is_bridge))
})

test_that("the cycle's widest routes are the unique hand-computed ones", {
    a <- mgcca_audit(ad_cycle(), L = 2)

    expect_equal(pair_of(a, "A", "C", "path"), "A -> B -> C")
    expect_equal(pair_of(a, "A", "C", "bottleneck"), 20)
    expect_equal(pair_of(a, "B", "D", "path"), "B -> C -> D")
    expect_equal(pair_of(a, "B", "D", "bottleneck"), 30)

    # the widest route between A and D goes the long way round: its weakest
    # edge, 20, beats the direct edge's 10.
    expect_equal(pair_of(a, "A", "D", "n_jk"), 10L)
    expect_equal(pair_of(a, "A", "D", "path"), "A -> B -> C -> D")
    expect_equal(pair_of(a, "A", "D", "bottleneck"), 20)
    expect_equal(pair_of(a, "A", "D", "s_e"), 10)        # s_e stays the direct count
})

test_that("the single-bridge routes bottleneck on the bridge", {
    a <- mgcca_audit(ad_bridge(), L = 2)

    expect_equal(pair_of(a, "A", "D", "path"), "A -> C -> D")
    expect_equal(pair_of(a, "A", "D", "bottleneck"), 20)
    expect_equal(pair_of(a, "A", "D", "bottleneck_edge"), "C--D")
    expect_equal(pair_of(a, "B", "D", "path"), "B -> C -> D")
    expect_equal(pair_of(a, "B", "D", "bottleneck"), 20)
})

test_that("support against L is reported as description, never as a proof", {
    a <- mgcca_audit(ad_chain(), L = 2)
    expect_true(all(a$provenance$pairs$L == 2L))
    expect_equal(pair_of(a, "A1", "C", "n_jk_ge_L"), FALSE)      # n_jk = 0
    expect_equal(pair_of(a, "A1", "C", "bottleneck_ge_L"), TRUE) # B = 40
    expect_match(a$provenance$rank_assumption, "descriptive support only")
    expect_match(a$provenance$rank_assumption,
                 "not empirically checked by this audit")
})


# ---- (3) representation ------------------------------------------------------

test_that("Kish ESS, max weight and proportions are the hand-computed values", {
    # Two strata of four, both blocks complete. Stratum A weights (1,1,1,5):
    #   ESS = 8^2 / (1+1+1+25) = 64/28. Stratum B weights (2,2,2,2): ESS = 4.
    av <- matrix(TRUE, 8, 2, dimnames = list(NULL, c("A", "B")))
    a <- mgcca_audit(av, L = 1, strata = rep(c("s1", "s2"), each = 4),
                     weights = c(1, 1, 1, 5, 2, 2, 2, 2))
    r <- a$representation

    expect_equal(nrow(r), 4L)
    expect_equal(r$proportion, rep(1, 4))
    expect_equal(r$n_obs, rep(4L, 4))
    expect_equal(r$ess[r$stratum == "s1"], rep(64 / 28, 2))
    expect_equal(r$ess[r$stratum == "s2"], rep(4, 2))
    expect_equal(r$max_weight[r$stratum == "s1"], rep(5, 2))
    expect_equal(r$max_weight[r$stratum == "s2"], rep(2, 2))
    # concentration = share of the nominal size lost to the weights' spread
    expect_equal(r$weight_concentration[r$stratum == "s1"], rep(1 - (64 / 28) / 4, 2))
    expect_equal(r$weight_concentration[r$stratum == "s2"], rep(0, 2))

    # the ESS-based pair support is ADDITIONAL: the raw count is untouched
    expect_equal(pair_of(a, "A", "B", "s_e"), 8)
    # all eight are in both blocks: (1+1+1+5+2+2+2+2)^2 / (3*1 + 25 + 4*4)
    expect_equal(pair_of(a, "A", "B", "s_e_ess"), 16^2 / 44)
})

test_that("ESS is absent, not zero, when no weights are declared", {
    av <- matrix(TRUE, 8, 2, dimnames = list(NULL, c("A", "B")))
    a <- mgcca_audit(av, L = 1, strata = rep(c("s1", "s2"), each = 4))

    expect_false(a$weights_declared)
    expect_true(all(is.na(a$representation$ess)))
    expect_true(all(is.na(a$representation$max_weight)))
    expect_true(all(is.na(a$representation$weight_concentration)))
    expect_true(all(is.na(a$provenance$pairs$s_e_ess)))
    expect_false(any(a$representation$below_ess_min))
    expect_false(any(a$representation$above_max_weight))
})

test_that("the availability tables carry the hand-computed counts", {
    a <- mgcca_audit(ad_bridge(), L = 2,
                     strata = rep(c("s1", "s2"), c(30, 20)))

    ab <- a$availability$by_block
    expect_equal(ab$n_j, c(30L, 30L, 50L, 20L))
    expect_equal(ab$pi_hat_j, c(0.6, 0.6, 1, 0.4))

    bs <- a$availability$by_stratum
    expect_equal(bs$n_jc[bs$stratum == "s1"], c(30L, 30L, 30L, 0L))
    expect_equal(bs$n_jc[bs$stratum == "s2"], c(0L, 0L, 20L, 20L))

    # per-stratum joint counts appear only when strata are declared
    expect_false(is.null(a$availability$pairs_by_stratum))
    expect_null(mgcca_audit(ad_bridge(), L = 2)$availability$pairs_by_stratum)
})

test_that("a declared stratum that empties a joint cell removes the edge", {
    # A and B are both complete in stratum s1 and A is absent from s2. Pooled,
    # the pair looks fine; within s2 the joint cell is empty.
    m <- matrix(FALSE, 100, 2, dimnames = list(NULL, c("A", "B")))
    m[1:50, c(1, 2)] <- TRUE
    m[51:100, 2] <- TRUE
    st <- rep(c("s1", "s2"), each = 50)

    pooled <- mgcca_audit(m, L = 2)
    expect_equal(pair_of(pooled, "A", "B", "code"), "DIRECT_SUPPORT")

    a <- mgcca_audit(m, L = 2, strata = st)
    expect_equal(pair_of(a, "A", "B", "n_jk"), 50L)     # not a structural zero
    expect_equal(pair_of(a, "A", "B", "n_jkc_min"), 0L)
    expect_equal(pair_of(a, "A", "B", "code"), "UNSUPPORTED")
})


# ---- (4) agreement with the prototype classifier ----------------------------
#
# The classifier this module replaces lives in the salvage prototype
# (mgcca_new_method/prototype/v1_method.R, v1_classify_graph): an edge wherever
# the joint availability clears 0.05 in EVERY declared stratum, then "direct" /
# "model_assisted" (same component) / "unidentified" (different components). It
# is reimplemented here rather than sourced, so the two implementations are
# genuinely independent, and the designs are regenerated from the prototype's
# own frozen RNG recipe (generator.R: sub_seed, gen_replicate substreams 2 and
# 4, draw_masks_disconnected, the `starved` mechanism).
#
# The frozen vocabulary splits the prototype's "unidentified" in two, by whether
# the joint cell is empty anywhere at all, so the agreement is checked under
# that documented mapping.

proto_components <- function(adj) {
    J <- nrow(adj); comp <- integer(J); cur <- 0L
    for (v in seq_len(J)) {
        if (comp[v] > 0L) next
        cur <- cur + 1L; stack <- v
        while (length(stack)) {
            u <- stack[1L]; stack <- stack[-1L]
            if (comp[u] > 0L) next
            comp[u] <- cur
            stack <- c(stack, setdiff(which(adj[u, ]), which(comp > 0L)))
        }
    }
    comp
}

proto_classify <- function(M, strata, eps_pos = 0.05) {
    J <- ncol(M); slev <- sort(unique(strata))
    pi_pair <- matrix(NA_real_, J, J)
    for (j in seq_len(J)) for (k in j:J) {
        both <- as.integer(M[, j] == 1L & M[, k] == 1L)
        pi_pair[j, k] <- pi_pair[k, j] <-
            min(vapply(slev, function(s) mean(both[strata == s]), 0))
    }
    adj <- pi_pair >= eps_pos
    diag(adj) <- TRUE
    comp <- proto_components(adj)
    out <- character(0)
    for (j in seq_len(J - 1L)) for (k in (j + 1L):J)
        out <- c(out, if (adj[j, k]) "direct"
                      else if (comp[j] == comp[k]) "model_assisted"
                      else "unidentified")
    list(class = out, n_components = max(comp))
}

frozen_to_proto <- function(code)
    c(DIRECT_SUPPORT = "direct",
      MODEL_TRANSFER_ELIGIBLE = "model_assisted",
      UNSUPPORTED = "unidentified",
      STRUCTURAL_ZERO = "unidentified")[code]

# generator.R: sub_seed(master, rep, stream, tag)
proto_seed <- function(master, rep, stream, tag = 0L)
    as.integer(((as.numeric(master) %% 1e5) + 7919 * rep + 104729 * stream +
                    1299709 * tag) %% 2147483647)

# gen_replicate() substream 2: n = 500, p_B = 0.30 exact, balanced by design
proto_g <- function(rep, n = 500L) {
    set.seed(proto_seed(20260906L, rep, 2L, 1L))
    g <- rep("A", n)
    g[sort(sample.int(n, as.integer(round(0.30 * n))))] <- "B"
    g
}

# KT-3a, draw_masks_disconnected(): patterns {1,2} and {3} only, frac12 = 0.5
proto_mask_kt3a <- function(rep, n = 500L, J = 3L) {
    set.seed(proto_seed(20260906L, rep, 4L, 1L))
    M <- matrix(0L, n, J)
    idx <- sample.int(n)
    n12 <- round(0.5 * n)
    M[idx[seq_len(n12)], c(1L, 2L)] <- 1L
    M[idx[(n12 + 1L):n], 3L] <- 1L
    dimnames(M) <- list(NULL, c("b1", "b2", "b3"))
    M
}

# KT-3b, the `starved` mechanism: w = exp(4 * 1{g = B}), pi_3 = 0.45, whole-row
# removal at a fixed count, then the zero-block repair.
proto_mask_kt3b <- function(rep, n = 500L, J = 3L) {
    g <- proto_g(rep, n)
    set.seed(proto_seed(20260906L, rep, 3L, 1L))
    invisible(rep_len(seq_len(6L), n)[sample.int(n)])   # substream 3 (k), unused
    set.seed(proto_seed(20260906L, rep, 4L, 1L))
    w <- exp(4.0 * as.numeric(g == "B"))
    piv <- c(0, 0, 0.45)
    M <- matrix(1L, n, J)
    for (j in seq_len(J)) {
        m <- round(piv[j] * n)
        if (m > 0) M[sample.int(n, m, prob = w), j] <- 0L
    }
    empty <- which(rowSums(M) == 0L)
    if (length(empty)) for (i in empty) M[i, sample.int(J, 1L)] <- 1L
    dimnames(M) <- list(NULL, c("b1", "b2", "b3"))
    M
}

test_that("the classifier agrees with the prototype on all 100 KT-3a designs", {
    ours <- proto <- character(0)
    ncomp_ours <- ncomp_proto <- integer(0)
    for (rep in 1:100) {
        g <- proto_g(rep)
        M <- proto_mask_kt3a(rep)
        a <- mgcca_audit(M, L = 2, strata = g)
        p <- proto_classify(M, g)
        ours <- c(ours, unname(frozen_to_proto(a$provenance$pairs$code)))
        proto <- c(proto, p$class)
        ncomp_ours <- c(ncomp_ours, a$graph$n_components)
        ncomp_proto <- c(ncomp_proto, p$n_components)
    }
    expect_equal(length(ours), 300L)
    expect_identical(ours, proto)
    expect_identical(ncomp_ours, ncomp_proto)

    # and the KT-3a design really is the one the prototype ran on: patterns
    # {1,2} and {3} only, so blocks 1 and 2 are direct and block 3 is cut off
    a <- mgcca_audit(proto_mask_kt3a(1L), L = 2, strata = proto_g(1L))
    expect_equal(pair_of(a, "b1", "b2", "code"), "DIRECT_SUPPORT")
    expect_equal(pair_of(a, "b1", "b3", "code"), "STRUCTURAL_ZERO")
    expect_equal(pair_of(a, "b2", "b3", "code"), "STRUCTURAL_ZERO")
    expect_equal(pair_of(a, "b1", "b2", "n_jk"), 250L)
})

test_that("it also agrees on the KT-3b starved designs, where the cell is thin", {
    ours <- proto <- character(0)
    for (rep in 1:20) {
        g <- proto_g(rep)
        M <- proto_mask_kt3b(rep)
        a <- mgcca_audit(M, L = 2, strata = g)
        ours <- c(ours, unname(frozen_to_proto(a$provenance$pairs$code)))
        proto <- c(proto, proto_classify(M, g)$class)
    }
    expect_identical(ours, proto)

    # this design is the interesting half of the mapping: block 3 IS jointly
    # observed with the others overall, but not once inside stratum B, so the
    # pair is UNSUPPORTED rather than a structural zero.
    a <- mgcca_audit(proto_mask_kt3b(1L), L = 2, strata = proto_g(1L))
    expect_gt(pair_of(a, "b1", "b3", "n_jk"), 0L)
    expect_equal(pair_of(a, "b1", "b3", "n_jkc_min"), 0L)
    expect_equal(pair_of(a, "b1", "b3", "code"), "UNSUPPORTED")
})


# ---- (5) contract fields ----------------------------------------------------

test_that("structural_checks is present and says exactly not_evaluated", {
    for (av in list(ad_chain(), ad_cycle(), ad_disconnected(), ad_bridge()))
        expect_identical(mgcca_audit(av, L = 2)$structural_checks,
                         "not_evaluated")
})

test_that("only the frozen vocabulary is ever emitted", {
    codes <- c("DIRECT_SUPPORT", "MODEL_TRANSFER_ELIGIBLE", "UNSUPPORTED",
               "STRUCTURAL_ZERO")
    advis <- c("LOW_SUPPORT", "SINGLE_BRIDGE", "SEVERE_REPRESENTATION_WARNING")
    seen_codes <- seen_adv <- character(0)
    for (av in list(ad_chain(), ad_cycle(), ad_disconnected(), ad_bridge(),
                    ad_thin())) {
        a <- mgcca_audit(av, L = 2, weights = rep(c(1, 2), length.out = nrow(av)))
        seen_codes <- c(seen_codes, a$status$pairs$code)
        seen_adv <- c(seen_adv,
                      unlist(strsplit(c(a$status$pairs$advisory,
                                        a$status$strata$advisory), ";")))
    }
    expect_true(all(seen_codes %in% codes))
    expect_true(all(seen_adv[nzchar(seen_adv)] %in% advis))
    expect_setequal(unique(seen_codes), codes)

    # no IDENTIFIED_* code exists anywhere in what the object emits
    expect_false(any(grepl("IDENTIFIED", c(seen_codes, seen_adv))))
    expect_false(any(grepl("NOT_IDENTIFIED",
                           capture.output(print(mgcca_audit(ad_thin(), L = 2))))))
})

test_that("the audit never refuses, whatever the design", {
    # nothing observed at all
    empty <- matrix(FALSE, 10, 3, dimnames = list(NULL, c("A", "B", "C")))
    a <- expect_no_error(mgcca_audit(empty, L = 5))
    expect_s3_class(a, "mgcca_audit")
    expect_true(all(a$provenance$pairs$code == "STRUCTURAL_ZERO"))
    expect_equal(a$graph$n_components, 3L)
    expect_true(all(a$availability$by_block$n_j == 0L))

    # one individual, everything present
    one <- matrix(TRUE, 1, 2, dimnames = list(NULL, c("A", "B")))
    b <- expect_no_error(mgcca_audit(one, L = 1))
    expect_equal(pair_of(b, "A", "B", "code"), "DIRECT_SUPPORT")
    expect_true(b$status$pairs$low_support)          # n_jk = 1 < ess_min = 30

    # every declared weight zero: ESS is 0, and it is still reported, not refused
    c0 <- expect_no_error(mgcca_audit(matrix(TRUE, 6, 2,
                                             dimnames = list(NULL, c("A", "B"))),
                                      L = 1, weights = rep(0, 6)))
    expect_equal(c0$representation$ess, rep(0, 2))
    expect_true(all(grepl("SEVERE_REPRESENTATION_WARNING",
                          c0$representation$advisory)))

    # and the whole family still prints and summarises
    expect_output(print(a), "structural_checks: not_evaluated")
    expect_output(summary(c0), "evidences positivity")
})

test_that("the printed output keeps the licensed wording", {
    out <- capture.output(print(mgcca_audit(ad_chain(), L = 2)))
    expect_match(paste(out, collapse = " "), "design-support fragility")
    expect_match(paste(out, collapse = " "),
                 "identified\\s+conditional on the declared rank-L model")
    expect_match(paste(out, collapse = " "), "declared lines with margin")
    expect_match(paste(out, collapse = " "), "refuse")

    s <- capture.output(summary(mgcca_audit(ad_chain(), L = 2)))
    expect_match(paste(s, collapse = " "), "does not prove numerical rank")
    expect_match(paste(s, collapse = " "),
                 "evidences positivity for that cell")
    expect_match(paste(s, collapse = " "),
                 "design-support fragility, not identification robustness")
})

test_that("bad inputs are refused as inputs, with a message that says why", {
    av <- ad_chain()
    expect_error(mgcca_audit(av), "declared rank")
    expect_error(mgcca_audit(av, L = 0), "declared rank")
    expect_error(mgcca_audit(av, L = 2.5), "declared rank")
    expect_error(mgcca_audit(av, L = 2, strata = rep(1, 3)), "one entry per individual")
    expect_error(mgcca_audit(av, L = 2, weights = rep(1, 3)), "one entry per individual")
    expect_error(mgcca_audit(av, L = 2, weights = rep(-1, 100)), "non-negative")
    expect_error(mgcca_audit(av, L = 2, design_lines = list(nope = 1)), "unknown design line")
    expect_error(mgcca_audit(av, L = 2, design_lines = list(ess_min = -1)), "non-negative")
    expect_error(mgcca_audit(av[, 1, drop = FALSE], L = 2), "at least two blocks")
    bad <- av; bad[1, 1] <- NA
    expect_error(mgcca_audit(bad, L = 2), "must not contain NA")
    expect_error(mgcca_audit(matrix(2, 4, 2), L = 2), "logical or 0/1")
})


# ---- (6) edge cases ---------------------------------------------------------

test_that("J = 2 works and has exactly one pair", {
    av <- matrix(TRUE, 40, 2, dimnames = list(NULL, c("A", "B")))
    a <- mgcca_audit(av, L = 2)
    expect_equal(nrow(a$provenance$pairs), 1L)
    expect_equal(a$graph$n_components, 1L)
    expect_equal(pair_of(a, "A", "B", "kappa"), 1L)     # one edge is a bridge
    expect_true(pair_of(a, "A", "B", "single_bridge"))
    expect_true(a$graph$complete)
})

test_that("an all-complete design is DIRECT_SUPPORT throughout", {
    av <- matrix(TRUE, 40, 3, dimnames = list(NULL, c("A", "B", "C")))
    a <- mgcca_audit(av, L = 2)
    expect_true(a$graph$complete)
    expect_true(all(a$provenance$pairs$code == "DIRECT_SUPPORT"))
    expect_true(all(a$provenance$pairs$kappa == 2L))    # a triangle
    expect_equal(nrow(a$graph$bridges), 0L)
    expect_false(any(nzchar(a$status$pairs$advisory)))
    expect_false(any(nzchar(a$status$strata$advisory)))
})

test_that("an empty overlap is a structural zero", {
    av <- matrix(FALSE, 100, 2, dimnames = list(NULL, c("A", "B")))
    av[1:50, 1] <- TRUE
    av[51:100, 2] <- TRUE
    a <- mgcca_audit(av, L = 2)
    expect_equal(pair_of(a, "A", "B", "n_jk"), 0L)
    expect_equal(pair_of(a, "A", "B", "code"), "STRUCTURAL_ZERO")
    expect_equal(a$graph$n_components, 2L)
})

test_that("a single declared stratum reproduces the undeclared case", {
    av <- ad_bridge()
    a <- mgcca_audit(av, L = 2)
    b <- mgcca_audit(av, L = 2, strata = rep("only", nrow(av)))
    expect_identical(a$provenance$pairs, b$provenance$pairs)
    expect_identical(a$graph$kappa, b$graph$kappa)
    expect_equal(a$representation$proportion, b$representation$proportion)
    expect_equal(a$strata_levels, "(all)")
    expect_equal(b$strata_levels, "only")
    expect_false(a$strata_declared)
    expect_true(b$strata_declared)
})

test_that("degenerate weights are described, not hidden", {
    av <- matrix(TRUE, 6, 2, dimnames = list(NULL, c("A", "B")))
    # one individual carries all the weight: ESS falls to nearly 1
    w <- c(1000, 1, 1, 1, 1, 1)
    a <- mgcca_audit(av, L = 1, weights = w)
    expect_equal(a$representation$ess, rep(sum(w)^2 / sum(w^2), 2))
    expect_lt(a$representation$ess[1], 1.1)
    expect_equal(a$representation$max_weight, rep(1000, 2))
    expect_true(all(a$representation$above_max_weight))
    expect_true(all(grepl("SEVERE_REPRESENTATION_WARNING",
                          a$representation$advisory)))

    # equal weights: ESS is exactly the count, concentration exactly zero
    b <- mgcca_audit(av, L = 1, weights = rep(7, 6))
    expect_equal(b$representation$ess, rep(6, 2))
    expect_equal(b$representation$weight_concentration, rep(0, 2))
})


# ---- (7) input polymorphism -------------------------------------------------

test_that("availability from blocks, from HDF5 and by hand are identical", {
    ids <- sprintf("ind%03d", 1:60)
    blocks <- list(
        meth = matrix(rnorm(30 * 4), 30, 4,
                      dimnames = list(ids[1:30], paste0("m", 1:4))),
        expr = matrix(rnorm(60 * 3), 60, 3,
                      dimnames = list(ids, paste0("e", 1:3))),
        prot = matrix(rnorm(30 * 2), 30, 2,
                      dimnames = list(ids[31:60], paste0("p", 1:2))))

    hand <- matrix(FALSE, 60, 3,
                   dimnames = list(ids, c("meth", "expr", "prot")))
    hand[1:30, 1] <- TRUE
    hand[, 2] <- TRUE
    hand[31:60, 3] <- TRUE

    h5 <- tempfile(fileext = ".h5")
    on.exit({
        try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
        unlink(h5)
    }, add = TRUE)
    desc <- mgcca_import_hdf5(blocks, filename = h5, overwriteFile = TRUE)
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)

    from_mem <- mgcca_audit(blocks, L = 2)
    from_h5 <- mgcca_audit(desc, L = 2)
    from_path <- mgcca_audit(h5, L = 2)
    from_hand <- mgcca_audit(hand, L = 2)

    expect_identical(from_mem$table, hand)
    expect_identical(from_h5$table, hand)
    expect_identical(from_mem$individuals, ids)

    # A BARE PATH has no declared dataset order, so the file's own listing
    # decides it -- exactly as it does for mgcca() itself. Same design, blocks
    # in the file's order; pass the descriptor to fix the order.
    expect_identical(from_path$table, hand[, sort(colnames(hand))])
    expect_identical(from_path$individuals, ids)

    # and everything the audit derives from them agrees too
    expect_identical(from_mem$provenance$pairs, from_hand$provenance$pairs)
    expect_identical(from_h5$provenance$pairs, from_hand$provenance$pairs)
    expect_identical(from_h5$graph$kappa, from_hand$graph$kappa)
    expect_equal(from_mem$source, "in-memory blocks")
    expect_equal(from_h5$source, "HDF5")

    # the design itself: meth and prot share nobody, but both meet expr
    expect_equal(pair_of(from_h5, "meth", "prot", "code"),
                 "MODEL_TRANSFER_ELIGIBLE")
    expect_equal(pair_of(from_h5, "meth", "prot", "path"),
                 "meth -> expr -> prot")
})

test_that("blocks without individual IDs are refused, as mgcca refuses them", {
    blocks <- list(a = matrix(1, 4, 2), b = matrix(1, 4, 2))
    expect_error(mgcca_audit(blocks, L = 1), "no rownames")
    expect_error(mgcca_audit(42, L = 1), "unsupported input type")
    expect_error(mgcca_audit("no-such-file.h5", L = 1), "does not exist")
})


# ---- (8) round trip through HDF5 --------------------------------------------

test_that("a saved audit comes back whole", {
    ids <- sprintf("ind%03d", 1:60)
    blocks <- list(
        meth = matrix(rnorm(30 * 4), 30, 4,
                      dimnames = list(ids[1:30], paste0("m", 1:4))),
        expr = matrix(rnorm(60 * 3), 60, 3,
                      dimnames = list(ids, paste0("e", 1:3))),
        prot = matrix(rnorm(30 * 2), 30, 2,
                      dimnames = list(ids[31:60], paste0("p", 1:2))))
    h5 <- tempfile(fileext = ".h5")
    on.exit({
        try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
        unlink(h5)
    }, add = TRUE)
    desc <- mgcca_import_hdf5(blocks, filename = h5, overwriteFile = TRUE)
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)

    # a fixture that exercises every component: a transfer route with its
    # bottleneck, advisory codes, declared strata AND declared weights
    st <- rep(c("s1", "s2"), each = 30)
    w <- rep(c(1, 2, 3), each = 20)
    a <- mgcca_audit(desc, L = 2, strata = st, weights = w, save_hdf5 = TRUE)

    expect_false(is.null(a$saved))
    expect_equal(a$saved$group, "availability_audit")
    expect_true(all(c("table", "pairs", "paths", "representation",
                      "pair_codes", "advisory_codes", "strata", "weights",
                      "adjacency", "kappa") %in% a$saved$tables))

    b <- mgcca_audit_load(h5)
    expect_s3_class(b, "mgcca_audit")
    expect_identical(a, b)

    # every combination of the two optional inputs, since each changes which
    # columns carry values and which carry NA
    for (cs in list(list(strata = NULL, weights = NULL, L = 3L),
                    list(strata = rep(c("a", "b", "c"), each = 20),
                         weights = NULL, L = 1L),
                    list(strata = NULL, weights = c(rep(1, 59), 500), L = 2L))) {
        aa <- mgcca_audit(desc, L = cs$L, strata = cs$strata,
                          weights = cs$weights, save_hdf5 = TRUE)
        expect_identical(aa, mgcca_audit_load(h5))
    }

    a2 <- mgcca_audit(desc, L = 3, save_hdf5 = TRUE)
    b2 <- mgcca_audit_load(h5)
    expect_identical(a2, b2)
    expect_false(b2$strata_declared)
    expect_false(b2$weights_declared)
    expect_identical(mgcca_audit_load(h5, check_inputs = FALSE), b2)

    # what is on disk is inspectable, not a blob: the legend names the codes
    leg <- BigDataStatMeth::hdf5_matrix(h5, "availability_audit/pair_codes")
    lrn <- rownames(leg)
    if (is.function(leg$close)) leg$close()
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
    expect_equal(lrn, c("DIRECT_SUPPORT", "MODEL_TRANSFER_ELIGIBLE",
                        "UNSUPPORTED", "STRUCTURAL_ZERO"))

    # the routes are stored as sequences of integer block indices
    pth <- BigDataStatMeth::hdf5_matrix(h5, "availability_audit/paths")
    pm <- as.matrix(pth)
    if (is.function(pth$close)) pth$close()
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
    expect_equal(nrow(pm), nrow(a2$provenance$pairs))
    expect_true(all(pm == round(pm)))
    expect_true(all(pm >= 0 & pm <= a2$J))
    # pairs are in combn() order, so row 2 is meth--prot: meth -> expr -> prot
    expect_equal(unname(pm[2, seq_len(3)]), c(1, 2, 3))
})

test_that("saving needs somewhere to save to, and loading needs an audit there", {
    expect_error(mgcca_audit(ad_chain(), L = 2, save_hdf5 = TRUE),
                 "needs an HDF5-backed design")
    expect_error(mgcca_audit(ad_chain(), L = 2, save_hdf5 = "yes"),
                 "must be TRUE or FALSE")

    h5 <- tempfile(fileext = ".h5")
    on.exit({
        try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
        unlink(h5)
    }, add = TRUE)
    hm <- BigDataStatMeth::hdf5_create_matrix(
        h5, "MGCCA_IN/x", data = matrix(1, 2, 2,
                                        dimnames = list(c("a", "b"), c("u", "v"))),
        dtype = "double", overwrite = TRUE)
    if (is.function(hm$close)) hm$close()
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)

    expect_error(mgcca_audit_load(h5), "no mgcca audit is stored")
    expect_error(mgcca_audit_load("no-such-file.h5"), "does not exist")
})

test_that("the loader reads what is stored rather than recomputing it", {
    # The proof: change a stored number and the loaded object must change with
    # it -- and say so. A loader that recomputed from the stored inputs would
    # hand back the original value and notice nothing.
    ids <- sprintf("ind%03d", 1:40)
    blocks <- list(
        a = matrix(rnorm(20 * 2), 20, 2,
                   dimnames = list(ids[1:20], c("a1", "a2"))),
        b = matrix(rnorm(40 * 2), 40, 2,
                   dimnames = list(ids, c("b1", "b2"))),
        c = matrix(rnorm(20 * 2), 20, 2,
                   dimnames = list(ids[21:40], c("c1", "c2"))))
    h5 <- tempfile(fileext = ".h5")
    on.exit({
        try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
        unlink(h5)
    }, add = TRUE)
    desc <- mgcca_import_hdf5(blocks, filename = h5, overwriteFile = TRUE)
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
    a <- mgcca_audit(desc, L = 2, save_hdf5 = TRUE)

    hm <- BigDataStatMeth::hdf5_matrix(h5, "availability_audit/pairs")
    pm <- as.matrix(hm)
    cn <- colnames(hm)
    if (is.function(hm$close)) hm$close()
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)
    colnames(pm) <- cn
    pm[1L, "n_jk"] <- 12345
    hm <- BigDataStatMeth::hdf5_create_matrix(h5, "availability_audit/pairs",
                                              data = pm, dtype = "double",
                                              overwrite = TRUE)
    if (is.function(hm$close)) hm$close()
    try(BigDataStatMeth::hdf5_close_all(), silent = TRUE)

    expect_warning(b <- mgcca_audit_load(h5), "does not match what its own")
    expect_equal(b$provenance$pairs$n_jk[1L], 12345L)
    expect_false(identical(a$provenance$pairs, b$provenance$pairs))

    # and with the cross-check switched off it is silent, but still the stored
    # value that comes back
    expect_no_warning(d <- mgcca_audit_load(h5, check_inputs = FALSE))
    expect_equal(d$provenance$pairs$n_jk[1L], 12345L)
})

# The module's own code, with string literals blanked out: a package name
# mentioned inside an error message is not a call to it, and these two tests are
# about what the code CALLS.
audit_module_code <- function() {
    ns <- asNamespace("mgcca")
    nms <- c("mgcca_audit", "mgcca_audit_load",
             grep("^\\.mgcca_av_", ls(ns, all.names = TRUE), value = TRUE),
             "print.mgcca_audit", "summary.mgcca_audit", "plot.mgcca_audit")
    src <- unlist(lapply(nms, function(n) deparse(get(n, envir = ns))))
    gsub('"[^"]*"', '""', src)
}

test_that("nothing in the module is serialised", {
    src <- audit_module_code()
    for (bad in c("serialize(", "unserialize(", "saveRDS(", "readRDS(",
                  "rawToChar(", "charToRaw(", "base::serialize"))
        expect_false(any(grepl(bad, src, fixed = TRUE)))
})

test_that("the module reads and writes HDF5 through BigDataStatMeth, never rhdf5", {
    src <- audit_module_code()
    expect_false(any(grepl("rhdf5", src, fixed = TRUE)))
    expect_true(any(grepl("BigDataStatMeth::hdf5_matrix", src, fixed = TRUE)))
    expect_true(any(grepl("BigDataStatMeth::hdf5_create_matrix", src,
                          fixed = TRUE)))
})


# ---- (9) the visual layer ----------------------------------------------------
#
# The frozen visual rules V1-V7 are design decisions, not incidental styling, so
# they are pinned here: a later tweak that reintroduces a traffic light, or that
# lets edge width mean something other than co-observation, has to break a test
# to get in.

# The vignette's design: four blocks, one of them reachable only across a single
# bridge, with declared strata and declared weights.
ad_vignette <- function() {
    ids <- sprintf("ind%02d", 1:85)
    av <- cbind(methylation = seq_along(ids) %in% 1:60,
                expression  = seq_along(ids) %in% 1:60,
                proteins    = rep(TRUE, length(ids)),
                metabolites = seq_along(ids) %in% 61:85)
    rownames(av) <- ids
    mgcca_audit(av, L = 2,
                strata = rep(c("cohort A", "cohort B"), length.out = 85),
                weights = ifelse(seq_along(ids) %in% 1:60, 1, 3))
}

# the module's internal constants, reached the same way the earlier tests reach
# its functions (no `:::`, which R CMD check would rather we did not use)
av_internal <- function(nm) get(nm, envir = asNamespace("mgcca"))

# perceived luminance, for "which of these is the strongest colour"
lum <- function(hex) {
    v <- grDevices::col2rgb(hex)
    as.numeric(0.299 * v[1, ] + 0.587 * v[2, ] + 0.114 * v[3, ]) / 255
}
hue <- function(hex) {
    h <- grDevices::rgb2hsv(grDevices::col2rgb(hex))
    list(h = 360 * as.numeric(h["h", ]), s = as.numeric(h["s", ]))
}


test_that("every plot type draws, and only overview draws a composite", {
    a <- ad_vignette()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (ty in c("matrix", "graph", "availability", "representation",
                 "patterns")) {
        p <- plot(a, type = ty)
        expect_s3_class(p, "ggplot")
        expect_no_error(print(p))
    }
    expect_no_error(print(plot(a, type = "graph",
                               pair = c("methylation", "metabolites"))))

    # overview draws the composite itself and hands back its panels
    panels <- plot(a)
    expect_type(panels, "list")
    expect_length(panels, 4L)
    expect_named(panels, c("matrix", "graph", "availability",
                           "representation"))
    for (p in panels) expect_s3_class(p, "ggplot")

    # bare plot() IS overview, and it returns invisibly
    expect_invisible(plot(a))
    expect_identical(names(plot(a, type = "overview")), names(panels))
})

test_that("V1: no traffic-light semantics anywhere in the palettes", {
    for (col in c(av_internal(".mgcca_av_class_fill"),
                  av_internal(".mgcca_av_redundancy"),
                  av_internal(".mgcca_av_component_pal"),
                  av_internal(".mgcca_av_stratum_pal"),
                  av_internal(".mgcca_av_ramp_high"),
                  av_internal(".mgcca_av_accent"))) {
        z <- hue(col)
        if (z$s < 0.2) next                     # greys carry no hue
        expect_false(z$h <= 20 || z$h >= 345)   # not red
        expect_false(z$h >= 90 && z$h <= 160)   # not green
    }
})

test_that("V2 and V6: DIRECT is strongest, and the two neutrals differ", {
    fills <- av_internal(".mgcca_av_class_fill")
    l <- lum(fills)
    names(l) <- names(fills)
    # darkest = visually strongest, and that must be DIRECT_SUPPORT
    expect_identical(names(which.min(l)), "DIRECT_SUPPORT")
    expect_lt(l[["DIRECT_SUPPORT"]], l[["MODEL_TRANSFER_ELIGIBLE"]])
    # V6: the two unsupported classes are far apart in lightness, and the
    # structural zero additionally carries its own glyph
    expect_gt(abs(l[["UNSUPPORTED"]] - l[["STRUCTURAL_ZERO"]]), 0.15)

    dis <- mgcca_audit(ad_disconnected(), L = 2)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plot(dis, type = "matrix")
    shapes <- vapply(p$layers, function(l)
        if (is.null(l$aes_params$shape)) NA_integer_
        else as.integer(l$aes_params$shape), integer(1))
    expect_true(4L %in% shapes)                 # the X glyph layer is there
    expect_no_error(print(p))
})

test_that("V3: edge width is co-observation count, and says so", {
    a <- ad_vignette()
    p <- plot(a, type = "graph")
    nm <- vapply(p$scales$scales, function(s)
        if ("linewidth" %in% s$aesthetics) as.character(s$name) else NA_character_,
        character(1))
    expect_true("jointly observed n" %in% nm)
    expect_match(p$labels$caption,
                 "Edge width denotes raw co-observation count only.",
                 fixed = TRUE)
})

test_that("V4: the graph legend says redundancy, not robustness", {
    p <- plot(ad_vignette(), type = "graph")
    nm <- unlist(lapply(p$scales$scales, function(s) as.character(s$name)))
    expect_true("Design-support redundancy" %in% nm)
    txt <- paste(c(unlist(p$labels), nm), collapse = " ")
    expect_false(grepl("robust", txt, ignore.case = TRUE))
    expect_false(grepl("reliab", txt, ignore.case = TRUE))
})

test_that("V5: declared lines are labelled as declared, in a neutral style", {
    a <- ad_vignette()
    p <- plot(a, type = "representation")
    vl <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
    expect_length(vl, 1L)
    expect_identical(vl[[1L]]$aes_params$linetype, "dashed")
    expect_identical(vl[[1L]]$aes_params$colour, "grey45")
    expect_true(any(vapply(p$layers, function(l)
        !is.null(l$data) && is.data.frame(l$data) &&
            "label" %in% names(l$data) &&
            any(l$data$label == "declared design line"), logical(1))))
    expect_match(p$labels$caption, "declared design lines with margin",
                 fixed = TRUE)
    # two aligned facets, never one dual axis
    expect_s3_class(p$facet, "FacetWrap")
    expect_match(p$labels$subtitle, "not comparable", fixed = TRUE)
})

test_that("V7: a provenance overlay always carries the rank-L condition", {
    a <- ad_vignette()
    p <- plot(a, type = "graph", pair = c("methylation", "metabolites"))
    expect_identical(p$labels$subtitle,
                     "Model-transfer-eligible under declared rank-L assumption")
    expect_match(p$labels$caption, "methylation -> proteins -> metabolites",
                 fixed = TRUE)
    expect_match(p$labels$caption, "L = 2 as declared", fixed = TRUE)
    expect_false(grepl("identified path", p$labels$subtitle, fixed = TRUE))

    # when the audit stored alternative routes, the subtitle says so rather
    # than letting one route look privileged
    cyc <- mgcca_audit(ad_cycle(), L = 2)
    p2 <- plot(cyc, type = "graph", pair = c("A", "C"))
    expect_true(pair_of(cyc, "A", "C", "alternative_paths"))
    expect_match(p2$labels$subtitle, "alternative routes exist", fixed = TRUE)
})

test_that("the provenance overlay refuses pairs it cannot draw a route for", {
    a <- ad_vignette()
    # DIRECT: carried by its own participants, so there is no route to draw
    expect_error(plot(a, type = "graph",
                      pair = c("methylation", "expression")),
                 "coded DIRECT_SUPPORT")
    expect_error(plot(a, type = "graph",
                      pair = c("methylation", "expression")),
                 "no transfer route to draw")
    # STRUCTURAL_ZERO and UNSUPPORTED
    expect_error(plot(mgcca_audit(ad_disconnected(), L = 2), type = "graph",
                      pair = c("A", "C")),
                 "coded STRUCTURAL_ZERO")
    expect_error(plot(mgcca_audit(ad_thin(), L = 2), type = "graph",
                      pair = c("A", "B")),
                 "coded UNSUPPORTED")
    # and the plain input errors
    expect_error(plot(a, type = "graph", pair = c("nope", "expression")),
                 "unknown block name")
    expect_error(plot(a, type = "graph", pair = c("proteins", "proteins")),
                 "two different blocks")
    expect_error(plot(a, type = "graph", pair = "proteins"),
                 "must be two block names")
    expect_error(plot(a, type = "matrix", pair = c("methylation", "proteins")),
                 "needs type = \"graph\"")
})

test_that("the matrix prints the support quantity its class calls for", {
    a <- ad_vignette()
    p <- plot(a, type = "matrix")
    txt <- Filter(function(l) inherits(l$geom, "GeomText") &&
                      "code" %in% names(l$data), p$layers)[[1L]]$data
    lab <- stats::setNames(txt$label, paste(txt$code))
    expect_true(all(lab[names(lab) == "DIRECT_SUPPORT"] %in%
                        as.character(a$provenance$pairs$n_jk)))
    expect_true(all(grepl("^B=", lab[names(lab) == "MODEL_TRANSFER_ELIGIBLE"])))

    z <- plot(mgcca_audit(ad_disconnected(), L = 2), type = "matrix")
    zt <- Filter(function(l) inherits(l$geom, "GeomText") &&
                     "code" %in% names(l$data), z$layers)[[1L]]$data
    expect_identical(unique(zt$label[zt$code == "STRUCTURAL_ZERO"]), "0")

    u <- plot(mgcca_audit(ad_thin(), L = 2), type = "matrix")
    ut <- Filter(function(l) inherits(l$geom, "GeomText") &&
                     "code" %in% names(l$data), u$layers)[[1L]]$data
    expect_identical(unique(ut$label[ut$code == "UNSUPPORTED"]), "")
})

test_that("the graph node order is component, then degree, then name", {
    # Columns are deliberately in a different order from the answer. b--c--d is
    # a triangle, e hangs off b alone, a is isolated -- so degrees are
    # b 3, c 2, d 2, e 1, a 0, and c/d tie on degree.
    m <- matrix(FALSE, 75, 5,
                dimnames = list(NULL, c("e", "d", "c", "b", "a")))
    m[1:30, c("b", "c", "d")] <- TRUE
    m[31:55, c("b", "e")] <- TRUE
    m[56:75, "a"] <- TRUE
    a <- mgcca_audit(m, L = 2)
    expect_equal(a$graph$n_components, 2L)

    node_order <- av_internal(".mgcca_av_node_order")
    ord <- node_order(a)
    got <- a$blocks[ord]

    # component first (and components keep their own numbering), then degree
    # descending -- b leads its component although its name does not -- then
    # the name, which is what separates the tied c and d.
    expect_identical(got, c("b", "c", "d", "e", "a"))
    expect_identical(unname(a$graph$component[got]),
                     sort(unname(a$graph$component[got])))
    expect_identical(unname(rowSums(a$graph$adjacency)[got]),
                     c(3, 2, 2, 1, 0))

    # and it is a pure function of the audit: same input, same order
    expect_identical(node_order(a), ord)
})

test_that("pattern ordering is frequency, then lexicographic, deterministically", {
    a <- ad_vignette()
    p <- plot(a, type = "patterns")
    bars <- Filter(function(l) inherits(l$geom, "GeomCol"), p$layers)[[1L]]$data
    # levels run bottom-to-top, so the plotted order top-down is their reverse
    top_down <- rev(levels(bars$pattern))
    expect_identical(top_down, c("1110", "0011"))
    expect_identical(bars$value[match(top_down, as.character(bars$pattern))],
                     c(60L, 25L))

    # an exact tie falls back to the pattern key, ascending
    tie <- rbind(cbind(A = TRUE, B = TRUE, C = FALSE),
                 cbind(A = TRUE, B = TRUE, C = FALSE),
                 cbind(A = FALSE, B = TRUE, C = TRUE),
                 cbind(A = FALSE, B = TRUE, C = TRUE))
    pt <- plot(mgcca_audit(tie, L = 1), type = "patterns")
    tb <- Filter(function(l) inherits(l$geom, "GeomCol"), pt$layers)[[1L]]$data
    expect_identical(rev(levels(tb$pattern)), c("011", "110"))

    # the counts are the availability table's own, aggregated at plot time
    expect_equal(sum(bars$value), a$n)
    expect_null(a[["patterns"]])          # nothing was added to the object
})

test_that("the availability panel keeps its fixed [0, 1] scale and block order", {
    a <- ad_vignette()
    p <- plot(a, type = "availability")
    fs <- Filter(function(s) "fill" %in% s$aesthetics, p$scales$scales)[[1L]]
    expect_equal(fs$limits, c(0, 1))
    expect_identical(levels(p$data$block), a$blocks)
    expect_identical(levels(p$data$stratum), rev(a$strata_levels))

    # the ramp must NOT be the classification blue: the two panels are read
    # side by side in the overview, and a shared hue would invite the reader to
    # connect a pale tile with a DIRECT_SUPPORT cell, which means nothing
    expect_false(identical(av_internal(".mgcca_av_ramp_high"),
                           av_internal(".mgcca_av_class_fill")[["DIRECT_SUPPORT"]]))
})

test_that("each recoloured panel encodes a variable, not decoration", {
    a <- ad_vignette()

    # graph nodes: the connected component
    g <- plot(a, type = "graph")
    fs <- Filter(function(s) "fill" %in% s$aesthetics, g$scales$scales)
    expect_length(fs, 1L)
    expect_identical(as.character(fs[[1L]]$name), "Connected component")
    nodes <- Filter(function(l) inherits(l$geom, "GeomPoint"), g$layers)[[1L]]$data
    expect_identical(levels(nodes$component),
                     as.character(sort(unique(a$graph$component))))

    # two components must come out as two colours
    m <- matrix(FALSE, 60, 3, dimnames = list(NULL, c("A", "B", "C")))
    m[1:30, c("A", "B")] <- TRUE
    m[31:60, "C"] <- TRUE
    g2 <- plot(mgcca_audit(m, L = 2), type = "graph")
    nd <- Filter(function(l) inherits(l$geom, "GeomPoint"), g2$layers)[[1L]]$data
    expect_length(levels(nd$component), 2L)

    # representation bars: the declared stratum, same colours in both facets
    r <- plot(a, type = "representation")
    rs <- Filter(function(s) "fill" %in% s$aesthetics, r$scales$scales)
    expect_length(rs, 1L)
    expect_identical(as.character(rs[[1L]]$name), "Declared stratum")
    expect_identical(levels(r$data$stratum), a$strata_levels)
    for (m2 in levels(r$data$measure))
        expect_identical(as.character(r$data$stratum[r$data$measure == m2]),
                         as.character(r$data$stratum[
                             r$data$measure == levels(r$data$measure)[1L]]))

    # with no strata declared the single-level key is suppressed, not shown
    plain <- plot(mgcca_audit(ad_bridge(), L = 2), type = "representation")
    expect_identical(plain$guides$guides$fill, "none")
})

test_that("the representation panel falls back to counts with no weights", {
    a <- mgcca_audit(ad_bridge(), L = 2)
    expect_false(a$weights_declared)
    p <- plot(a, type = "representation")
    expect_true(any(grepl("no weights declared", levels(p$data$measure))))
    expect_false(any(grepl("^ESS", p$data$label)))
})


test_that("the overview's graph panel keeps its subtitle clear of the legend", {
    a <- ad_vignette()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    panels <- plot(a)
    standalone <- plot(a, type = "graph")

    # the compact panel sits beside a legend column, so it takes the short form
    expect_identical(panels$graph$labels$subtitle,
                     "an edge is a pair supported directly")
    # the standalone keeps the whole sentence
    expect_gt(nchar(standalone$labels$subtitle),
              nchar(panels$graph$labels$subtitle))
    expect_match(standalone$labels$subtitle, "node order is component",
                 fixed = TRUE)
    # and V3's frozen caption is never what gets shortened
    for (p in list(panels$graph, standalone))
        expect_match(p$labels$caption,
                     "Edge width denotes raw co-observation count only.",
                     fixed = TRUE)
})


test_that("every in-panel size is a multiple of base_size, in both variants", {
    a <- ad_vignette()

    # Absolute sizes only: rel() values are multiples of the base size by
    # construction, so they cannot be the ones that drift. A hard-coded point
    # value, wherever it hides, shows up here as a size that fails to double
    # when the base size doubles.
    absolute_sizes <- function(p) {
        lay <- unlist(lapply(p$layers, function(l) l$aes_params$size))
        th <- unlist(lapply(
            c("text", "legend.text", "legend.title", "plot.title",
              "plot.subtitle", "plot.caption", "axis.text", "axis.text.x"),
            function(el) {
                e <- p$theme[[el]]
                if (is.null(e) || is.null(e$size) ||
                    inherits(e$size, "rel")) NULL else e$size
            }))
        as.numeric(c(lay, th))
    }

    for (ty in c("matrix", "graph", "availability", "representation",
                 "patterns")) {
        draw <- av_internal(paste0(".mgcca_av_plot_", ty))
        for (compact in c(FALSE, TRUE)) {
            s1 <- absolute_sizes(draw(a, base_size = 10, compact = compact))
            s2 <- absolute_sizes(draw(a, base_size = 20, compact = compact))
            expect_gt(length(s1), 0L)
            expect_identical(length(s1), length(s2))
            expect_equal(s2, 2 * s1)
        }
    }
})

test_that("the composite grows with base_size instead of ignoring it", {
    a <- ad_vignette()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    # the documented default, and the fraction the four shared panels take
    expect_identical(formals(av_internal("plot.mgcca_audit"))$base_size, 13)
    scale <- av_internal(".mgcca_av_overview_scale")

    small <- plot(a, base_size = 12)
    big <- plot(a, base_size = 24)
    for (nm in names(small)) {
        expect_equal(big[[nm]]$theme$text$size, 24 * scale)
        expect_equal(small[[nm]]$theme$text$size, 12 * scale)
    }
})

test_that("the patterns panel names every block, however many there are", {
    # Regression: the two facets share one x scale, and the break function
    # told them apart by the sign of the upper limit. The right-hand expansion
    # is a fraction of the span, so from six blocks on it pushed the block
    # panel's limit to zero, the panel tested as the count panel, and every
    # block name vanished.
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    block_labels <- function(a) {
        b <- ggplot2::ggplot_build(plot(a, type = "patterns"))
        b$layout$panel_params[[1L]]$x$get_labels()
    }

    expect_identical(block_labels(ad_vignette()), ad_vignette()$blocks)

    n <- 140L
    i <- seq_len(n)
    wide <- cbind(genomics    = i %in% c(1:50, 91:100),
                  methylation = i %in% 1:50,
                  expression  = i %in% c(1:50, 51:90),
                  proteins    = i %in% c(51:90, 91:100),
                  metabolites = i %in% 101:140,
                  microbiome  = i %in% 101:140)
    rownames(wide) <- sprintf("id%03d", i)
    six <- mgcca_audit(wide, L = 2)
    expect_identical(block_labels(six), six$blocks)

    # and the count panel still gets numbers, not block names
    b <- ggplot2::ggplot_build(plot(six, type = "patterns"))
    counts <- b$layout$panel_params[[2L]]$x$get_labels()
    expect_false(any(counts %in% six$blocks))
})

test_that("a single-block pattern is drawn in its place, not at the top", {
    # Regression, found on the real HELIX design (5 blocks, 19 patterns, one
    # single-block pattern of 204 individuals drawn above the pattern of 881).
    # The discrete y scale used to collect its level order layer by layer, and
    # the first layer -- the connector segments -- has no row for a pattern
    # that carries only one block, because there is nothing to connect. Those
    # patterns were therefore appended after every level the segment layer had
    # already trained, i.e. drawn at the TOP of the plot, out of frequency
    # order. The scale's limits are pinned explicitly now.
    #
    # `panel_params$y$get_limits()` returns the levels bottom to top, so the
    # LAST one is the top row and the plotted order is their reverse.
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    y_limits <- function(a) {
        b <- ggplot2::ggplot_build(plot(a, type = "patterns"))
        b$layout$panel_params[[1L]]$y$get_limits()
    }

    # four blocks; two of the patterns carry a single block, and one of those
    # is neither the most nor the least frequent, so a scale that appends it
    # cannot land it in the right row by luck
    n <- 107L
    i <- seq_len(n)
    av <- cbind(A = i %in% c(1:40, 41:65),
                B = i %in% c(1:40, 41:65),
                C = i %in% c(1:40, 66:83, 96:107),
                D = i %in% c(1:40, 84:95, 96:107))
    rownames(av) <- sprintf("id%03d", i)
    a <- mgcca_audit(av, L = 2)

    # 1111 = 40, 1100 = 25, 0010 = 18, 0001 = 12, 0011 = 12: frequency
    # descending, and the exact tie broken by the pattern key ascending
    expect_identical(rev(y_limits(a)),
                     c("1111", "1100", "0010", "0001", "0011"))

    # read off the plot itself: top to bottom, the counts never increase
    bars <- Filter(function(l) inherits(l$geom, "GeomCol"),
                   plot(a, type = "patterns")$layers)[[1L]]$data
    top_down <- rev(y_limits(a))
    counts <- bars$value[match(top_down, as.character(bars$pattern))]
    expect_identical(counts, c(40L, 25L, 18L, 12L, 12L))
    expect_true(all(diff(counts) <= 0L))

    # the same, with six blocks, where the block-panel break rule also has to
    # keep naming every block (the earlier regression in this file)
    m <- 105L
    j <- seq_len(m)
    wide <- cbind(genomics    = j %in% c(1:45, 46:70),
                  methylation = j %in% c(1:45, 46:70, 96:105),
                  expression  = j %in% c(1:45, 46:70),
                  proteins    = j %in% c(1:45, 71:85),
                  metabolites = j %in% c(1:45, 86:95),
                  microbiome  = j %in% c(1:45, 86:95))
    rownames(wide) <- sprintf("id%03d", j)
    six2 <- mgcca_audit(wide, L = 2)

    # 111111 = 45, 111000 = 25, 000100 = 15, 000011 = 10, 010000 = 10
    six_top_down <- rev(y_limits(six2))
    expect_identical(six_top_down,
                     c("111111", "111000", "000100", "000011", "010000"))
    six_bars <- Filter(function(l) inherits(l$geom, "GeomCol"),
                       plot(six2, type = "patterns")$layers)[[1L]]$data
    six_counts <- six_bars$value[match(six_top_down,
                                       as.character(six_bars$pattern))]
    expect_identical(six_counts, c(45L, 25L, 15L, 10L, 10L))
    expect_true(all(diff(six_counts) <= 0L))

    bb <- ggplot2::ggplot_build(plot(six2, type = "patterns"))
    expect_identical(bb$layout$panel_params[[1L]]$x$get_labels(), six2$blocks)
})
