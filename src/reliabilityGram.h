// reliabilityGram.h -- K1: present-only, per-feature-standardized, feature-STREAMED
// participant Gram  G = Z^T Z  (n_pr x n_pr), Z_{p x n_pr} standardized per feature (row)
// over the PRESENT participants only, with the (n_pr - 1) sample-SD denominator.
//
// This is the GENERIC block primitive of the mgcca reliability-layer C++/HDF5 port
// (PORT_FIRST_KERNEL.md, ChatGPT r148). It is written with a clean interface and NO
// mgcca-specific dependencies so it is a candidate to lift into BigDataStatMeth later
// (port-bdsm-reuse-architecture) -- but is NOT upstreamed during K1 certification.
//
// ORACLE (numerical source of truth): reliability/helix/34_methyl_gram_dense.R
//   mu = rowMeans(M); Mc = M - mu; sdv = sqrt(rowSums(Mc^2)/(n-1)); Z = Mc/sdv;
//   G = crossprod(Z) = Z^T Z.  Standardize on the OBSERVED participants only.
//
// FROZEN CONVENTIONS (PORT_FIRST_KERNEL.md):
//  * Logical layout p x n_pr (features x present-participants); the HDF5 dataset is
//    read in R-view (nrows_r = features, ncols_r = all participants). Present columns
//    are selected by ID; absent participants NEVER enter centring/scaling/crossprod.
//  * Zero-variance rule COPIED from the oracle: STOP-and-report on s_f <= var_eps
//    (script 34 l.45 stopifnot(all(sdv > 1e-8))). No VAR_EPS injection, no drop, no
//    rescale, no exclusion mode -> p_eff = p on every successful run. The expected
//    algebraic identity tr(G) = p_eff (n_pr - 1) holds mathematically; its finite-
//    precision residual is RECORDED (e_trace) and checked against a tolerance, NOT
//    asserted bitwise (ChatGPT r150 §4).
//  * Failure is ATOMIC (r150 §5): a degenerate feature throws AFTER some chunks have
//    accumulated, so NO partial Gram is returned; the error reports the ABSOLUTE
//    logical feature index and its id when available.
//  * Design A accumulation (r150 §3): the FULL chunk cross-products are accumulated
//    (both triangles populated), so sym_err = max|G - G^T| is a GENUINE asymmetry of
//    the accumulated matrix; it is recorded BEFORE the deterministic symmetrization.
//    Single, deterministic outer feature-chunk loop; dense within-chunk work; serial
//    accumulation into ONE n_pr x n_pr Gram; NO parallel cross-chunk reduction; NO
//    standardized p x n_pr dataset materialized. Effective Eigen thread count recorded.
//
// Include AFTER BigDataStatMeth.hpp.
#ifndef MGCCA_RELIABILITY_GRAM_HPP
#define MGCCA_RELIABILITY_GRAM_HPP

#include <BigDataStatMeth.hpp>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <stdexcept>

namespace mgcca {
namespace reliability {

// Per-chunk ledger record (independent of the final matrix comparison).
struct ChunkRecord {
    long   f_start;          // logical feature start index (0-based)
    long   f_end;            // logical feature end index (0-based, inclusive)
    long   n_features;       // features read + accepted in this chunk
    double max_abs_row_sum;  // max_f | sum_i Z_{f i} |  (post-standardization; ~0)
    double max_abs_ss_dev;   // max_f | ||Z_f||^2 - (n_pr - 1) |                (~0)
    double trace_contrib;    // sum over this chunk's features of ||Z_f||^2  (~ fb*(n_pr-1))
};

// Completion diagnostics + result.
struct GramResult {
    Eigen::MatrixXd          G;            // n_pr x n_pr, symmetrized
    std::vector<std::string> ids;          // present participant ids, frozen order
    long   p_eff;                          // features actually entering G
    long   n_pr;                           // present participants
    long   N_all;                          // participants in the source dataset
    long   chunk;                          // feature chunk size used
    int    n_chunks;
    long   last_chunk;                     // size of the final (possibly partial) chunk
    double e_trace;                        // |tr(G) - p_eff(n_pr-1)| / (p_eff(n_pr-1))
    double e_center;                       // ||G 1|| / (||G||_F sqrt(n_pr))
    double sym_err;                        // max |G - G^T| BEFORE symmetrizing
    double neg_mass;                       // sum(min(eig,0)) / lam_max
    double min_eig_ratio;                  // min(eig) / lam_max
    long   num_rank;                       // # eig > 1e-8 * lam_max
    int    eigen_threads;                  // effective Eigen thread count (r150 §1)
    std::vector<ChunkRecord> ledger;
};

// Read one R-view dimname vector ("rownames"=features or "colnames"=participants).
inline std::vector<std::string> read_dimname(BigDataStatMeth::hdf5Dataset* d,
                                             const char* which) {
    Rcpp::List dn = d->readDimnames();
    if (!dn.containsElementNamed(which))
        throw std::runtime_error(std::string("dataset has no ") + which);
    Rcpp::CharacterVector cn = dn[which];
    std::vector<std::string> out(cn.size());
    for (R_xlen_t i = 0; i < cn.size(); ++i) out[i] = Rcpp::as<std::string>(cn[i]);
    return out;
}

// K1 primitive. `present_ids` is the FROZEN present-only participant set/order.
inline GramResult present_only_gram(const std::string& file,
                                    const std::string& group,
                                    const std::string& ds,
                                    const std::vector<std::string>& present_ids,
                                    long chunk,
                                    double var_eps = 1e-8) {
    if (chunk < 1) throw std::runtime_error("chunk must be >= 1");
    const long n_pr = (long)present_ids.size();
    if (n_pr < 2) throw std::runtime_error("n_pr must be >= 2 (present participants)");

    // --- open source dataset (R-view p x N: features x all-participants) ---
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(file, group, ds, false));
    d->openDataset();
    if (d->getDatasetptr() == nullptr)
        throw std::runtime_error("cannot open " + group + "/" + ds);
    const long p = (long)d->nrows_r();       // features
    const long N = (long)d->ncols_r();       // all participants in the dataset
    if (p < 1) throw std::runtime_error("dataset has no features");

    // --- read dimnames: rownames=features (for error reporting), colnames=participants ---
    std::vector<std::string> featnm = read_dimname(d.get(), "rownames");
    std::vector<std::string> colnm  = read_dimname(d.get(), "colnames");
    if ((long)featnm.size() != p)
        throw std::runtime_error("rownames length != nrows (feature dimnames mismatch)");
    if ((long)colnm.size() != N)
        throw std::runtime_error("colnames length != ncols (dimnames/dataset mismatch)");
    std::unordered_map<std::string, long> pos;
    pos.reserve(colnm.size() * 2);
    for (long c = 0; c < N; ++c) {
        if (!pos.emplace(colnm[c], c).second)
            throw std::runtime_error("duplicate participant id in dataset colnames: " + colnm[c]);
    }
    std::vector<long> src(n_pr);
    {
        std::unordered_map<std::string, char> seen;
        for (long c = 0; c < n_pr; ++c) {
            auto it = pos.find(present_ids[c]);
            if (it == pos.end())
                throw std::runtime_error("present id not found in dataset: " + present_ids[c]);
            if (!seen.emplace(present_ids[c], 1).second)
                throw std::runtime_error("duplicate present id supplied: " + present_ids[c]);
            src[c] = it->second;
        }
    }

    // --- stream feature-row chunks; standardize per feature; accumulate G ---
    GramResult R;
    R.G = Eigen::MatrixXd::Zero(n_pr, n_pr);
    R.ids = present_ids;
    R.p_eff = p;                 // all features enter (STOP on degenerate)
    R.n_pr = n_pr; R.N_all = N; R.chunk = chunk;
    R.n_chunks = 0; R.last_chunk = 0;
    const double nm1 = (double)(n_pr - 1);

    for (long f0 = 0; f0 < p; f0 += chunk) {
        const long fb = std::min<long>(chunk, p - f0);
        // read fb feature-rows x all N participant-cols (R-view fb x N)
        Eigen::MatrixXd blk(fb, N);
        std::vector<hsize_t> off = {0, (hsize_t)f0}, cnt = {(hsize_t)N, (hsize_t)fb},
                             st = {1, 1}, bl = {1, 1};
        d->readDatasetBlock(off, cnt, st, bl, blk.data());

        // present-only selection -> fb x n_pr
        Eigen::MatrixXd sub(fb, n_pr);
        for (long c = 0; c < n_pr; ++c) sub.col(c) = blk.col(src[c]);

        // finiteness guard (distinct from degeneracy): explicit std::isfinite element
        // scan -- Eigen's allFinite() is not reliable for NaN under this toolchain, so
        // do not depend on it. Report the offending feature (and participant).
        for (long c = 0; c < n_pr; ++c) {
            for (long r = 0; r < fb; ++r) {
                if (!std::isfinite(sub(r, c)))
                    throw std::runtime_error(
                        "non-finite value (NA/NaN/Inf) in present block at feature logical index " +
                        std::to_string(f0 + r) + " (id '" + featnm[f0 + r] +
                        "'), participant '" + present_ids[c] +
                        "' -- run ABORTED, no partial Gram returned");
            }
        }

        // per-feature (per-row) centre + (n_pr - 1) SD standardization
        Eigen::VectorXd rmean = sub.rowwise().mean();
        sub.colwise() -= rmean;                                  // subtract per-row mean
        Eigen::VectorXd ss = sub.rowwise().squaredNorm();        // fb
        Eigen::VectorXd sd = (ss.array() / nm1).sqrt();          // fb
        for (long r = 0; r < fb; ++r) {
            if (!(sd(r) > var_eps))
                throw std::runtime_error(
                    "degenerate feature (sd <= var_eps) at logical index " +
                    std::to_string(f0 + r) + " (id '" + featnm[f0 + r] +
                    "'); upstream QC must remove constant features (oracle script 34 l.45)"
                    " -- run ABORTED, no partial Gram returned");
        }
        Eigen::VectorXd inv = sd.array().inverse();
        sub.array().colwise() *= inv.array();                    // scale each row by 1/sd

        // accumulate the participant Gram
        R.G.noalias() += sub.transpose() * sub;                  // n_pr x n_pr

        // ledger for this chunk (post-standardization diagnostics)
        ChunkRecord rec;
        rec.f_start = f0; rec.f_end = f0 + fb - 1; rec.n_features = fb;
        rec.max_abs_row_sum = sub.rowwise().sum().cwiseAbs().maxCoeff();
        rec.max_abs_ss_dev  = (sub.rowwise().squaredNorm().array() - nm1).abs().maxCoeff();
        rec.trace_contrib   = sub.squaredNorm();
        R.ledger.push_back(rec);
        ++R.n_chunks; R.last_chunk = fb;
    }

    // --- symmetry diagnostic BEFORE mirroring, then symmetrize ---
    R.sym_err = (R.G - R.G.transpose()).cwiseAbs().maxCoeff();
    R.G = 0.5 * (R.G + R.G.transpose());

    // --- completion invariants ---
    const double tr = R.G.trace();
    const double den = (double)R.p_eff * nm1;
    R.e_trace = std::abs(tr - den) / den;
    Eigen::VectorXd g1 = R.G.rowwise().sum();
    const double fro = R.G.norm();
    R.e_center = g1.norm() / (fro * std::sqrt((double)n_pr));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(R.G, Eigen::EigenvaluesOnly);
    Eigen::VectorXd ev = es.eigenvalues();                       // ascending
    const double lam_max = ev.maxCoeff();
    const double scale = (lam_max > 0 ? lam_max : 1.0);
    double negm = 0.0; long rank = 0;
    for (long i = 0; i < ev.size(); ++i) {
        if (ev(i) < 0) negm += ev(i);
        if (ev(i) > 1e-8 * scale) ++rank;
    }
    R.neg_mass = negm / scale;
    R.min_eig_ratio = ev.minCoeff() / scale;
    R.num_rank = rank;
    R.eigen_threads = Eigen::nbThreads();
    return R;
}

}  // namespace reliability
}  // namespace mgcca
#endif  // MGCCA_RELIABILITY_GRAM_HPP
