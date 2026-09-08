// reliabilitySensitivity.h -- K2 production route: the streamed methylation SENSITIVITY
// accumulators reduced to a PARTICIPANT-SPACE contraction that REUSES the sealed K1 Gram
// (PORT_SECOND_KERNEL.md §14, ChatGPT r168/r170; identity confirmed in tests/harness/31).
//
// Since feature chunks partition the columns of Z, Sum_c Z_c Z_c^T = Z Z^T = Q, so with
//   H_o = alpha_m ( Sm C^T Rr^T + Rr C Sm^T )   (n_fit x n_fit, symmetric)   and D_o = P_m H_o Z,
// the exact script-39 accumulators are
//   accX_o = J   S_pr Q H_o^T P_m ,   accB_o = K_g S_pr Q H_o^T P_m ,
// and Q = E G_K1 E^T where G_K1 is the SEALED K1 present-only Gram, E embeds mids -> fit$ids.
// Minimal calc (r170 §4): T_o = (Q[prm,] H_o) with methyl-absent COLUMNS zeroed; then
//   accX = T_o - colMeans(T_o) ;  accB = cohortRowMean(T_o) - colMeans(T_o).
// Gradient norm (r170 §5): ||D_o||^2 = tr(A_o Q A_o^T), A_o = P_m H_o (raw, before sqrt).
// NO methyl feature stream, NO chunk: all work is n x n. Include AFTER BigDataStatMeth.hpp + mgcca_io.h.
#ifndef MGCCA_RELIABILITY_SENSITIVITY_HPP
#define MGCCA_RELIABILITY_SENSITIVITY_HPP

#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace mgcca {
namespace reliability {

struct SensResult {
    double trbH, trtH, normD2, sym_err;
    Eigen::MatrixXd accB, accX;             // n_pr x n_fit
    long n_fit, n_pr, L, n_meth, n_grp;
    double e_qmm;                           // ||Q[mids,mids] - G_K1|| / ||G_K1||  (== 0 by construction)
};

// K2 production kernel. G_K1 (n_meth x n_meth) is read from HDF5; mids_to_fit[i] is the 0-based
// fit-row index of K1 participant i (K1 column order). present (n_fit), prm (n_pr, 0-based fit
// indices), grp (n_pr, 0-based cohort group id). Rr n_fit x L, Sm n_fit x (n_fit-L), C L x (n_fit-L).
inline SensResult sensitivity_from_gram(const std::string& file, const std::string& group,
                                        const std::string& ds,
                                        const std::vector<long>& mids_to_fit,
                                        const Eigen::MatrixXd& Rr, const Eigen::MatrixXd& Sm,
                                        const Eigen::MatrixXd& C, double alpha,
                                        const std::vector<char>& present,
                                        const std::vector<long>& prm,
                                        const std::vector<int>& grp,
                                        double neg_tol = 1e-8) {
    // --- read the sealed K1 Gram ---
    Eigen::MatrixXd G = mgcca::read_full(file, group, ds);        // n_meth x n_meth
    const long n_meth = (long)G.rows();
    if (G.cols() != n_meth) throw std::runtime_error("K1 Gram is not square");
    if (!G.allFinite()) throw std::runtime_error("K1 Gram has non-finite entries");

    const long n_fit = (long)Rr.rows(), L = (long)Rr.cols();
    if (Sm.rows() != n_fit || Sm.cols() != n_fit - L)
        throw std::runtime_error("Sm dims != n_fit x (n_fit-L)");
    if (C.rows() != L || C.cols() != n_fit - L)
        throw std::runtime_error("C dims != L x (n_fit-L)");
    if ((long)present.size() != n_fit) throw std::runtime_error("present length != n_fit");
    if ((long)mids_to_fit.size() != n_meth) throw std::runtime_error("mids_to_fit length != n_meth");
    if (!Rr.allFinite() || !Sm.allFinite() || !C.allFinite() || !std::isfinite(alpha))
        throw std::runtime_error("non-finite Rr/Sm/C/alpha");
    const long n_pr = (long)prm.size();
    if (n_pr < 2 || (long)grp.size() != n_pr) throw std::runtime_error("prm/grp length mismatch (n_pr>=2)");

    // --- embed Q = E G_K1 E^T (n_fit x n_fit), zero outside mids ---
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n_fit, n_fit);
    for (long i = 0; i < n_meth; ++i) {
        const long fi = mids_to_fit[i];
        if (fi < 0 || fi >= n_fit) throw std::runtime_error("mids_to_fit index out of range");
        for (long j = 0; j < n_meth; ++j) Q(fi, mids_to_fit[j]) = G(i, j);
    }

    // --- H_o via A0 + A0^T (ChatGPT r172 §2: structural symmetry by construction, one product chain,
    //     coefficient applied ONCE) --- A0 = Rr C Sm^T ; H_o = alpha (A0 + A0^T). ---
    Eigen::MatrixXd A0 = Rr * C * Sm.transpose();                 // n_fit x n_fit
    Eigen::MatrixXd Ho = alpha * (A0 + A0.transpose());
    double sym = (Ho - Ho.transpose()).cwiseAbs().maxCoeff();     // ~0 by construction

    // --- T_o = (Q[prm,] H_o) with methyl-absent COLUMNS zeroed (= S_pr Q H_o^T P_m) ---
    Eigen::MatrixXd Qprm(n_pr, n_fit);
    for (long r = 0; r < n_pr; ++r) {
        const long fr = prm[r];
        if (fr < 0 || fr >= n_fit) throw std::runtime_error("prm index out of range");
        Qprm.row(r) = Q.row(fr);
    }
    Eigen::MatrixXd T = Qprm * Ho;                                // n_pr x n_fit  (H_o symmetric)
    for (long c = 0; c < n_fit; ++c) if (!present[c]) T.col(c).setZero();

    // --- accX = T - colMeans ; accB = cohortRowMean(T) - colMeans (r170 §4; no J/K_g matrices) ---
    Eigen::RowVectorXd colmean = T.colwise().mean();             // over the n_pr rows
    long n_grp = 0; for (int g : grp) n_grp = std::max<long>(n_grp, g + 1);
    Eigen::MatrixXd gsum = Eigen::MatrixXd::Zero(n_grp, n_fit);
    std::vector<long> gcnt(n_grp, 0);
    for (long r = 0; r < n_pr; ++r) { gsum.row(grp[r]) += T.row(r); gcnt[grp[r]]++; }
    for (long g = 0; g < n_grp; ++g) { if (gcnt[g] == 0) throw std::runtime_error("empty cohort group");
        gsum.row(g) /= (double)gcnt[g]; }
    SensResult R;
    R.accX = T.rowwise() - colmean;
    R.accB = Eigen::MatrixXd(n_pr, n_fit);
    for (long r = 0; r < n_pr; ++r) R.accB.row(r) = gsum.row(grp[r]) - colmean;

    // --- scalars ---
    R.trbH = R.accB.squaredNorm();
    R.trtH = R.accX.squaredNorm();
    Eigen::MatrixXd A = Ho;                                       // A_o = P_m H_o (zero absent ROWS)
    for (long r = 0; r < n_fit; ++r) if (!present[r]) A.row(r).setZero();
    R.normD2 = (A * Q).cwiseProduct(A).sum();                     // tr(A Q A^T), raw before sqrt
    if (!std::isfinite(R.normD2)) throw std::runtime_error("||D||^2 non-finite");
    const double scale = std::max(1.0, std::abs(R.normD2));
    if (R.normD2 < -neg_tol * scale)
        throw std::runtime_error("||D||^2 materially negative (not FP noise)");

    // Q[mids,mids] == G_K1 by construction (diagnostic)
    double num = 0, den = 0;
    for (long i = 0; i < n_meth; ++i) for (long j = 0; j < n_meth; ++j) {
        const double q = Q(mids_to_fit[i], mids_to_fit[j]);
        num += (q - G(i, j)) * (q - G(i, j)); den += G(i, j) * G(i, j); }
    R.e_qmm = std::sqrt(num) / std::max(std::sqrt(den), 1e-12);

    R.sym_err = sym; R.n_fit = n_fit; R.n_pr = n_pr; R.L = L; R.n_meth = n_meth; R.n_grp = n_grp;
    return R;
}

}  // namespace reliability
}  // namespace mgcca
#endif  // MGCCA_RELIABILITY_SENSITIVITY_HPP
