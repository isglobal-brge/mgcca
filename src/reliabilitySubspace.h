// reliabilitySubspace.h -- K3a production kernel: the top-L (Ky Fan) balanced subspace
// REFIT from block Grams (PORT_THIRD_KERNEL.md, sealed ChatGPT r192). Faithful transcription
// of the R oracle `assemble5` = fit_pgcca_from_grams + topL_weights + rebalance + per-block
// geometry (reliability/helix/institution_balancing.R + lam_from_grams in script 37).
//
// HDF5 is the input contract (Grams read with mgcca::read_full); the n×n refit runs on a
// DENSE Eigen core (SelfAdjointEigenSolver, as run_eigen/run_XKX already do for small m×m).
//
// Exact op order (r190 E): B_j = G_j R_j (NOT symmetrised per term) -> sum -> diagonal-scale
// D·M·D -> symmetrise -> eigensolve; each A_j symmetrised before its own eigensolve.
// Distinct projectors (r190 A): unbalanced P0 (S0 = D(Σ B_j)D), balanced P (S_bal = D(Σ w_j B_j)D).
// D is a DIAGONAL operator (r190 B): d_i = K_i^{-1/2}, K_i = Σ_j m_{ij}.
// Solve audit (r192): C_j = G_j+λ_j I ; ρ_j = ‖C_jR_j−I‖_F/(‖C_j‖_F‖R_j‖_F+‖I‖_F) ; r^max_j.
// Include AFTER BigDataStatMeth.hpp + mgcca_io.h.
#ifndef MGCCA_RELIABILITY_SUBSPACE_HPP
#define MGCCA_RELIABILITY_SUBSPACE_HPP

#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mgcca {
namespace reliability {

struct SubspaceFit {
    bool eligible = false;
    long n = 0, J = 0, L = 0;
    Eigen::VectorXd D;                     // n         (d_i = K_i^{-1/2})
    Eigen::VectorXd lambda;                // J         (regularisation, RAW-G_j spectrum)
    Eigen::VectorXi rank_num;              // J         (numeric rank of A_j)
    Eigen::VectorXd cL;                    // J         (Ky Fan top-L capacity)
    Eigen::VectorXd w;                     // J         (mean-one weights; empty if ineligible)
    Eigen::VectorXd q;                     // J         (contribution shares; empty if ineligible)
    Eigen::VectorXd mu0, mu;               // n         (unbalanced / balanced spectra, descending)
    Eigen::MatrixXd V0, V;                 // n×L       (leading eigenvectors)
    Eigen::MatrixXd P0, P;                 // n×n       (top-L projectors)
    double gap0 = 0, gap = 0;              // μ_L − μ_{L+1}
    // solve audit (r192) + block diagnostics
    Eigen::VectorXd solve_rho, solve_rmax, symB;   // J
    std::vector<char> solve_ok;                    // J
    // input-Gram + lambda diagnostics (r196 §2/§6)
    Eigen::VectorXd symG;                          // J  input Gram asymmetry max|G-Gᵀ|
    Eigen::VectorXd eig_min, eig_max;              // J  raw G_j spectrum extremes
    Eigen::VectorXd lam_sum_raw, lam_sum_clip;     // J  Ση vs Σpmax(η,0)  (difference = clip effect)
    Eigen::VectorXi n_clipped;                     // J  # raw eigenvalues < 0 truncated by pmax
    // structural invariants (r190 F)
    double symS0 = 0, symSbal = 0, P_sym = 0, P_idem = 0, trP = 0;
    double sumq = 0, meanw = 0;
    std::vector<Eigen::MatrixXd> B;        // J × (n×n)  (persisted reusable estimator quantity)
};

// descending symmetric eigendecomposition: values[0] >= values[1] >= ...
inline void eig_desc(const Eigen::MatrixXd& Sin, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs) {
    Eigen::MatrixXd S = 0.5 * (Sin + Sin.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    if (es.info() != Eigen::Success) throw std::runtime_error("subspace: eigensolver failed");
    vals = es.eigenvalues().reverse();                       // ascending -> descending
    vecs = es.eigenvectors().rowwise().reverse();            // columns reversed to match
}

// K3a refit engine. gram_ds[j] = HDF5 dataset name (in `group`) of the n×n block Gram G_j
// (padded absent rows/cols = 0). masks[j] = length-n presence (1/0). L, gamma frozen params.
inline SubspaceFit fit_from_grams(const std::string& file, const std::string& group,
                                  const std::vector<std::string>& gram_ds,
                                  const std::vector<std::vector<char> >& masks,
                                  long L, double gamma,
                                  double eps = 1e-8, double rank_tol = 1e-8,
                                  double sym_tol = 1e-8, double pad_tol = 1e-8) {
    SubspaceFit R;
    const long J = (long)gram_ds.size();
    if (J < 1) throw std::runtime_error("subspace: need >= 1 block");
    if ((long)masks.size() != J) throw std::runtime_error("subspace: masks length != J");
    if (L < 1) throw std::runtime_error("subspace: L < 1");
    // r198 hardening: validate numerical parameters up front (NaN tolerances would disable guards, since
    // every comparison against NaN is false).
    if (!(std::isfinite(gamma) && gamma > 0)) throw std::runtime_error("subspace: gamma must be finite > 0");
    if (!(std::isfinite(eps) && eps > 0)) throw std::runtime_error("subspace: eps must be finite > 0");
    if (!(std::isfinite(rank_tol) && rank_tol > 0 && rank_tol < 1)) throw std::runtime_error("subspace: rank_tol must be finite in (0,1)");
    if (!(std::isfinite(sym_tol) && sym_tol >= 0)) throw std::runtime_error("subspace: sym_tol must be finite >= 0");
    if (!(std::isfinite(pad_tol) && pad_tol >= 0)) throw std::runtime_error("subspace: pad_tol must be finite >= 0");

    // --- read Grams; validate square / finite / equal dims (r190 D block-contract) ---
    std::vector<Eigen::MatrixXd> G(J);
    long n = -1;
    R.symG = Eigen::VectorXd(J);
    for (long j = 0; j < J; ++j) {
        G[j] = mgcca::read_full(file, group, gram_ds[j]);
        if (G[j].rows() < 1 || G[j].cols() < 1) throw std::runtime_error("subspace: empty Gram");  // r198: before maxCoeff()
        if (G[j].rows() != G[j].cols()) throw std::runtime_error("subspace: Gram not square");
        if (n < 0) n = (long)G[j].rows();
        else if ((long)G[j].rows() != n) throw std::runtime_error("subspace: Gram dims differ across blocks");
        if (!G[j].allFinite()) throw std::runtime_error("subspace: non-finite Gram entry");
        if ((long)masks[j].size() != n) throw std::runtime_error("subspace: mask length != n");
        const double gscale = std::max(1.0, G[j].cwiseAbs().maxCoeff());
        // r196 §6: input Grams MUST be symmetric within tol -- a materially asymmetric Gram is an INPUT ERROR,
        // NOT silently symmetrised. (The estimator only ever forms symmetric block Grams; the eigensolver below
        // receives sym(G_j) so its triangle choice is irrelevant, but a bad input must fail closed here.)
        R.symG(j) = (G[j] - G[j].transpose()).cwiseAbs().maxCoeff();
        if (R.symG(j) > sym_tol * gscale)
            throw std::runtime_error("subspace: input Gram materially asymmetric (block " + std::to_string((long long)j) + ")");
        // r196 §5: absent (mask==0) participants MUST be zero-padded rows AND columns (fail closed otherwise).
        for (long i = 0; i < n; ++i) if (!masks[j][i]) {
            if (G[j].row(i).cwiseAbs().maxCoeff() > pad_tol * gscale ||
                G[j].col(i).cwiseAbs().maxCoeff() > pad_tol * gscale)
                throw std::runtime_error("subspace: absent participant has non-zero Gram row/col (bad zero-padding)");
        }
    }
    if (n < 1 || L > n) throw std::runtime_error("subspace: bad n / L>n");
    R.n = n; R.J = J; R.L = L;

    // --- K_i = Σ_j m_{ij} ; D = 1/sqrt(K) (fail-closed if any K<1, r190 B/D) ---
    Eigen::VectorXd D(n);
    for (long i = 0; i < n; ++i) {
        long Ki = 0; for (long j = 0; j < J; ++j) Ki += masks[j][i] ? 1 : 0;
        if (Ki < 1) throw std::runtime_error("subspace: participant present in ZERO blocks (K<1)");
        D(i) = 1.0 / std::sqrt((double)Ki);
    }
    if (!D.allFinite()) throw std::runtime_error("subspace: non-finite D");
    R.D = D;

    // --- per block: λ_j (RAW G_j spectrum) ; R_j = (G_j+λ_jI)^{-1} (SPD) ; B_j = G_j R_j ; solve audit ---
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    R.lambda = Eigen::VectorXd(J);
    R.B.assign(J, Eigen::MatrixXd());
    R.solve_rho = Eigen::VectorXd(J); R.solve_rmax = Eigen::VectorXd(J);
    R.symB = Eigen::VectorXd(J); R.solve_ok.assign(J, 0);
    R.eig_min = Eigen::VectorXd(J); R.eig_max = Eigen::VectorXd(J);
    R.lam_sum_raw = Eigen::VectorXd(J); R.lam_sum_clip = Eigen::VectorXd(J); R.n_clipped = Eigen::VectorXi(J);
    for (long j = 0; j < J; ++j) {
        // λ_j = γ · Σ_r ẽ_r / #{ẽ_r > rank_tol·max ẽ}  with ẽ = pmax(η,0) (RAW G_j spectrum, negatives
        // truncated to 0 BEFORE sum/rank -- verbatim lam_from_grams in script 37; topL below does NOT truncate).
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eg(0.5 * (G[j] + G[j].transpose()),
                                                          Eigen::EigenvaluesOnly);
        if (eg.info() != Eigen::Success) throw std::runtime_error("subspace: G_j eigenvalues failed");
        Eigen::VectorXd eta = eg.eigenvalues();
        double mx = 0, s = 0; long nclip = 0;
        for (long r = 0; r < eta.size(); ++r) { double e = eta(r) > 0 ? eta(r) : 0.0; if (e > mx) mx = e; s += e; if (eta(r) < 0) nclip++; }
        if (!std::isfinite(mx) || mx <= 0) throw std::runtime_error("subspace: G_j max eigenvalue non-positive/non-finite");
        long rk = 0;
        for (long r = 0; r < eta.size(); ++r) { double e = eta(r) > 0 ? eta(r) : 0.0; if (e > rank_tol * mx) rk++; }
        if (rk < 1) throw std::runtime_error("subspace: G_j numeric rank 0");
        double lam = gamma * s / (double)rk;
        R.eig_min(j) = eta.minCoeff(); R.eig_max(j) = eta.maxCoeff();          // r196 §2 λ audit
        R.lam_sum_raw(j) = eta.sum(); R.lam_sum_clip(j) = s; R.n_clipped(j) = (int)nclip;
        if (!std::isfinite(lam) || lam <= 0) throw std::runtime_error("subspace: λ_j non-positive/non-finite");
        R.lambda(j) = lam;
        // R_j via SPD cholesky solve of C_j = G_j + λ_j I ; B_j = G_j R_j
        Eigen::MatrixXd Cj = G[j] + lam * I;
        Eigen::LLT<Eigen::MatrixXd> llt(0.5 * (Cj + Cj.transpose()));
        R.solve_ok[j] = (llt.info() == Eigen::Success) ? 1 : 0;
        if (!R.solve_ok[j]) throw std::runtime_error("subspace: SPD factorisation of (G_j+λ_jI) failed");
        Eigen::MatrixXd Rj = llt.solve(I);
        if (!Rj.allFinite()) throw std::runtime_error("subspace: R_j non-finite");           // r198 finiteness
        R.B[j] = G[j] * Rj;
        if (!R.B[j].allFinite()) throw std::runtime_error("subspace: B_j non-finite");
        // solve audit (r192): ρ_j, r^max_j, symmetry residual of B_j
        Eigen::MatrixXd resid = Cj * Rj - I;
        double denom = Cj.norm() * Rj.norm() + std::sqrt((double)n);   // ‖I‖_F = sqrt(n)
        R.solve_rho(j)  = resid.norm() / std::max(denom, 1e-300);
        R.solve_rmax(j) = resid.cwiseAbs().maxCoeff();
        R.symB(j)       = (R.B[j] - R.B[j].transpose()).cwiseAbs().maxCoeff();
        if (!std::isfinite(R.solve_rho(j)) || !std::isfinite(R.solve_rmax(j)) || !std::isfinite(R.symB(j)))
            throw std::runtime_error("subspace: non-finite solve diagnostic");
    }

    // --- unbalanced fit -> P0 (r190 A): M0 = Σ B_j ; S0 = D M0 D ; sym ; eigen ; P0 ---
    Eigen::MatrixXd M0 = Eigen::MatrixXd::Zero(n, n);
    for (long j = 0; j < J; ++j) M0 += R.B[j];
    Eigen::MatrixXd S0 = D.asDiagonal() * M0 * D.asDiagonal();
    R.symS0 = (S0 - S0.transpose()).cwiseAbs().maxCoeff();
    eig_desc(S0, R.mu0, R.V0);
    if (!R.mu0.allFinite() || !R.V0.allFinite()) throw std::runtime_error("subspace: non-finite unbalanced spectrum");
    Eigen::MatrixXd V0L = R.V0.leftCols(L);
    R.P0 = V0L * V0L.transpose();
    R.V0 = V0L;
    R.gap0 = R.mu0(L - 1) - (L < n ? R.mu0(L) : 0.0);

    // --- top-L Ky Fan weights: A_j = D B_j D (sym) ; eigen ; c_{j,L} ; eligibility ; w_j ---
    std::vector<Eigen::MatrixXd> A(J);
    R.cL = Eigen::VectorXd(J); R.rank_num = Eigen::VectorXi(J);
    bool eligible = true;
    for (long j = 0; j < J; ++j) {
        A[j] = D.asDiagonal() * R.B[j] * D.asDiagonal();
        A[j] = 0.5 * (A[j] + A[j].transpose());
        Eigen::VectorXd ev; Eigen::MatrixXd evec;
        eig_desc(A[j], ev, evec);                          // descending
        double mx = ev(0);
        long rk = 0; for (long r = 0; r < ev.size(); ++r) if (ev(r) > rank_tol * mx) rk++;
        R.rank_num(j) = (int)rk;
        double c = 0; for (long l = 0; l < L; ++l) c += ev(l);
        R.cL(j) = c;
        if (rk < L || !(std::isfinite(c) && c > eps)) eligible = false;  // r190 status change; non-finite/<=eps -> ineligible
    }
    R.eligible = eligible;
    if (!eligible) return R;                                // caller decides (mirrors assemble5 eligible=FALSE)

    // w_j = J·(1/c_{j,L}) / Σ_k (1/c_{k,L})   (mean-one)
    Eigen::VectorXd inv = R.cL.cwiseInverse();
    R.w = (double)J * inv / inv.sum();
    if (!R.w.allFinite() || (R.w.array() <= 0.0).any()) throw std::runtime_error("subspace: non-finite/non-positive weights");
    R.meanw = R.w.mean();

    // --- balanced refit -> final P (r190 A): S_bal = D(Σ w_j B_j)D ; sym ; eigen ; P ; gap ---
    Eigen::MatrixXd Mw = Eigen::MatrixXd::Zero(n, n);
    for (long j = 0; j < J; ++j) Mw += R.w(j) * R.B[j];
    Eigen::MatrixXd Sbal = D.asDiagonal() * Mw * D.asDiagonal();
    R.symSbal = (Sbal - Sbal.transpose()).cwiseAbs().maxCoeff();
    eig_desc(Sbal, R.mu, R.V);
    if (!R.mu.allFinite() || !R.V.allFinite()) throw std::runtime_error("subspace: non-finite balanced spectrum");
    Eigen::MatrixXd VL = R.V.leftCols(L);
    R.P = VL * VL.transpose();
    R.V = VL;
    R.gap = R.mu(L - 1) - (L < n ? R.mu(L) : 0.0);

    // --- contribution shares (r190 C): q_j = tr(P·w_j A_j) / tr(P·S_bal) ---
    Eigen::VectorXd num(J);
    for (long j = 0; j < J; ++j) num(j) = R.w(j) * R.P.cwiseProduct(A[j]).sum();   // A_j symmetric
    double den = R.P.cwiseProduct(Sbal).sum();                                     // tr(P S_bal)
    if (!std::isfinite(den) || std::abs(den) < 1e-300) throw std::runtime_error("subspace: tr(P S_bal) ~ 0");
    R.q = num / den;
    if (!R.q.allFinite()) throw std::runtime_error("subspace: non-finite q");
    R.sumq = R.q.sum();

    // --- structural invariants (r190 F) ---
    R.P_sym  = (R.P - R.P.transpose()).cwiseAbs().maxCoeff();
    R.P_idem = (R.P * R.P - R.P).cwiseAbs().maxCoeff();
    R.trP    = R.P.trace();
    return R;
}

}  // namespace reliability
}  // namespace mgcca
#endif  // MGCCA_RELIABILITY_SUBSPACE_HPP
