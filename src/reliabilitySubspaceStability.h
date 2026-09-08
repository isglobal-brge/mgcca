// reliabilitySubspaceStability.h -- K3b: the subspace-STABILITY layer on the sealed K3a engine.
// -----------------------------------------------------------------------------------------------
// The genuinely-new numerical content of K3b (PORT_FOURTH_KERNEL.md, ChatGPT r222-r232). Built ON the
// BigDataStatMeth C++/HDF5 API, exactly like the sealed K1/K2/K3a kernels: BigDataStatMeth owns the scalable
// substrate (HDF5 create/read/write, feature-streamed exact-A methyl Gram = sealed K1, small-block Grams via
// BigDataStatMeth crossprod/tcrossprod, participant-space Gram persistence, and the sealed K3a refit). K3b adds
// only the BOUNDED participant-space dense LEAVES that the estimator port already does in Eigen (the run_eigen /
// run_XKX precedent: SelfAdjointEigenSolver + D*M*D congruence on the small m×m/n×n; a possible post-JSS port of
// reusable dense helpers into BigDataStatMeth is future engineering, NOT a condition of K3b validity -- ChatGPT
// r232 approved the scale-based boundary):
//   (1) gramc_block_hdf5  -- variant-C observed double-centring  H*G[O,O]*H, embedded with zero rows/cols for
//        block-absent participants (== script-37 recenter + embed). Reads the already-reduced participant Gram in
//        FULL from HDF5 (BigDataStatMeth-backed mgcca::read_full), centres it in Eigen, writes it through the
//        BigDataStatMeth-backed I/O (mgcca::write_full_create). The centring itself is a bounded dense
//        participant-space op (Eigen), NOT a block-wise streamed operation.
//   (2) subspace_overlap  -- ovl_rows: mean squared cosine of principal angles between two L-frames (QR-
//        orthonormalise each, SVD of the cross-product), rank guard. Small dense leaf (n×L, L=3) -> Eigen.
//   (3) proj_overlap      -- trPP: tr(Pa Pb)/L for rank-L projectors.
//
// Metric contract transcribed VERBATIM from script 37 (frozen, PORT_FOURTH_KERNEL.md §metric):
//   recenter(G) = H G H,  H = I - 11'/m   (m = nrow of the OBSERVED sub-block)
//   ovl_rows(Va,Vb) = if (qr(Va)$rank<ncol || qr(Vb)$rank<ncol) NA
//                     else mean( svd( t(qr.Q(qr(Va))) %*% qr.Q(qr(Vb)) )$d^2 )
//   trPP(Pa,Pb) = sum(Pa*Pb)/L
// Rank guard (r224 §5.2 / r230 §3): frozen relative threshold K3B_RANK_TOL = 1e-7 (mirrors R's qr() default) via
// ColPivHouseholderQR; rank<ncol -> NaN (the C++ INTERNAL sentinel; the Rcpp wrapper maps it to R NA_real_ --
// NaN is is.na()-true but not literally NA, r232 §E). The R oracle returned finite overlaps for all 22 fixtures,
// establishing full-rank status UNDER THE R CRITERION; the K3b DRIVER records rho(V)=sigma_min/sigma_max for every
// certified frame and MUST demonstrate each lies safely above K3B_RANK_TOL before the C++ real-fixture rank
// decision is certified (r238 §2 -- NOT claimed as completed evidence here). The NaN branch is unit-tested both
// sides of the threshold. Finiteness uses an explicit std::isfinite scan (Eigen::allFinite misses NaN under this
// toolchain -- the K1 lesson).
//
// Include AFTER BigDataStatMeth.hpp + mgcca_io.h.

#ifndef MGCCA_RELIABILITY_SUBSPACE_STABILITY_HPP
#define MGCCA_RELIABILITY_SUBSPACE_STABILITY_HPP

#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace mgcca {
namespace reliability {

// frozen single source of truth for the rank threshold (r232 §D).
inline constexpr double K3B_RANK_TOL = 1e-7;
inline constexpr double K3B_SYM_TOL  = 1e-8;   // relative symmetry criterion for input Grams and projectors

// explicit finiteness scan (Eigen::allFinite misses NaN under this toolchain -- K1 lesson).
inline bool all_finite(const Eigen::MatrixXd& M) {
    const double* p = M.data();
    for (Eigen::Index i = 0; i < M.size(); ++i) if (!std::isfinite(p[i])) return false;
    return true;
}
inline void assert_sym(const Eigen::MatrixXd& G, const char* who) {
    const double sc = std::max(1.0, G.cwiseAbs().maxCoeff());
    if ((G - G.transpose()).cwiseAbs().maxCoeff() > K3B_SYM_TOL * sc)
        throw std::runtime_error(std::string(who) + ": input matrix materially asymmetric");
}

// core observed double-centring + embed (Eigen), shared by the HDF5 and in-memory entry points.
inline Eigen::MatrixXd gramc_core(const Eigen::MatrixXd& Gsid, const std::vector<char>& present, const char* who) {
    const long ns = (long)Gsid.rows();
    if (ns < 2) throw std::runtime_error(std::string(who) + ": Gram must contain at least 2 participants");  // r234 §2: before assert_sym (empty maxCoeff)
    if (Gsid.cols() != ns) throw std::runtime_error(std::string(who) + ": Gsid not square");
    if ((long)present.size() != ns) throw std::runtime_error(std::string(who) + ": present length != nrow(Gsid)");
    if (!all_finite(Gsid)) throw std::runtime_error(std::string(who) + ": non-finite Gram entry");
    assert_sym(Gsid, who);
    std::vector<long> O; O.reserve(ns);
    for (long i = 0; i < ns; ++i) { char v = present[i];
        if (v != 0 && v != 1) throw std::runtime_error(std::string(who) + ": present must be 0/1");
        if (v == 1) O.push_back(i); }
    const long m = (long)O.size();
    if (m < 2) throw std::runtime_error(std::string(who) + ": fewer than 2 observed participants (ineligible block)");
    Eigen::MatrixXd sub(m, m);
    for (long a = 0; a < m; ++a) for (long b = 0; b < m; ++b) sub(a, b) = Gsid(O[a], O[b]);
    Eigen::VectorXd rmean = sub.rowwise().mean();
    Eigen::RowVectorXd cmean = sub.colwise().mean();
    const double gmean = sub.mean();
    Eigen::MatrixXd C = sub; C.colwise() -= rmean; C.rowwise() -= cmean; C.array() += gmean;
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(ns, ns);
    for (long a = 0; a < m; ++a) for (long b = 0; b < m; ++b) out(O[a], O[b]) = C(a, b);
    return out;
}

// (1) variant-C observed double-centring, embedded -- HDF5-in / HDF5-out through the BigDataStatMeth substrate.
inline void gramc_block_hdf5(const std::string& file,
                             const std::string& in_group, const std::string& ds_in,
                             const std::string& out_group, const std::string& ds_out,
                             const std::vector<char>& present, int comp = 0) {
    Eigen::MatrixXd Gsid = mgcca::read_full(file, in_group, ds_in);      // BigDataStatMeth read
    Eigen::MatrixXd out = gramc_core(Gsid, present, "gramc_block_hdf5");
    mgcca::write_full_create(file, out_group, ds_out, out, comp);        // BigDataStatMeth write (fail-closed create)
}

// in-memory variant-C double-centring (matrix in/out) -- for local unit smoke of the centring math.
inline Eigen::MatrixXd gramc_block_mem(const Eigen::MatrixXd& Gsid, const std::vector<char>& present) {
    return gramc_core(Gsid, present, "gramc_block_mem");
}

// thin orthonormal basis of the column space of V (n x k), via Householder QR (economy Q).
inline Eigen::MatrixXd thin_Q(const Eigen::MatrixXd& V) {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(V);
    return qr.householderQ() * Eigen::MatrixXd::Identity(V.rows(), V.cols());
}
// numeric column rank of V at the FROZEN K3B_RANK_TOL (r240 §1: no mutable threshold anywhere -- the constant is
// the single executable source of truth; the removed test-only wrapper was the only reason to parameterise it).
inline long qr_rank(const Eigen::MatrixXd& V) {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(V); qr.setThreshold(K3B_RANK_TOL); return (long)qr.rank();
}

// (2) ovl_rows: mean squared cosine of the principal angles between the column spaces of Va, Vb.
inline double subspace_overlap(const Eigen::MatrixXd& Va, const Eigen::MatrixXd& Vb) {
    if (Va.rows() != Vb.rows()) throw std::runtime_error("subspace_overlap: row count differs");
    if (Va.cols() != Vb.cols()) throw std::runtime_error("subspace_overlap: col count differs");
    if (Va.cols() < 1) throw std::runtime_error("subspace_overlap: need >= 1 column");
    if (!all_finite(Va) || !all_finite(Vb)) throw std::runtime_error("subspace_overlap: non-finite frame entry");
    if (qr_rank(Va) < Va.cols() || qr_rank(Vb) < Vb.cols())
        return std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd Qa = thin_Q(Va), Qb = thin_Q(Vb);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Qa.transpose() * Qb);
    return svd.singularValues().array().square().mean();
}

// (3) trPP: tr(Pa Pb)/L for rank-L projectors.
inline double proj_overlap(const Eigen::MatrixXd& Pa, const Eigen::MatrixXd& Pb, long L) {
    if (Pa.rows() != Pa.cols() || Pb.rows() != Pb.cols()) throw std::runtime_error("proj_overlap: projectors must be square");
    if (Pa.rows() != Pb.rows() || Pa.cols() != Pb.cols()) throw std::runtime_error("proj_overlap: dim mismatch");
    if (L < 1 || L > Pa.rows()) throw std::runtime_error("proj_overlap: need 1 <= L <= nrow(P)");
    if (!all_finite(Pa) || !all_finite(Pb)) throw std::runtime_error("proj_overlap: non-finite projector entry");
    // r234 §3: self-defend as the other primitives do -- projectors must be symmetric (the sealed K3a output P is
    // symmetric-idempotent; the driver additionally binds provenance/order/idempotence). trPP would silently
    // return a meaningless number on a non-projector, so guard rather than trust the caller.
    assert_sym(Pa, "proj_overlap(Pa)"); assert_sym(Pb, "proj_overlap(Pb)");
    return Pa.cwiseProduct(Pb).sum() / (double)L;
}

} // namespace reliability
} // namespace mgcca
#endif
