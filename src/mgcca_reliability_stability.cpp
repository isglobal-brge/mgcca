// mgcca_reliability_stability.cpp -- thin [[Rcpp::export]] wrappers over the K3b subspace-stability
// primitives (mgcca::reliability:: gramc_block_hdf5 / subspace_overlap / proj_overlap). Built ON the
// BigDataStatMeth C++/HDF5 API like the sealed K1/K2/K3a kernels: the Gram-C construction is HDF5-in/HDF5-out
// through BigDataStatMeth (mgcca::read_full / write_full_create); the bounded participant-space dense leaves
// (centring congruence, principal-angle metric) are Eigen, the run_eigen / run_XKX precedent (ChatGPT r232
// approved the scale-based boundary; post-JSS migration of reusable dense helpers into BigDataStatMeth is
// future work). PORT_FOURTH_KERNEL.md. Reproduces script-37's recenter / ovl_rows / trPP.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include "reliabilitySubspaceStability.h"
#include <cmath>
using namespace Rcpp;

static Eigen::MatrixXd asEigen(const Rcpp::NumericMatrix& M) {
    return Eigen::Map<const Eigen::MatrixXd>(M.begin(), M.nrow(), M.ncol());
}

// r232 §B: validate the ORIGINAL R integer BEFORE narrowing to char -- a direct (char) cast can turn
// NA_INTEGER / 256 / 257 / ... into an accepted 0/1, silently bypassing the binary check (fail-open).
static std::vector<char> checked_binary_indicator(const Rcpp::IntegerVector& present) {
    std::vector<char> out(static_cast<std::size_t>(present.size()));
    for (R_xlen_t i = 0; i < present.size(); ++i) {
        const int value = present[i];
        if (value == NA_INTEGER || (value != 0 && value != 1))
            throw std::runtime_error("present must contain only non-missing 0/1 values");
        out[static_cast<std::size_t>(i)] = static_cast<char>(value);
    }
    return out;
}

//' K3b variant-C Gram-C block (HDF5-in/HDF5-out via BigDataStatMeth): H*G[O,O]*H embedded
//' @return A list with the \code{file} and the \code{path} of the centred
//'   Gram-C block written into the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_gramc_block_hdf5(std::string file, std::string in_group, std::string ds_in,
                                        std::string out_group, std::string ds_out,
                                        Rcpp::IntegerVector present, int comp = 0) {
    try {
        H5::Exception::dontPrint();
        std::vector<char> pr = checked_binary_indicator(present);
        mgcca::reliability::gramc_block_hdf5(file, in_group, ds_in, out_group, ds_out, pr, comp);
        return Rcpp::List::create(Rcpp::Named("file") = file,
                                  Rcpp::Named("path") = out_group + "/" + ds_out);
    } catch (H5::Exception& e) { Rf_error("reliability_gramc_block_hdf5 HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) { Rf_error("reliability_gramc_block_hdf5 error: %s", e.what()); }
    return R_NilValue;
}

//' K3b variant-C Gram-C block (in-memory, for unit smoke of the centring math)
//' @return A numeric matrix: the centred Gram-C block \eqn{H G[O,O] H} embedded
//'   back into the full participant indexing.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix reliability_gramc_block_mem(Rcpp::NumericMatrix Gsid, Rcpp::IntegerVector present) {
    try {
        std::vector<char> pr = checked_binary_indicator(present);
        Eigen::MatrixXd C = mgcca::reliability::gramc_block_mem(asEigen(Gsid), pr);
        Rcpp::NumericMatrix M(C.rows(), C.cols());
        if (C.size()) std::copy(C.data(), C.data() + C.size(), M.begin());
        return M;
    } catch (std::exception& e) { Rf_error("reliability_gramc_block_mem error: %s", e.what()); }
    return R_NilValue;
}

//' K3b subspace overlap ovl_rows (frozen rank_tol): mean cos^2 of principal angles; NA if rank-deficient
//' @return A single number: the mean squared cosine of the principal angles
//'   between the two subspaces, or \code{NA} if either basis is rank-deficient.
//' @keywords internal
// [[Rcpp::export]]
double reliability_subspace_overlap(Rcpp::NumericMatrix Va, Rcpp::NumericMatrix Vb) {
    try {
        double v = mgcca::reliability::subspace_overlap(asEigen(Va), asEigen(Vb));   // uses K3B_RANK_TOL
        return std::isnan(v) ? NA_REAL : v;                                          // r232 §E: NaN -> R NA
    } catch (std::exception& e) { Rf_error("reliability_subspace_overlap error: %s", e.what()); }
    return NA_REAL;
}

// NOTE (r234 §1 / r238 §1 / r242 §1): the production R-facing entry exposes NO rank-threshold argument. Boundary
// tests call THIS same production entry with synthetic frames whose rank margins lie on either side of the frozen
// K3B_RANK_TOL. No mutable threshold exists in the wrapper, primitive, or test harness -- K3B_RANK_TOL is the
// single executable source of truth.

//' K3b projector overlap trPP = tr(Pa Pb)/L
//' @return A single number: \eqn{tr(P_a P_b)/L}, the projector overlap of the
//'   two subspaces.
//' @keywords internal
// [[Rcpp::export]]
double reliability_proj_overlap(Rcpp::NumericMatrix Pa, Rcpp::NumericMatrix Pb, int L) {
    try { return mgcca::reliability::proj_overlap(asEigen(Pa), asEigen(Pb), (long)L);
    } catch (std::exception& e) { Rf_error("reliability_proj_overlap error: %s", e.what()); }
    return NA_REAL;
}

//' K3b variant-A small-block Gram via the BigDataStatMeth C++ tcrossprod API: G = Z Z' (participant Gram).
//' Z (in_group/in_ds) = the subset-standardised block table (participants x features, R view). G is written to
//' out_group/out_ds (n x n), isSymmetric -- VERBATIM the sealed estimator run_svd/run_XKX pattern
//' (mgcca_phases.h). BigDataStatMeth::tcrossprod is ADAPTIVE (PATH 1 in-memory preload if A <= ~20% RAM /
//' PATH 2 block-wise streaming otherwise, r244 §1), not unconditionally out-of-core. No R-side BigDataStatMeth
//' algebra dispatch enters the certified small-block Gram path (r242 §2 / r248 §3, PORT_FOURTH_KERNEL §H).
//' Returns file, path, nrow, ncol, block_size, threads_requested (n = participants; for binding + the parity report).
//' @return A list with \code{file}, the \code{path} of the Gram written into the
//'   HDF5 file, its \code{nrow} and \code{ncol}, the \code{block_size} and
//'   \code{threads_requested} (for binding and the parity report).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_tcrossprod_hdf5(std::string file, std::string in_group, std::string in_ds,
                                       std::string out_group, std::string out_ds,
                                       Rcpp::Nullable<int> threads = R_NilValue) {
    try {
        H5::Exception::dontPrint();
        // r252 §1: fail-closed thread contract validated BEFORE any HDF5 handle is opened -- NULL = automatic;
        // if supplied it must be -1 (automatic) or a positive integer (reject NA / 0 / < -1). Deterministic
        // error precedence: an invalid argument causes NO file access. The K3b cert run passes 1L.
        int th = -1;   // -1 == automatic (NULL)
        if (threads.isNotNull()) { th = Rcpp::as<int>(threads);
            if (th == NA_INTEGER || th == 0 || th < -1)
                throw std::runtime_error("reliability_tcrossprod_hdf5: threads must be -1 (automatic) or a positive integer"); }
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dZ(
            new BigDataStatMeth::hdf5Dataset(file, in_group, in_ds, false)); dZ->openDataset();
        if (dZ->getDatasetptr() == nullptr) throw std::runtime_error("reliability_tcrossprod_hdf5: cannot open " + in_group + "/" + in_ds);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dZ2(
            new BigDataStatMeth::hdf5Dataset(file, in_group, in_ds, false)); dZ2->openDataset();
        // overwrite=false => createDataset creates if absent, THROWS if the output already exists (fail-closed,
        // r244 §2.8) -- unlike the estimator's tmp-overwriting pattern; the K3b certified path never overwrites.
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dG(
            new BigDataStatMeth::hdf5Dataset(file, out_group, out_ds, false));
        dG->setCompressionLevel(dZ->getCompressionLevel());
        int blk = BigDataStatMeth::getMaxBlockSize(dZ->nrows(), dZ->ncols(), dZ2->nrows(), dZ2->ncols(), 2, R_NilValue);
        BigDataStatMeth::tcrossprod(dZ.get(), dZ2.get(), dG.get(), true, blk, 0, false, true, threads);  // G = Z Z' (m x m)
        // r250 §1: report R-oriented dims from the OUTPUT dataset (nrows_r/ncols_r), NOT the input's HDF5-native
        // nrows() (which is transposed). For an R input Z (n x p), G is n x n.
        return Rcpp::List::create(Rcpp::Named("file") = file, Rcpp::Named("path") = out_group + "/" + out_ds,
                                  Rcpp::Named("nrow") = (double)dG->nrows_r(), Rcpp::Named("ncol") = (double)dG->ncols_r(),
                                  Rcpp::Named("block_size") = (double)blk, Rcpp::Named("threads_requested") = th);
    } catch (H5::Exception& e) { Rf_error("reliability_tcrossprod_hdf5 HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) { Rf_error("reliability_tcrossprod_hdf5 error: %s", e.what()); }
    return R_NilValue;
}
