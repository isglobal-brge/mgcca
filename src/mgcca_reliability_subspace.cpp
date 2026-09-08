// mgcca_reliability_subspace.cpp -- thin [[Rcpp::export]] wrapper over the K3a production
// kernel mgcca::reliability::fit_from_grams (top-L balanced subspace refit from block Grams).
// PORT_THIRD_KERNEL.md (sealed ChatGPT r192). R supplies the J Gram HDF5 dataset names + the
// per-block presence masks + frozen params (L, gamma); C++ reads the Grams and returns the FULL
// fit object (r190 §2 / r192): D, lambda, B_j, unbalanced (P0,V0,mu0), balanced (P,V,mu), w, q,
// cL, gaps, structural invariants, and the solve audit (rho_j, r^max_j, sym(B_j)).
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include "reliabilitySubspace.h"
using namespace Rcpp;

//' K3a: top-L balanced subspace refit from block Grams (participant-space, dense Eigen core)
//' @return A list with the refitted basis \code{V} and projector \code{P}, the
//'   retained rank \code{L}, the balancing weights \code{lambda}, the spectrum
//'   \code{mu} and its \code{gap}, the \code{eligible} participants and the
//'   numerical validity flags (\code{solve_ok}, \code{symSbal}, \code{P_idem},
//'   \code{eig_min}, \code{n_clipped}, ...).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_subspace_fit(std::string file, std::string group,
                                    Rcpp::CharacterVector gram_ds,
                                    Rcpp::IntegerMatrix masks,   // n x J presence (1/0), column j = block j
                                    int L, double gamma,
                                    double eps = 1e-8, double rank_tol = 1e-8,
                                    double sym_tol = 1e-8, double pad_tol = 1e-8,
                                    bool return_B = true) {   // r198 §7: canonical full-fit includes B_j
    try {
        H5::Exception::dontPrint();
        const long J = gram_ds.size(), n = masks.nrow();
        if (masks.ncol() != J) throw std::runtime_error("masks must be n x J (one column per block)");
        std::vector<std::string> ds(J); for (long j=0;j<J;++j) ds[j] = Rcpp::as<std::string>(gram_ds[j]);
        std::vector<std::vector<char> > mk(J, std::vector<char>(n));
        for (long j=0;j<J;++j) for (long i=0;i<n;++i) {
            int v = masks(i,j);                                   // r196 §5: masks MUST be strictly binary
            if (v != 0 && v != 1) throw std::runtime_error("masks must be binary (0/1)");
            mk[j][i] = (char)v;
        }

        mgcca::reliability::SubspaceFit R = mgcca::reliability::fit_from_grams(
            file, group, ds, mk, (long)L, gamma, eps, rank_tol, sym_tol, pad_tol);

        auto toM = [](const Eigen::MatrixXd& E){ Rcpp::NumericMatrix M(E.rows(), E.cols());
            if (E.size()) std::copy(E.data(), E.data()+E.size(), M.begin()); return M; };
        auto toV = [](const Eigen::VectorXd& v){ return Rcpp::NumericVector(v.data(), v.data()+v.size()); };
        auto toIV = [](const Eigen::VectorXi& v){ return Rcpp::IntegerVector(v.data(), v.data()+v.size()); };
        Rcpp::LogicalVector sok(R.solve_ok.size());
        for (R_xlen_t i=0;i<(R_xlen_t)R.solve_ok.size();++i) sok[i] = (R.solve_ok[i] != 0);

        Rcpp::List out = Rcpp::List::create(
            Rcpp::Named("eligible") = R.eligible,
            Rcpp::Named("n") = (double)R.n, Rcpp::Named("J") = (double)R.J, Rcpp::Named("L") = (double)R.L,
            Rcpp::Named("D") = toV(R.D), Rcpp::Named("lambda") = toV(R.lambda),
            Rcpp::Named("rank_num") = toIV(R.rank_num), Rcpp::Named("cL") = toV(R.cL),
            Rcpp::Named("w") = toV(R.w), Rcpp::Named("q") = toV(R.q),
            Rcpp::Named("mu0") = toV(R.mu0), Rcpp::Named("mu") = toV(R.mu),
            Rcpp::Named("V0") = toM(R.V0), Rcpp::Named("V") = toM(R.V),
            Rcpp::Named("P0") = toM(R.P0), Rcpp::Named("P") = toM(R.P),
            Rcpp::Named("gap0") = R.gap0, Rcpp::Named("gap") = R.gap);
        Rcpp::List diag = Rcpp::List::create(
            Rcpp::Named("solve_rho") = toV(R.solve_rho), Rcpp::Named("solve_rmax") = toV(R.solve_rmax),
            Rcpp::Named("symB") = toV(R.symB), Rcpp::Named("solve_ok") = sok,
            Rcpp::Named("symS0") = R.symS0, Rcpp::Named("symSbal") = R.symSbal,
            Rcpp::Named("P_sym") = R.P_sym, Rcpp::Named("P_idem") = R.P_idem,
            Rcpp::Named("trP") = R.trP, Rcpp::Named("sumq") = R.sumq, Rcpp::Named("meanw") = R.meanw,
            Rcpp::Named("symG") = toV(R.symG), Rcpp::Named("eig_min") = toV(R.eig_min),
            Rcpp::Named("eig_max") = toV(R.eig_max), Rcpp::Named("lam_sum_raw") = toV(R.lam_sum_raw),
            Rcpp::Named("lam_sum_clip") = toV(R.lam_sum_clip), Rcpp::Named("n_clipped") = toIV(R.n_clipped));
        out["diagnostics"] = diag;
        if (return_B) {
            Rcpp::List Bl(R.B.size());
            for (std::size_t j=0;j<R.B.size();++j) Bl[j] = toM(R.B[j]);
            out["B"] = Bl;
        }
        return out;
    } catch (H5::Exception& e) {
        Rf_error("reliability_subspace_fit HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("reliability_subspace_fit error: %s", e.what());
    }
    return R_NilValue;
}

//' K3a orientation check: read an HDF5 dataset via the EXACT helper the kernel uses (mgcca::read_full)
//' so an asymmetric fixture can prove no hidden transpose (r198 §4).
//' @return A numeric matrix: the HDF5 dataset read back in the R view through the
//'   kernel's own reader.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix reliability_read_full_test(std::string file, std::string group, std::string dataset) {
    try {
        H5::Exception::dontPrint();
        Eigen::MatrixXd M = mgcca::read_full(file, group, dataset);
        Rcpp::NumericMatrix out(M.rows(), M.cols());
        if (M.size()) std::copy(M.data(), M.data() + M.size(), out.begin());
        return out;
    } catch (H5::Exception& e) {
        Rf_error("reliability_read_full_test HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("reliability_read_full_test error: %s", e.what());
    }
    return R_NilValue;
}
