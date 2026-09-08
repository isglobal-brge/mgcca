// Rcpp entry points for the reliability layer's participant-space phases.
// Implementation in src/reliabilityPhases.h. Family B (path-based) like
// mgcca_rcpp: receives paths, opens and creates through BigDataStatMeth, returns
// results. Errors follow the house rule -- throw std::runtime_error inside the
// try, Rf_error in the catch, never Rcpp::stop.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "reliabilityPhases.h"
using namespace Rcpp;

static Eigen::MatrixXd asEigenM(const Rcpp::NumericMatrix& M) {
    Eigen::MatrixXd E(M.nrow(), M.ncol());
    std::copy(M.begin(), M.end(), E.data());
    return E;
}
static Rcpp::NumericMatrix asRcppM(const Eigen::MatrixXd& E) {
    Rcpp::NumericMatrix M(E.rows(), E.cols());
    std::copy(E.data(), E.data() + E.size(), M.begin());
    return M;
}

//' Present-only participant Gram of one block, out-of-core (reliability phase R1)
//'
//' @param block_size 0 lets BigDataStatMeth choose; a positive value FORCES that
//'   block size, so the blocked path can be exercised on small fixtures.
//' @return A list with \code{filename}, the \code{path} of the Gram written into
//'   the HDF5 file, the participant \code{ids}, the \code{block_size} actually
//'   used and whether the block was \code{standardized}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_gram_block(std::string filename, std::string in_group,
                                  std::string dataset, std::string out_group,
                                  int block_size = 0,
                                  Rcpp::Nullable<int> threads = R_NilValue,
                                  bool standardize = true) {
    try {
        H5::Exception::dontPrint();
        Rcpp::CharacterVector ids = mgcca::reliability::gram_block(
            filename, in_group, dataset, out_group, block_size, threads, standardize);
        return Rcpp::List::create(
            Rcpp::Named("filename") = filename,
            Rcpp::Named("path")     = out_group + "/" + dataset,
            Rcpp::Named("ids")      = ids,
            Rcpp::Named("block_size") = block_size,
            Rcpp::Named("standardized") = standardize);
    } catch (H5::Exception& e) {
        Rf_error("reliability_gram_block HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("reliability_gram_block error: %s", e.what());
    }
    return R_NilValue;
}

//' Weighted operator, eigenbasis and per-block resolvents (reliability phase R2)
//' @return A list with the retained eigenbasis \code{V}, the eigenvalues
//'   \code{D}, the spectral \code{gap} at the retained rank and that rank
//'   \code{L}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_reference_fit(Rcpp::List Glist, Rcpp::List present,
                                     Rcpp::NumericVector lambda, int L) {
    try {
        H5::Exception::dontPrint();
        const std::size_t J = (std::size_t)Glist.size();
        if ((std::size_t)present.size() != J || (std::size_t)lambda.size() != J)
            Rf_error("Glist, present and lambda must have the same length");
        std::vector<Eigen::MatrixXd> G(J);
        std::vector<std::vector<char> > P(J);
        std::vector<double> lam(J);
        for (std::size_t j = 0; j < J; ++j) {
            G[j] = asEigenM(Rcpp::as<Rcpp::NumericMatrix>(Glist[j]));
            Rcpp::LogicalVector pv = Rcpp::as<Rcpp::LogicalVector>(present[j]);
            P[j].resize((std::size_t)pv.size());
            for (R_xlen_t i = 0; i < pv.size(); ++i) P[j][(std::size_t)i] = (pv[i] != 0);
            lam[j] = lambda[j];
        }
        mgcca::reliability::RefFit F =
            mgcca::reliability::reference_fit(G, P, lam, (long)L);

        Rcpp::List Rl(J);
        for (std::size_t j = 0; j < J; ++j) Rl[j] = asRcppM(F.Rlist[j]);
        Rcpp::NumericVector mu(F.mu.size()), D(F.D.size());
        std::copy(F.mu.data(), F.mu.data() + F.mu.size(), mu.begin());
        std::copy(F.D.data(),  F.D.data()  + F.D.size(),  D.begin());
        return Rcpp::List::create(
            Rcpp::Named("S") = asRcppM(F.S), Rcpp::Named("V") = asRcppM(F.V),
            Rcpp::Named("mu") = mu, Rcpp::Named("D") = D,
            Rcpp::Named("Rlist") = Rl, Rcpp::Named("gap") = F.gap,
            Rcpp::Named("L") = L);
    } catch (std::exception& e) {
        Rf_error("reliability_reference_fit error: %s", e.what());
    }
    return R_NilValue;
}

//' Per-block inputs for the sealed sensitivity kernel (phase R3b)
//'
//' Rr_j = R_j D V_L and Sm_j = R_j D V_R. Kept in C++ so the algebra lives in one
//' place; R only dispatches.
//' @return A list with the two per-block input matrices \code{Rr}
//'   (\eqn{R_j D V_L}) and \code{Sm} (\eqn{R_j D V_R}).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_block_inputs(Rcpp::NumericMatrix Rj, Rcpp::NumericVector D,
                                    Rcpp::NumericMatrix V, int L) {
    try {
        Eigen::MatrixXd R = asEigenM(Rj), Vm = asEigenM(V);
        Eigen::VectorXd d(D.size());
        std::copy(D.begin(), D.end(), d.data());
        const long n = Vm.rows();
        if (L < 1 || L >= n) Rf_error("L must satisfy 1 <= L < n");
        Eigen::MatrixXd DV = d.asDiagonal() * Vm;
        return Rcpp::List::create(
            Rcpp::Named("Rr") = asRcppM(R * DV.leftCols(L)),
            Rcpp::Named("Sm") = asRcppM(R * DV.rightCols(n - L)));
    } catch (std::exception& e) {
        Rf_error("reliability_block_inputs error: %s", e.what());
    }
    return R_NilValue;
}

//' Alignment functional and first-order perturbation coefficients (phase R3)
//' @return A list with the alignment functional \code{T}, the first-order
//'   perturbation coefficients \code{C} and the per-component weights \code{w}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_query(Rcpp::NumericMatrix V, Rcpp::NumericVector mu,
                             Rcpp::NumericVector zstar, int L) {
    try {
        mgcca::reliability::RefFit F;
        F.V = asEigenM(V);
        F.mu = Eigen::VectorXd(mu.size());
        std::copy(mu.begin(), mu.end(), F.mu.data());
        Eigen::VectorXd zs(zstar.size());
        std::copy(zstar.begin(), zstar.end(), zs.data());
        // w = V_L' z*, the query's coordinates ON the retained axes. Returned so the
        // map view can draw where a query points without recomputing anything.
        Eigen::VectorXd w = F.V.leftCols(L).transpose() * zs;
        Rcpp::NumericVector wr(w.size());
        std::copy(w.data(), w.data() + w.size(), wr.begin());
        return Rcpp::List::create(
            Rcpp::Named("T") = mgcca::reliability::alignment_T(F, zs, (long)L),
            Rcpp::Named("C") = asRcppM(mgcca::reliability::perturbation_C(F, zs, (long)L)),
            Rcpp::Named("w") = wr);
    } catch (std::exception& e) {
        Rf_error("reliability_query error: %s", e.what());
    }
    return R_NilValue;
}
