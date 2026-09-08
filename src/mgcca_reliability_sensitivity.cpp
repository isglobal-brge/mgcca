// mgcca_reliability_sensitivity.cpp -- thin [[Rcpp::export]] wrapper over the K2 production
// kernel mgcca::reliability::sensitivity_from_gram (participant-space contraction reusing the
// sealed K1 Gram). PORT_SECOND_KERNEL.md §14. R supplies the p-independent frozen-fit +
// phenotype constants; C++ reads G_K1 from HDF5 and returns the methyl-only accumulators/scalars.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include "reliabilitySensitivity.h"
using namespace Rcpp;

//' K2: methylation sensitivity accumulators from the sealed K1 Gram (participant-space)
//' @return A list with the sensitivity accumulators \code{trtH} and \code{accX},
//'   the symmetry residual \code{sym_err} and the sizes \code{n_pr},
//'   \code{n_meth}, \code{n_grp}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_sensitivity_gram(std::string file, std::string group, std::string dataset,
                                        Rcpp::IntegerVector mids_to_fit,
                                        Rcpp::NumericMatrix Rr, Rcpp::NumericMatrix Sm,
                                        Rcpp::NumericMatrix C, double alpha,
                                        Rcpp::LogicalVector present,
                                        Rcpp::IntegerVector prm, Rcpp::IntegerVector grp,
                                        double neg_tol = 1e-8) {
    try {
        H5::Exception::dontPrint();
        auto toE = [](const Rcpp::NumericMatrix& M) {
            Eigen::MatrixXd E(M.nrow(), M.ncol());
            std::copy(M.begin(), M.end(), E.data()); return E; };
        std::vector<long> m2f(mids_to_fit.size()); for (R_xlen_t i=0;i<mids_to_fit.size();++i) m2f[i]=mids_to_fit[i];
        std::vector<char> pres(present.size());   for (R_xlen_t i=0;i<present.size();++i)     pres[i]=(present[i]!=0);
        std::vector<long> prmv(prm.size());        for (R_xlen_t i=0;i<prm.size();++i)          prmv[i]=prm[i];
        std::vector<int>  grpv(grp.size());        for (R_xlen_t i=0;i<grp.size();++i)          grpv[i]=grp[i];

        mgcca::reliability::SensResult R = mgcca::reliability::sensitivity_from_gram(
            file, group, dataset, m2f, toE(Rr), toE(Sm), toE(C), alpha, pres, prmv, grpv, neg_tol);

        const int npr = (int)R.n_pr, nf = (int)R.n_fit;
        Rcpp::NumericMatrix accB(npr, nf), accX(npr, nf);
        std::copy(R.accB.data(), R.accB.data() + (std::size_t)npr*nf, accB.begin());
        std::copy(R.accX.data(), R.accX.data() + (std::size_t)npr*nf, accX.begin());
        return Rcpp::List::create(
            Rcpp::Named("trbH") = R.trbH, Rcpp::Named("trtH") = R.trtH,
            Rcpp::Named("normD2") = R.normD2, Rcpp::Named("sym_err") = R.sym_err,
            Rcpp::Named("e_qmm") = R.e_qmm, Rcpp::Named("accB") = accB, Rcpp::Named("accX") = accX,
            Rcpp::Named("n_fit") = (double)R.n_fit, Rcpp::Named("n_pr") = (double)R.n_pr,
            Rcpp::Named("L") = (double)R.L, Rcpp::Named("n_meth") = (double)R.n_meth,
            Rcpp::Named("n_grp") = (double)R.n_grp);
    } catch (H5::Exception& e) {
        Rf_error("reliability_sensitivity_gram HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("reliability_sensitivity_gram error: %s", e.what());
    }
    return R_NilValue;
}
