// mgcca_svd_rcpp — thin [[Rcpp::export]] wrapper over mgcca::run_svd (dual route:
// per-table projector M_j from the thin SVD of X_j, no p x p). See mgcca_phases.hpp.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.hpp"
using namespace Rcpp;

//' MGCCA dual SVD stage: M_j = U diag(w) U' from SVD(X_j)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_svd_rcpp(std::string filename, std::string tmp_group,
                          std::vector<std::string> datasets, int inv_method,
                          Rcpp::Nullable<std::vector<double>> lambda = R_NilValue,
                          std::string svd_method = "full",
                          Rcpp::Nullable<int> threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        std::vector<double> lam;
        if (inv_method == 2) {
            if (lambda.isNull()) throw std::runtime_error("penalized needs lambda");
            lam = Rcpp::as<std::vector<double>>(lambda);
        }
        mgcca::run_svd(filename, tmp_group, datasets, inv_method, lam, svd_method, threads);
        return Rcpp::List::create(
            Rcpp::Named("filename") = filename,
            Rcpp::Named("tmp_group") = tmp_group,
            Rcpp::Named("datasets") = datasets,
            Rcpp::Named("inv_method") = inv_method,
            Rcpp::Named("svd_method") = svd_method);
    } catch (H5::Exception& e) {
        Rf_error("mgcca_svd_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_svd_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
