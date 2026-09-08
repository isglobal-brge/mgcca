// mgcca_scores_svd_rcpp — thin wrapper over mgcca::run_scores_svd (dual route:
// A_j = V diag(w_A) U'Y from SVD factors). Needs SVD/<ds>/{u,v,d} (run_svd first).
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.h"
using namespace Rcpp;

//' MGCCA dual scores stage (A=V diag(w_A) U'Y, weights, scores) over HDF5
//' @return A list with the stage descriptor (\code{filename}), or \code{NULL} on
//'   error. The per-table weights and scores are written into \code{final_group}
//'   of the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_scores_svd_rcpp(std::string filename, std::string tmp_group,
                                 std::vector<std::string> datasets, int inv_method,
                                 int nfac,
                                 Rcpp::Nullable<std::vector<double>> lambda = R_NilValue,
                                 std::string final_group = "FINAL_RESULTS")
{
    try {
        H5::Exception::dontPrint();
        std::vector<double> lam;
        if (inv_method == 2) {
            if (lambda.isNull()) throw std::runtime_error("penalized needs lambda");
            lam = Rcpp::as<std::vector<double>>(lambda);
        }
        mgcca::run_scores_svd(filename, tmp_group, datasets, inv_method, lam, nfac, final_group);
        return Rcpp::List::create(Rcpp::Named("filename") = filename);
    } catch (H5::Exception& e) {
        Rf_error("mgcca_scores_svd_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_scores_svd_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
