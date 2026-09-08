// mgcca_XKX_rcpp — thin [[Rcpp::export]] wrapper over mgcca::run_XKX.
// See src/mgcca_phases.h for the implementation.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.h"
using namespace Rcpp;

//' MGCCA XKX stage (Mgram=X'X, inverse, Mi=X xkx X') over HDF5
//' @return A list with the stage descriptor (\code{filename}, \code{tmp_group},
//'   \code{datasets}, \code{inv_method}), or \code{NULL} on error. The results
//'   themselves (per-table \code{Mgram}, its inverse and \code{Mi}) are written
//'   into \code{tmp_group} of the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_XKX_rcpp(std::string filename, std::string tmp_group,
                          std::vector<std::string> datasets, int inv_method,
                          Rcpp::Nullable<std::vector<double>> lambda = R_NilValue,
                          Rcpp::Nullable<int> threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        std::vector<double> lam;
        if (inv_method == 2) {
            if (lambda.isNull()) throw std::runtime_error("penalized needs lambda");
            lam = Rcpp::as<std::vector<double>>(lambda);
        }
        mgcca::run_XKX(filename, tmp_group, datasets, inv_method, lam, threads);
        return Rcpp::List::create(
            Rcpp::Named("filename") = filename,
            Rcpp::Named("tmp_group") = tmp_group,
            Rcpp::Named("datasets") = datasets,
            Rcpp::Named("inv_method") = inv_method);
    } catch (H5::Exception& e) {
        Rf_error("mgcca_XKX_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_XKX_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
