// mgcca_eigen_rcpp — thin [[Rcpp::export]] wrapper over mgcca::run_eigen.
// See src/mgcca_phases.h for the implementation.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.h"
using namespace Rcpp;

//' MGCCA eigen stage (Ksum, MKsum05, eigen, Y) over HDF5
//' @return A list with the stage descriptor (\code{filename}, \code{nfac},
//'   \code{eig_values}, \code{Y_path}), or \code{NULL} on error. The shared
//'   components \code{Y} are written into \code{final_group} of the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_eigen_rcpp(std::string filename, std::string tmp_group,
                            std::vector<std::string> datasets, int nfac,
                            std::string final_group = "FINAL_RESULTS",
                            Rcpp::Nullable<int> threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        Rcpp::NumericVector ev = mgcca::run_eigen(filename, tmp_group, datasets,
                                                  nfac, final_group, threads);
        return Rcpp::List::create(
            Rcpp::Named("filename")   = filename,
            Rcpp::Named("nfac")       = nfac,
            Rcpp::Named("eig_values") = ev,
            Rcpp::Named("Y_path")     = final_group + "/Y");
    } catch (H5::Exception& e) {
        Rf_error("mgcca_eigen_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_eigen_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
