// mgcca_cor_ave_rcpp — thin [[Rcpp::export]] wrapper over mgcca::run_cor_ave.
// See src/mgcca_phases.hpp for the implementation.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.hpp"
using namespace Rcpp;

//' MGCCA corsY / p-values / AVE stage over HDF5
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_cor_ave_rcpp(std::string filename, std::string tmp_group,
                              std::vector<std::string> datasets, int nfac,
                              std::string final_group = "FINAL_RESULTS")
{
    try {
        H5::Exception::dontPrint();
        mgcca::run_cor_ave(filename, tmp_group, datasets, nfac, final_group);
        return Rcpp::List::create(Rcpp::Named("filename") = filename);
    } catch (H5::Exception& e) {
        Rf_error("mgcca_cor_ave_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_cor_ave_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
