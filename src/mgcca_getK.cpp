// mgcca_getK_rcpp — thin [[Rcpp::export]] wrapper over mgcca::run_getK.
// See src/mgcca_phases.h for the implementation.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.h"
using namespace Rcpp;

//' MGCCA getK stage (build padded X and indicator diagonal K) over HDF5
//' @return A list with the stage descriptor (\code{filename}, \code{out_group},
//'   \code{datasets}, \code{m}, \code{rn}), or \code{NULL} on error. The
//'   union-padded tables and their indicator diagonals are written into
//'   \code{out_group} of the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_getK_rcpp(std::string filename, std::string in_group,
                           std::vector<std::string> datasets,
                           std::string out_group = "MGCCA_TMP")
{
    try {
        H5::Exception::dontPrint();
        Rcpp::CharacterVector rn = mgcca::run_getK(filename, in_group, datasets, out_group);
        return Rcpp::List::create(
            Rcpp::Named("filename")  = filename,
            Rcpp::Named("out_group") = out_group,
            Rcpp::Named("datasets")  = datasets,
            Rcpp::Named("m")         = (int)rn.size(),
            Rcpp::Named("rn")        = rn);
    } catch (H5::Exception& e) {
        Rf_error("mgcca_getK_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_getK_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
