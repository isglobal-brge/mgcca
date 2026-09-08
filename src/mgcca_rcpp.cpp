// mgcca_rcpp — single-call orchestrator for the full MGCCA pipeline over HDF5.
// Runs all stages in one C++ call (no inter-stage handle churn): getK -> XKX ->
// eigen/Y -> corsY/pval/AVE -> (optional) scores. Data must already be in the
// HDF5 file under `in_group` (individuals x variables, rownames = individual IDs)
// -- use mgcca_import_hdf5() on the R side to write any input type there.
// Results are written under `final_group`; intermediates under `tmp_group`.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_phases.h"
using namespace Rcpp;

//' MGCCA over HDF5 (single-call C++ orchestrator)
//'
//' @param filename HDF5 file with input tables already written under \code{in_group}.
//' @param in_group Input group (default \code{"MGCCA_IN"}).
//' @param datasets Character vector of dataset names inside \code{in_group}.
//' @param nfac Number of shared components.
//' @param inv_method 1 = solve (SPD Cholesky), 2 = penalized (+lambda I),
//'   3 = geninv/ginv (pseudoinverse).
//' @param lambda Numeric vector (length = #tables) for \code{inv_method = 2}.
//' @param scores If TRUE, also compute per-table weights and scores.
//' @param scale If TRUE (default), column-center+scale each table in-HDF5
//'   (out-of-core, base-R \code{scale()} semantics) before the analysis, reading
//'   the raw tables from \code{in_group} and never loading them fully into RAM.
//' @param route Per-table algebra route: \code{"auto"} (default) picks per table
//'   by \code{min(n, p)} -- covariance \code{X'X} (p x p) when \code{p < n},
//'   Gram \code{XX'} (n x n) when \code{p >= n}; \code{"cov"} forces the
//'   covariance route for all tables; \code{"dual"} forces the Gram route.
//' @param tmp_group Group for intermediates (default \code{"MGCCA_TMP"}).
//' @param final_group Group for results (default \code{"FINAL_RESULTS"}).
//' @param threads Optional thread count.
//' @return A descriptor list (filename, datasets, nfac, m, eigenvalues, route,
//'   route_dual, final_group). Results are written under \code{final_group}.
//' @seealso \code{\link{mgcca}}, \code{\link{mgcca_results}}
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_rcpp(std::string filename, std::string in_group,
                      std::vector<std::string> datasets, int nfac, int inv_method,
                      Rcpp::Nullable<std::vector<double>> lambda = R_NilValue,
                      bool scores = false, bool scale = true,
                      std::string route = "auto",
                      std::string tmp_group = "MGCCA_TMP",
                      std::string final_group = "FINAL_RESULTS",
                      Rcpp::Nullable<int> threads = R_NilValue)
{
    try {
        H5::Exception::dontPrint();
        if (route != "auto" && route != "cov" && route != "dual")
            throw std::runtime_error("route must be 'auto', 'cov' or 'dual'");
        std::vector<double> lam;
        if (inv_method == 2) {
            if (lambda.isNull())
                throw std::runtime_error("penalized (inv_method=2) requires lambda");
            lam = Rcpp::as<std::vector<double>>(lambda);
            if (lam.size() != datasets.size())
                throw std::runtime_error("lambda length must equal the number of tables");
        }

        // Optional column-wise scaling, done entirely in HDF5 on the RAW input
        // (present individuals only) BEFORE getK pads missing individuals with
        // zeros. Scaled tables land in tmp/IN; getK then reads from there.
        std::string getK_group = in_group;
        if (scale) {
            const std::string scaled_group = tmp_group + "/IN";
            mgcca::run_normalize(filename, in_group, scaled_group,
                                 final_group + "/scaling", datasets,
                                 /*center*/true, /*scale*/true, threads);
            getK_group = scaled_group;
        }

        Rcpp::CharacterVector rn =
            mgcca::run_getK(filename, getK_group, datasets, tmp_group);

        // Decide the route per table (on the padded X_j: n = #individuals = rows,
        // p = #variables = cols) and split the tables into two work lists. Mi lands
        // in tmp/Mi keyed by dataset regardless of source, so run_eigen/run_cor_ave
        // stay route-agnostic; only the XKX/scores stages dispatch.
        std::vector<std::string> cov_ds, dual_ds;
        std::vector<double> cov_lam, dual_lam;
        std::vector<char> is_dual(datasets.size(), 0);   // route flag per table
        for (std::size_t j = 0; j < datasets.size(); ++j) {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/X", datasets[j], false));
            dX->openDataset();
            const std::size_t n = dX->nrows_r(), p = dX->ncols_r();
            bool dual = (route == "dual") || (route == "auto" && p >= n);
            is_dual[j] = dual ? 1 : 0;
            if (dual) { dual_ds.push_back(datasets[j]); if (inv_method == 2) dual_lam.push_back(lam[j]); }
            else      { cov_ds.push_back(datasets[j]);  if (inv_method == 2) cov_lam.push_back(lam[j]); }
        }

        if (!cov_ds.empty())
            mgcca::run_XKX(filename, tmp_group, cov_ds, inv_method, cov_lam, threads);
        if (!dual_ds.empty())
            mgcca::run_svd(filename, tmp_group, dual_ds, inv_method, dual_lam, "auto", threads);

        Rcpp::NumericVector ev =
            mgcca::run_eigen(filename, tmp_group, datasets, nfac, final_group, threads);
        mgcca::run_cor_ave(filename, tmp_group, datasets, nfac, final_group);
        if (scores) {
            if (!cov_ds.empty())
                mgcca::run_scores(filename, tmp_group, cov_ds, nfac, final_group);
            if (!dual_ds.empty())
                mgcca::run_scores_svd(filename, tmp_group, dual_ds, inv_method,
                                      dual_lam, nfac, final_group);
        }

        Rcpp::LogicalVector route_dual(datasets.size());
        for (std::size_t j = 0; j < datasets.size(); ++j) route_dual[j] = (bool)is_dual[j];
        route_dual.names() = Rcpp::wrap(datasets);

        return Rcpp::List::create(
            Rcpp::Named("filename")    = filename,
            Rcpp::Named("datasets")    = datasets,
            Rcpp::Named("nfac")        = nfac,
            Rcpp::Named("m")           = (int)rn.size(),
            Rcpp::Named("eig_values")  = ev,
            Rcpp::Named("scores")      = scores,
            Rcpp::Named("route")       = route,
            Rcpp::Named("route_dual")  = route_dual,
            Rcpp::Named("final_group") = final_group);

    } catch (H5::Exception& e) {
        Rf_error("mgcca_rcpp HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("mgcca_rcpp error: %s", e.what());
    }
    return R_NilValue;
}
