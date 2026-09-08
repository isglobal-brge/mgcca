// mgcca_reliability_gram.cpp -- thin [[Rcpp::export]] wrapper over the K1 primitive
// mgcca::reliability::present_only_gram (PORT_FIRST_KERNEL.md). The kernel COMPUTES and
// RETURNS the n_pr x n_pr Gram (with participant dimnames) and OPTIONALLY persists it to
// HDF5 entirely in C++ via mgcca::write_full_create (creates the file if missing, the
// hdf5File::createFile pattern) -- no R round-trip.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"
#include "reliabilityGram.h"
using namespace Rcpp;

//' K1: present-only, feature-streamed participant block Gram (HDF5)
//' @return A list with the participant Gram \code{G} and its \code{ids}, the
//'   effective sizes (\code{p_eff}, \code{n_pr}, \code{N_all}), the streaming
//'   layout (\code{chunk}, \code{n_chunks}, \code{last_chunk},
//'   \code{eigen_threads}) and the \code{invariants} / \code{ledger}
//'   diagnostics of the accumulation.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List reliability_gram_hdf5(std::string file, std::string group, std::string dataset,
                                 std::vector<std::string> present_ids, int chunk,
                                 double var_eps = 1e-8,
                                 std::string out_file = "", std::string out_group = "",
                                 std::string out_dataset = "", int out_compression = 0) {
    try {
        H5::Exception::dontPrint();
        mgcca::reliability::GramResult R = mgcca::reliability::present_only_gram(
            file, group, dataset, present_ids, (long)chunk, var_eps);

        // G -> NumericMatrix with participant dimnames
        const int n = (int)R.n_pr;
        Rcpp::NumericMatrix G(n, n);
        std::copy(R.G.data(), R.G.data() + (std::size_t)n * n, G.begin());
        Rcpp::CharacterVector ids(R.ids.begin(), R.ids.end());
        G.attr("dimnames") = Rcpp::List::create(ids, ids);

        // optional C++ persistence (creates the file if missing; participant dimnames)
        if (!out_file.empty() && !out_dataset.empty())
            mgcca::write_full_create(out_file, out_group, out_dataset, R.G,
                                     out_compression, ids, ids);

        // per-chunk ledger -> data.frame-ready columns
        const R_xlen_t nc = (R_xlen_t)R.ledger.size();
        Rcpp::IntegerVector f_start(nc), f_end(nc), n_features(nc);
        Rcpp::NumericVector max_abs_row_sum(nc), max_abs_ss_dev(nc), trace_contrib(nc);
        for (R_xlen_t i = 0; i < nc; ++i) {
            const auto& r = R.ledger[i];
            f_start[i] = (int)r.f_start; f_end[i] = (int)r.f_end;
            n_features[i] = (int)r.n_features;
            max_abs_row_sum[i] = r.max_abs_row_sum;
            max_abs_ss_dev[i]  = r.max_abs_ss_dev;
            trace_contrib[i]   = r.trace_contrib;
        }
        Rcpp::DataFrame ledger = Rcpp::DataFrame::create(
            Rcpp::Named("f_start") = f_start, Rcpp::Named("f_end") = f_end,
            Rcpp::Named("n_features") = n_features,
            Rcpp::Named("max_abs_row_sum") = max_abs_row_sum,
            Rcpp::Named("max_abs_ss_dev") = max_abs_ss_dev,
            Rcpp::Named("trace_contrib") = trace_contrib);

        Rcpp::List invariants = Rcpp::List::create(
            Rcpp::Named("e_trace") = R.e_trace,
            Rcpp::Named("e_center") = R.e_center,
            Rcpp::Named("sym_err") = R.sym_err,
            Rcpp::Named("neg_mass") = R.neg_mass,
            Rcpp::Named("min_eig_ratio") = R.min_eig_ratio,
            Rcpp::Named("num_rank") = (double)R.num_rank);

        return Rcpp::List::create(
            Rcpp::Named("G") = G,
            Rcpp::Named("ids") = ids,
            Rcpp::Named("p_eff") = (double)R.p_eff,
            Rcpp::Named("n_pr") = (double)R.n_pr,
            Rcpp::Named("N_all") = (double)R.N_all,
            Rcpp::Named("chunk") = (double)R.chunk,
            Rcpp::Named("n_chunks") = R.n_chunks,
            Rcpp::Named("last_chunk") = (double)R.last_chunk,
            Rcpp::Named("eigen_threads") = R.eigen_threads,
            Rcpp::Named("invariants") = invariants,
            Rcpp::Named("ledger") = ledger);
    } catch (H5::Exception& e) {
        Rf_error("reliability_gram_hdf5 HDF5 error: %s", e.getDetailMsg().c_str());
    } catch (std::exception& e) {
        Rf_error("reliability_gram_hdf5 error: %s", e.what());
    }
    return R_NilValue;
}
