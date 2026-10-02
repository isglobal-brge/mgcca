// mgcca_read_datasets.cpp — readers that answer for EVERY dataset named, in
// ONE file open.
//
// Driven from R, a question about J datasets of one file was J questions: the
// R-level reader opens the file, hands back a value, and closes it again, so a
// loop over J blocks opened and shut the same file J times, and on Windows
// every one of those reopens is another chance for the file to be refused.
// Each entry point below opens the file once, binds every dataset it needs to
// that one handle, and closes it once at the end -- the same discipline as the
// audit writer in mgcca_audit_save.cpp.
//
// NOTHING HERE COMPUTES. These readers move bytes: they return values,
// dimnames and dimensions exactly as they are stored, and every arithmetic
// step stays where it was, in R.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "mgcca_io.h"

using namespace Rcpp;

namespace {

// Open the file for a whole batch of reads, and refuse a group the file does
// not carry before any dataset is bound to the handle.
//
// READ-WRITE, although nothing here writes. HDF5 refuses a second open of a
// file this process already holds in the other mode, and BigDataStatMeth's own
// reader opens a file read-write even to read it. Matching that mode is what
// keeps a read here from being refused a file that is perfectly available.
std::unique_ptr<BigDataStatMeth::hdf5File> open_for_reading(
        const std::string& filename, const std::string& group,
        const char* who) {
    std::unique_ptr<BigDataStatMeth::hdf5File> f(
        new BigDataStatMeth::hdf5File(filename, false));
    f->openFile("rw");
    mgcca::check_readable(f.get(), group, "", who);
    return f;
}

}  // namespace

//' Read several datasets of one HDF5 file in a single open
//'
//' Returns the values and the participant row names of every dataset named in
//' \code{datasets}, read through one open of \code{filename}.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group holding the datasets.
//' @param datasets Names of the datasets to read, in the order they are wanted.
//' @return A list with \code{values} (one matrix per dataset) and
//'   \code{rownames} (one character vector per dataset, empty where the file
//'   carries none), both named after \code{datasets}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_read_blocks_rcpp(std::string filename, std::string group,
                                  std::vector<std::string> datasets) {
    H5::Exception::dontPrint();
    Rcpp::List values(datasets.size()), rownames(datasets.size());
    try {
        std::unique_ptr<BigDataStatMeth::hdf5File> file =
            open_for_reading(filename, group, "mgcca_read_blocks_rcpp");
        for (std::size_t k = 0; k < datasets.size(); ++k) {
            values[k] = mgcca::read_matrix(file.get(), group, datasets[k]);
            rownames[k] = mgcca::read_rownames(file.get(), group, datasets[k]);
        }
        file.reset();
    } catch (H5::Exception& ex) {
        Rf_error("mgcca_read_blocks_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_read_blocks_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_read_blocks_rcpp: unknown error");
    }
    Rcpp::CharacterVector nm(datasets.begin(), datasets.end());
    values.names() = nm;
    rownames.names() = nm;
    return Rcpp::List::create(Rcpp::Named("values") = values,
                              Rcpp::Named("rownames") = rownames);
}

//' Read the row names of several datasets of one HDF5 file in a single open
//'
//' Reads the participant identifiers of every dataset named in
//' \code{datasets} without reading one value of any of them.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group holding the datasets.
//' @param datasets Names of the datasets whose row names are wanted.
//' @return A list of character vectors named after \code{datasets}, each empty
//'   where the file carries no row names for that dataset.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_read_rownames_rcpp(std::string filename, std::string group,
                                    std::vector<std::string> datasets) {
    H5::Exception::dontPrint();
    Rcpp::List out(datasets.size());
    try {
        std::unique_ptr<BigDataStatMeth::hdf5File> file =
            open_for_reading(filename, group, "mgcca_read_rownames_rcpp");
        for (std::size_t k = 0; k < datasets.size(); ++k)
            out[k] = mgcca::read_rownames(file.get(), group, datasets[k]);
        file.reset();
    } catch (H5::Exception& ex) {
        Rf_error("mgcca_read_rownames_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_read_rownames_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_read_rownames_rcpp: unknown error");
    }
    out.names() = Rcpp::CharacterVector(datasets.begin(), datasets.end());
    return out;
}

//' Read the dimensions of several datasets of one HDF5 file in a single open
//'
//' Reports the shape each dataset has as R sees it. A dataset the file does not
//' carry is reported as \code{0} by \code{0} rather than raised as an error, so
//' a caller sizing a choice over a set of datasets is not stopped by one that
//' is absent.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group holding the datasets.
//' @param datasets Names of the datasets to measure.
//' @return A list with the numeric vectors \code{nrow} and \code{ncol}, one
//'   entry per dataset, named after \code{datasets}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_read_dimensions_rcpp(std::string filename, std::string group,
                                      std::vector<std::string> datasets) {
    H5::Exception::dontPrint();
    Rcpp::NumericVector nrow(datasets.size()), ncol(datasets.size());
    try {
        std::unique_ptr<BigDataStatMeth::hdf5File> file =
            open_for_reading(filename, group, "mgcca_read_dimensions_rcpp");
        for (std::size_t k = 0; k < datasets.size(); ++k) {
            if (!BigDataStatMeth::exists_HDF5_element(
                    file->getFileptr(), group + "/" + datasets[k]))
                continue;
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> d = mgcca::open_on(
                file.get(), group, datasets[k], "mgcca_read_dimensions_rcpp");
            nrow[k] = (double) d->nrows_r();
            ncol[k] = (double) d->ncols_r();
        }
        file.reset();
    } catch (H5::Exception& ex) {
        Rf_error("mgcca_read_dimensions_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_read_dimensions_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_read_dimensions_rcpp: unknown error");
    }
    Rcpp::CharacterVector nm(datasets.begin(), datasets.end());
    nrow.names() = nm;
    ncol.names() = nm;
    return Rcpp::List::create(Rcpp::Named("nrow") = nrow,
                              Rcpp::Named("ncol") = ncol);
}
