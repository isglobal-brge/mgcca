// mgcca_hold_file.cpp — hold ONE HDF5 file open for the whole of an analysis.
//
// Every entry point of this package opens the analysis file, does its whole
// job and closes it again, so two consecutive calls on one file are a close
// followed by a reopen -- and BigDataStatMeth runs a read-only pre-flight
// probe on each open, which is the call Windows refuses. That probe is skipped
// whenever the HDF5 instance doing the asking already holds the file open, so
// one handle owned here for the length of the analysis makes every open inside
// that window probe-free, leaving the entry points doing them untouched.
//
// NOTHING HERE COMPUTES, READS OR WRITES. The handle is opened, handed to R as
// an external pointer, and given up again.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>

using namespace Rcpp;

//' Open one HDF5 file and keep it open
//'
//' Opens \code{filename} and returns the open handle, which stays open until
//' it is released. The file is not read and not written.
//'
//' @param filename Path to the HDF5 file, exactly as every operation of the
//'   same analysis will name it.
//' @return An external pointer to the open file.
//' @keywords internal
// [[Rcpp::export]]
SEXP mgcca_hold_file_rcpp(std::string filename) {
    // The declared return type is SEXP rather than the pointer's own type
    // because the generated RcppExports.cpp repeats this signature and does
    // not include the BigDataStatMeth headers. What is returned is the
    // external pointer built at the end of this function.
    H5::Exception::dontPrint();
    // READ-WRITE, as every other open of this file in mgcca. HDF5 refuses a
    // second open of a file this process already holds in the other mode, and
    // the analyses that hold a handle do write to the file.
    std::unique_ptr<BigDataStatMeth::hdf5File> file(
        new BigDataStatMeth::hdf5File(filename, false));
    file->openFile("rw");
    if (file->getFileptr() == nullptr)
        Rf_error("mgcca_hold_file_rcpp: '%s' could not be opened",
                 filename.c_str());
    // The second argument registers the delete finalizer, so a handle that is
    // never released explicitly is still closed when R collects it.
    return Rcpp::XPtr<BigDataStatMeth::hdf5File>(file.release(), true);
}

//' Release a held HDF5 file handle
//'
//' Closes the file the handle holds and invalidates the handle. A handle that
//' holds nothing is left alone, so this can be called more than once.
//'
//' @param handle An external pointer returned by \code{mgcca_hold_file_rcpp}.
//' @return Invisibly \code{NULL}; called for the side effect of closing the
//'   file.
//' @keywords internal
// [[Rcpp::export]]
void mgcca_release_file_rcpp(SEXP handle) {
    if (handle == R_NilValue || TYPEOF(handle) != EXTPTRSXP) return;
    BigDataStatMeth::hdf5File* file =
        reinterpret_cast<BigDataStatMeth::hdf5File*>(R_ExternalPtrAddr(handle));
    if (file == nullptr) return;
    // Invalidated BEFORE the object goes, so a second call -- the garbage
    // collector's finalizer after an ordinary release, say -- finds nothing
    // left to give up.
    R_ClearExternalPtr(handle);
    try {
        delete file;   // ~hdf5File closes the file: openFile() took ownership
    } catch (...) {}
}

//' Whether a handle still holds an HDF5 file open
//'
//' @param handle An external pointer returned by \code{mgcca_hold_file_rcpp}.
//' @return \code{TRUE} while the handle holds an open file, \code{FALSE} once
//'   it has been released.
//' @keywords internal
// [[Rcpp::export]]
bool mgcca_file_is_open_rcpp(SEXP handle) {
    if (handle == R_NilValue || TYPEOF(handle) != EXTPTRSXP) return false;
    BigDataStatMeth::hdf5File* file =
        reinterpret_cast<BigDataStatMeth::hdf5File*>(R_ExternalPtrAddr(handle));
    return file != nullptr && file->getFileptr() != nullptr;
}
