// mgcca_io.h — small shared helpers to read/write whole HDF5 datasets in
// R-view order using the BigDataStatMeth C++ API. Include after BigDataStatMeth.hpp.
#ifndef MGCCA_IO_HPP
#define MGCCA_IO_HPP

#include <BigDataStatMeth.hpp>

namespace mgcca {

// Read a full dataset into an Eigen matrix in R-view (column-major) order.
inline Eigen::MatrixXd read_full(const std::string& fn, const std::string& grp,
                                 const std::string& ds) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, false));
    d->openDataset();
    if (d->getDatasetptr() == nullptr)
        throw std::runtime_error("cannot open " + grp + "/" + ds);
    const std::size_t nr = d->nrows_r(), nc = d->ncols_r();
    Eigen::MatrixXd M(nr, nc);
    std::vector<hsize_t> off = {0, 0}, cnt = {(hsize_t)nc, (hsize_t)nr},
                         st = {1, 1}, bl = {1, 1};
    d->readDatasetBlock(off, cnt, st, bl, M.data());
    return M;
}

// Read the rownames (individual IDs) of a dataset.
inline Rcpp::CharacterVector read_rownames(const std::string& fn,
                                           const std::string& grp,
                                           const std::string& ds) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, false));
    d->openDataset();
    Rcpp::List dn = d->readDimnames();
    return dn["rownames"];
}

// Write an Eigen matrix (R-view m x n) to an HDF5 dataset, with optional dimnames.
// overwrite_ds=false makes the write FAIL-CLOSED (errors if the dataset already exists).
inline void write_full(const std::string& fn, const std::string& grp,
                       const std::string& ds, const Eigen::MatrixXd& M, int comp = 0,
                       Rcpp::Nullable<Rcpp::CharacterVector> rn = R_NilValue,
                       Rcpp::Nullable<Rcpp::CharacterVector> cn = R_NilValue,
                       bool overwrite_ds = true) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, overwrite_ds));
    d->setCompressionLevel(comp);
    d->createDataset((int)M.rows(), (int)M.cols(), "real");
    d->writeDataset(Rcpp::wrap(M));
    if (rn.isNotNull() || cn.isNotNull()) {
        BigDataStatMeth::hdf5Dims dims(d.get());
        Rcpp::CharacterVector rncv = rn.isNull() ? Rcpp::CharacterVector(0)
                                                 : Rcpp::CharacterVector(rn.get());
        Rcpp::CharacterVector cncv = cn.isNull() ? Rcpp::CharacterVector(0)
                                                 : Rcpp::CharacterVector(cn.get());
        dims.writeDimnames(rncv, cncv);
    }
}

// Write an Eigen matrix to HDF5, CREATING the file first if it does not exist
// (via hdf5File::createFile, the rcpp_hdf5_create_matrix pattern), then delegating
// to write_full. Keeps persistence in C++ (no R round-trip) and works whether or not
// the target file already exists.
// FAIL-CLOSED by default (overwrite_ds=false): errors if the target dataset already
// exists, so a confirmatory run can never silently reuse or overwrite a prior artifact.
inline void write_full_create(const std::string& fn, const std::string& grp,
                              const std::string& ds, const Eigen::MatrixXd& M, int comp = 0,
                              Rcpp::Nullable<Rcpp::CharacterVector> rn = R_NilValue,
                              Rcpp::Nullable<Rcpp::CharacterVector> cn = R_NilValue,
                              bool overwrite_ds = false) {
    {   // ensure the file exists on disk; RAII closes the handle before write_full reopens it
        std::unique_ptr<BigDataStatMeth::hdf5File> f(
            new BigDataStatMeth::hdf5File(fn, false));      // overwrite=false: do NOT truncate
        int iRes = f->createFile();                          // EXEC_OK created / EXEC_WARNING exists
        if (iRes != EXEC_OK && iRes != EXEC_WARNING)          // global constants (BigDataStatMeth.hpp)
            throw std::runtime_error("write_full_create: cannot create/open " + fn);
    }
    write_full(fn, grp, ds, M, comp, rn, cn, overwrite_ds);
}

}  // namespace mgcca
#endif  // MGCCA_IO_HPP
