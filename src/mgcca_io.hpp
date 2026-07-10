// mgcca_io.hpp — small shared helpers to read/write whole HDF5 datasets in
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
inline void write_full(const std::string& fn, const std::string& grp,
                       const std::string& ds, const Eigen::MatrixXd& M, int comp = 0,
                       Rcpp::Nullable<Rcpp::CharacterVector> rn = R_NilValue,
                       Rcpp::Nullable<Rcpp::CharacterVector> cn = R_NilValue) {
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> d(
        new BigDataStatMeth::hdf5Dataset(fn, grp, ds, true));
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

}  // namespace mgcca
#endif  // MGCCA_IO_HPP
