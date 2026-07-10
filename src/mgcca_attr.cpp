// mgcca_attr.cpp — Rcpp wrappers for reading/writing HDF5 attributes on a
// dataset or a group, used to persist and reload the mgcca provenance manifest.
//
// The heavy lifting lives in the generic, provider-neutral layer
// (hdf5Attributes.hpp) and the mgccaDataset subclass; here we only resolve the
// target object (dataset vs group) and marshal the named R list.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "hdf5Attributes.hpp"
#include "mgccaDataset.hpp"

using namespace Rcpp;

//' Write HDF5 attributes to a dataset or group
//'
//' Writes each named element of \code{attrs} as an HDF5 attribute. When
//' \code{dataset} is the empty string the attributes are attached to the group
//' \code{group}; otherwise they are attached to the dataset \code{group/dataset}.
//' Supports scalar and vector character/integer/double/logical values (logical
//' is stored as integer 0/1, since HDF5 has no native boolean).
//'
//' @param filename Path to the HDF5 file.
//' @param group Group path (target group, or the parent group of the dataset).
//' @param dataset Dataset name, or "" to target the group itself.
//' @param attrs Named list of attribute values.
//' @keywords internal
// [[Rcpp::export]]
void mgcca_write_attrs_rcpp(std::string filename, std::string group,
                            std::string dataset, Rcpp::List attrs) {
    H5::Exception::dontPrint();
    try {
        if (attrs.size() == 0) return;
        Rcpp::RObject nmObj = attrs.names();
        if (nmObj.isNULL())
            throw std::runtime_error("attrs must be a named list");
        Rcpp::CharacterVector nm(nmObj);

        if (dataset != "") {
            std::unique_ptr<mgcca::mgccaDataset> ds(
                new mgcca::mgccaDataset(filename, group, dataset, false));
            ds->openDataset();
            for (int i = 0; i < attrs.size(); ++i)
                ds->writeAttribute(Rcpp::as<std::string>(nm[i]), attrs[i]);
        } else {
            std::unique_ptr<BigDataStatMeth::hdf5File> f(
                new BigDataStatMeth::hdf5File(filename, false));
            f->openFile("rw");
            if (f->getFileptr() == nullptr)
                throw std::runtime_error("cannot open file " + filename);
            H5::Group g = f->getFileptr()->openGroup(group);
            for (int i = 0; i < attrs.size(); ++i)
                mgcca::attr::writeAttribute(g, Rcpp::as<std::string>(nm[i]),
                                            attrs[i]);
            g.close();
        }
    } catch (std::exception& ex) {
        Rf_error("mgcca_write_attrs_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_write_attrs_rcpp: unknown error");
    }
}

//' Read all HDF5 attributes from a dataset or group
//'
//' Returns every attribute (except BigDataStatMeth's internal bookkeeping
//' attribute) as a named list. When \code{dataset} is the empty string the
//' attributes are read from the group \code{group}; otherwise from the dataset
//' \code{group/dataset}.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group path (source group, or the parent group of the dataset).
//' @param dataset Dataset name, or "" to read from the group itself.
//' @return Named list of attribute values.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_read_attrs_rcpp(std::string filename, std::string group,
                                 std::string dataset) {
    H5::Exception::dontPrint();
    Rcpp::List out;
    try {
        if (dataset != "") {
            std::unique_ptr<mgcca::mgccaDataset> ds(
                new mgcca::mgccaDataset(filename, group, dataset, false));
            ds->openDataset();
            std::vector<std::string> names = ds->listAttributeNames();
            for (const auto& n : names) out[n] = ds->readAttribute(n);
        } else {
            std::unique_ptr<BigDataStatMeth::hdf5File> f(
                new BigDataStatMeth::hdf5File(filename, false));
            f->openFile("r");
            if (f->getFileptr() == nullptr)
                throw std::runtime_error("cannot open file " + filename);
            H5::Group g = f->getFileptr()->openGroup(group);
            std::vector<std::string> names = mgcca::attr::listAttributeNames(g);
            for (const auto& n : names)
                out[n] = mgcca::attr::readAttribute(g, n);
            g.close();
        }
    } catch (std::exception& ex) {
        Rf_error("mgcca_read_attrs_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_read_attrs_rcpp: unknown error");
    }
    return out;
}

//' List the immediate children of an HDF5 group
//'
//' Returns the names of the objects (datasets and subgroups) directly under
//' \code{group}, without depending on the rhdf5 package. Used by
//' \code{mgcca_load} to auto-discover table names in files written before the
//' provenance manifest existed.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group whose children to list.
//' @return Character vector of child object names.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::CharacterVector mgcca_list_group_rcpp(std::string filename,
                                            std::string group) {
    H5::Exception::dontPrint();
    std::vector<std::string> names;
    try {
        std::unique_ptr<BigDataStatMeth::hdf5File> f(
            new BigDataStatMeth::hdf5File(filename, false));
        f->openFile("r");
        if (f->getFileptr() == nullptr)
            throw std::runtime_error("cannot open file " + filename);
        H5::Group g = f->getFileptr()->openGroup(group);
        hsize_t n = g.getNumObjs();
        names.reserve((size_t)n);
        for (hsize_t i = 0; i < n; ++i)
            names.push_back(g.getObjnameByIdx(i));
        g.close();
    } catch (std::exception& ex) {
        Rf_error("mgcca_list_group_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_list_group_rcpp: unknown error");
    }
    return Rcpp::wrap(names);
}
