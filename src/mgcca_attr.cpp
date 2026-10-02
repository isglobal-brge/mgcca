// mgcca_attr.cpp — Rcpp wrappers for reading/writing HDF5 attributes on a
// dataset or a group, used to persist and reload the mgcca provenance manifest.
//
// The heavy lifting lives in the generic, provider-neutral layer
// (hdf5Attributes.h) and the mgccaDataset subclass; here we only resolve the
// target object (dataset vs group) and marshal the named R list.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "hdf5Attributes.h"
#include "mgccaDataset.h"
#include "mgcca_io.h"

using namespace Rcpp;

namespace {

// The names of a list, refused when the list has values but no names.
Rcpp::CharacterVector named(const Rcpp::List& l, const char* what) {
    Rcpp::RObject nm = l.names();
    if (l.size() > 0 && nm.isNULL())
        throw std::runtime_error(std::string(what) + " must be a named list");
    return Rcpp::CharacterVector(nm);
}

// Write a named list as attributes of a group of an already open file.
void write_group_attrs(BigDataStatMeth::hdf5File* f, const std::string& group,
                       Rcpp::List attrs) {
    if (attrs.size() == 0) return;
    Rcpp::CharacterVector nm = named(attrs, "attrs");
    H5::Group g = f->getFileptr()->openGroup(group);
    try {
        for (int i = 0; i < attrs.size(); ++i)
            mgcca::attr::writeAttribute(g, Rcpp::as<std::string>(nm[i]),
                                        attrs[i]);
    } catch (...) {
        g.close();
        throw;
    }
    g.close();
}

// Every attribute of a group of an already open file, as a named list.
Rcpp::List read_group_attrs(BigDataStatMeth::hdf5File* f,
                            const std::string& group) {
    Rcpp::List out;
    H5::Group g = f->getFileptr()->openGroup(group);
    try {
        std::vector<std::string> names = mgcca::attr::listAttributeNames(g);
        for (const auto& n : names) out[n] = mgcca::attr::readAttribute(g, n);
    } catch (...) {
        g.close();
        throw;
    }
    g.close();
    return out;
}

}  // namespace

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
//' @return Invisibly \code{NULL}; called for the side effect of writing the
//'   attributes into the HDF5 file.
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

//' Write a provenance manifest into one HDF5 file in a single open
//'
//' Writes \code{attrs} as attributes of \code{group} and, for each element of
//' \code{dataset_attrs}, that element as attributes of the dataset named after
//' it under \code{dataset_group}. The run-level facts and the per-table facts
//' of one manifest therefore reach the file through a single open.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group the run-level attributes belong to.
//' @param attrs Named list of run-level attribute values.
//' @param dataset_group Group holding the datasets the per-table attributes
//'   belong to, or \code{""} when there are none.
//' @param dataset_attrs Named list, one named list of attribute values per
//'   dataset, named after the datasets.
//' @return Invisibly \code{NULL}; called for the side effect of writing the
//'   manifest into the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
void mgcca_write_manifest_rcpp(std::string filename, std::string group,
                               Rcpp::List attrs, std::string dataset_group,
                               Rcpp::List dataset_attrs) {
    H5::Exception::dontPrint();
    try {
        // Check the R side first, with nothing open: a list that cannot be
        // read is then refused from a frame that holds no handle.
        Rcpp::CharacterVector dnm = named(dataset_attrs, "dataset_attrs");
        std::vector<Rcpp::List> per((std::size_t)dataset_attrs.size());
        for (int i = 0; i < dataset_attrs.size(); ++i) {
            per[(std::size_t)i] = Rcpp::List(dataset_attrs[i]);
            named(per[(std::size_t)i], "dataset_attrs");
        }

        std::unique_ptr<BigDataStatMeth::hdf5File> file(
            new BigDataStatMeth::hdf5File(filename, false));
        file->openFile("rw");
        if (file->getFileptr() == nullptr)
            throw std::runtime_error("cannot open file " + filename);

        write_group_attrs(file.get(), group, attrs);

        for (int i = 0; i < dataset_attrs.size(); ++i) {
            const std::string ds = Rcpp::as<std::string>(dnm[i]);
            // Checked before the dataset object exists: its constructor would
            // CREATE the group it is given, and a manifest writer must never
            // add a group to the file it is describing.
            mgcca::check_readable(file.get(), dataset_group, ds,
                                  "mgcca_write_manifest_rcpp");
            std::unique_ptr<mgcca::mgccaDataset> d(
                new mgcca::mgccaDataset(file.get(), dataset_group, ds, false));
            d->openDataset();
            Rcpp::List a = per[(std::size_t)i];
            Rcpp::CharacterVector anm = a.names();
            for (int k = 0; k < a.size(); ++k)
                d->writeAttribute(Rcpp::as<std::string>(anm[k]), a[k]);
        }

        file.reset();

    } catch (H5::Exception& ex) {
        Rf_error("mgcca_write_manifest_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_write_manifest_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_write_manifest_rcpp: unknown error");
    }
}

//' Read a provenance manifest from one HDF5 file in a single open
//'
//' Reads the attributes of \code{group} and the attributes of every dataset
//' under \code{dataset_group}, through one open of \code{filename}. A
//' \code{dataset_group} the file does not carry yields no per-dataset
//' attributes rather than an error, since a file can hold a run-level manifest
//' and no per-table one.
//'
//' @param filename Path to the HDF5 file.
//' @param group Group the run-level attributes belong to.
//' @param dataset_group Group holding the datasets whose attributes are wanted,
//'   or \code{""} to read none.
//' @return A list with \code{group} (the run-level attributes) and
//'   \code{datasets} (one named list of attribute values per dataset, named
//'   after the datasets).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_read_manifest_rcpp(std::string filename, std::string group,
                                    std::string dataset_group) {
    H5::Exception::dontPrint();
    Rcpp::List run, per;
    try {
        std::unique_ptr<BigDataStatMeth::hdf5File> file(
            new BigDataStatMeth::hdf5File(filename, false));
        file->openFile("r");
        mgcca::check_readable(file.get(), group, "",
                              "mgcca_read_manifest_rcpp");

        run = read_group_attrs(file.get(), group);

        if (!dataset_group.empty() &&
            BigDataStatMeth::exists_HDF5_element(file->getFileptr(),
                                                 dataset_group)) {
            H5::Group g = file->getFileptr()->openGroup(dataset_group);
            std::vector<std::string> kids;
            try {
                hsize_t n = g.getNumObjs();
                for (hsize_t i = 0; i < n; ++i) {
                    const std::string nm = g.getObjnameByIdx(i);
                    // Only the datasets: the dimnames of a table are stored
                    // beside it as a group, and a group carries no manifest.
                    if (g.childObjType(nm) == H5O_TYPE_DATASET)
                        kids.push_back(nm);
                }
            } catch (...) {
                g.close();
                throw;
            }
            g.close();

            for (const auto& ds : kids) {
                std::unique_ptr<mgcca::mgccaDataset> d(
                    new mgcca::mgccaDataset(file.get(), dataset_group, ds,
                                            false));
                d->openDataset();
                Rcpp::List a;
                std::vector<std::string> names = d->listAttributeNames();
                for (const auto& n : names) a[n] = d->readAttribute(n);
                per[ds] = a;
            }
        }

        file.reset();

    } catch (H5::Exception& ex) {
        Rf_error("mgcca_read_manifest_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_read_manifest_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_read_manifest_rcpp: unknown error");
    }
    return Rcpp::List::create(Rcpp::Named("group") = run,
                              Rcpp::Named("datasets") = per);
}
