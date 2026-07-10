// mgccaDataset.hpp — a subclass of BigDataStatMeth::hdf5Dataset that adds
// attribute I/O by delegating to the generic (provider-neutral) mgcca::attr
// functions in hdf5Attributes.hpp.
//
// PURPOSE: demonstrate that hdf5Dataset can be extended with attribute support
// from an external package (mgcca), without modifying BigDataStatMeth. These
// methods are the exact hook that could later fold into hdf5Dataset inside BDSM.
//
// Include AFTER BigDataStatMeth.hpp.
#ifndef MGCCA_MGCCADATASET_HPP
#define MGCCA_MGCCADATASET_HPP

#include <string>
#include <vector>

#include <BigDataStatMeth.hpp>
#include "hdf5Attributes.hpp"

namespace mgcca {

class mgccaDataset : public BigDataStatMeth::hdf5Dataset {
public:
    using BigDataStatMeth::hdf5Dataset::hdf5Dataset;  // inherit constructors

    // Write / update a single attribute on this dataset.
    void writeAttribute(const std::string& name, Rcpp::RObject value) {
        ensureOpen();
        mgcca::attr::writeAttribute(*getDatasetptr(), name, value);
    }

    // Read a single attribute from this dataset.
    Rcpp::RObject readAttribute(const std::string& name) {
        ensureOpen();
        return mgcca::attr::readAttribute(*getDatasetptr(), name);
    }

    // Names of all attributes on this dataset (dropping BDSM's "internal").
    std::vector<std::string> listAttributeNames(bool skip_internal = true) {
        ensureOpen();
        return mgcca::attr::listAttributeNames(*getDatasetptr(), skip_internal);
    }

    // True if this dataset carries an attribute named `name`.
    bool hasAttribute(const std::string& name) {
        ensureOpen();
        return mgcca::attr::hasAttribute(*getDatasetptr(), name);
    }

private:
    // Attribute ops need an open H5::DataSet handle; caller must openDataset().
    void ensureOpen() {
        if (getDatasetptr() == nullptr)
            throw std::runtime_error(
                "mgccaDataset: dataset not open; call openDataset() first");
    }
};

}  // namespace mgcca
#endif  // MGCCA_MGCCADATASET_HPP
