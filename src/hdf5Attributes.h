// hdf5Attributes.h — generic HDF5 attribute I/O for any object that can carry
// attributes (dataset, group, or file), built on top of the HDF5 C++ (H5Cpp)
// API that BigDataStatMeth already links against (Rhdf5lib).
//
// DESIGN: these are provider-neutral free functions operating on H5::H5Object&,
// the common base of H5::DataSet, H5::Group and H5::H5File (all three expose
// createAttribute/openAttribute/getNumAttrs/attrExists in HDF5 1.10). They carry
// no mgcca-specific logic, so they are meant to migrate verbatim into
// BigDataStatMeth's hdf5Dataset/hdf5Group/hdf5File hierarchy in the future.
//
// Supported R value types (scalar and vector; a scalar is just a length-1
// vector — every attribute is stored rank-1 with dims = {length}):
//   character            -> variable-length HDF5 string (interop with rhdf5)
//   logical              -> NATIVE_INT (0/1); HDF5 has no native bool
//   integer              -> NATIVE_INT
//   numeric/double       -> NATIVE_DOUBLE
//
// Include AFTER BigDataStatMeth.hpp.
#ifndef MGCCA_HDF5ATTRIBUTES_HPP
#define MGCCA_HDF5ATTRIBUTES_HPP

#include <algorithm>
#include <string>
#include <vector>

#include <BigDataStatMeth.hpp>

namespace mgcca {
namespace attr {

// ---------------------------------------------------------------------------
// Write a single attribute, dispatching on the R type of `value`. If an
// attribute with the same name already exists it is removed first, so writes
// are idempotent and may change type/length.
// ---------------------------------------------------------------------------
inline void writeAttribute(H5::H5Object& obj, const std::string& name,
                           Rcpp::RObject value) {
    if (obj.attrExists(name)) obj.removeAttr(name);

    // Order matters: Rf_isNumeric() is TRUE for integers too, so check the
    // narrower predicates (logical, integer) before the numeric fallback.
    if (Rf_isString(value)) {
        Rcpp::CharacterVector v(value);
        hsize_t n = (hsize_t)v.size();
        H5::DataSpace space(1, &n);
        H5::StrType stype(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<std::string> owner(n);
        std::vector<const char*> cptr(n);
        for (hsize_t i = 0; i < n; ++i) {
            owner[i] = Rcpp::as<std::string>(v[i]);
            cptr[i] = owner[i].c_str();
        }
        H5::Attribute a = obj.createAttribute(name, stype, space);
        a.write(stype, cptr.data());
        a.close();

    } else if (Rf_isLogical(value)) {
        Rcpp::LogicalVector v(value);
        hsize_t n = (hsize_t)v.size();
        std::vector<int> buf(v.begin(), v.end());
        H5::DataSpace space(1, &n);
        H5::Attribute a =
            obj.createAttribute(name, H5::PredType::NATIVE_INT, space);
        a.write(H5::PredType::NATIVE_INT, buf.data());
        a.close();

    } else if (Rf_isInteger(value)) {
        Rcpp::IntegerVector v(value);
        hsize_t n = (hsize_t)v.size();
        std::vector<int> buf(v.begin(), v.end());
        H5::DataSpace space(1, &n);
        H5::Attribute a =
            obj.createAttribute(name, H5::PredType::NATIVE_INT, space);
        a.write(H5::PredType::NATIVE_INT, buf.data());
        a.close();

    } else if (Rf_isNumeric(value)) {
        Rcpp::NumericVector v(value);
        hsize_t n = (hsize_t)v.size();
        std::vector<double> buf(v.begin(), v.end());
        H5::DataSpace space(1, &n);
        H5::Attribute a =
            obj.createAttribute(name, H5::PredType::NATIVE_DOUBLE, space);
        a.write(H5::PredType::NATIVE_DOUBLE, buf.data());
        a.close();

    } else {
        throw std::runtime_error(
            "writeAttribute: unsupported R type for attribute '" + name + "'");
    }
}

// ---------------------------------------------------------------------------
// Read a single attribute into the matching R vector type. Length-1 attributes
// come back as length-1 vectors (natural R scalars). Handles both
// variable-length and fixed-length HDF5 strings.
// ---------------------------------------------------------------------------
inline Rcpp::RObject readAttribute(H5::H5Object& obj, const std::string& name) {
    H5::Attribute a = obj.openAttribute(name);
    H5::DataType dtype = a.getDataType();
    H5T_class_t cls = dtype.getClass();
    H5::DataSpace space = a.getSpace();
    hsize_t n = (hsize_t)space.getSimpleExtentNpoints();

    if (cls == H5T_STRING) {
        H5::StrType stype = a.getStrType();
        Rcpp::CharacterVector out(n);
        if (stype.isVariableStr()) {
            std::vector<char*> buf(n, nullptr);
            a.read(stype, buf.data());
            for (hsize_t i = 0; i < n; ++i)
                out[i] = buf[i] ? buf[i] : "";
            H5::DataSet::vlenReclaim(buf.data(), stype, space);
        } else {
            size_t sz = stype.getSize();
            std::vector<char> buf(n * sz, '\0');
            a.read(stype, buf.data());
            for (hsize_t i = 0; i < n; ++i) {
                const char* p = &buf[i * sz];
                size_t len = 0;
                while (len < sz && p[len] != '\0') ++len;
                out[i] = std::string(p, len);
            }
        }
        a.close();
        return out;

    } else if (cls == H5T_INTEGER) {
        std::vector<int> buf(n);
        a.read(H5::PredType::NATIVE_INT, buf.data());
        a.close();
        return Rcpp::IntegerVector(buf.begin(), buf.end());

    } else if (cls == H5T_FLOAT) {
        std::vector<double> buf(n);
        a.read(H5::PredType::NATIVE_DOUBLE, buf.data());
        a.close();
        return Rcpp::NumericVector(buf.begin(), buf.end());

    } else {
        a.close();
        throw std::runtime_error(
            "readAttribute: unsupported HDF5 type for attribute '" + name + "'");
    }
}

// ---------------------------------------------------------------------------
// Names of all attributes attached to `obj`. `skip_internal` drops the
// "internal" bookkeeping attribute that BigDataStatMeth stamps on every dataset
// it creates, so callers get a clean, user-facing manifest.
// ---------------------------------------------------------------------------
inline std::vector<std::string> listAttributeNames(H5::H5Object& obj,
                                                   bool skip_internal = true) {
    int num = obj.getNumAttrs();
    std::vector<std::string> names;
    names.reserve(num > 0 ? (size_t)num : 0);
    for (int i = 0; i < num; ++i) {
        H5::Attribute a = obj.openAttribute((unsigned int)i);
        std::string nm = a.getName();
        a.close();
        if (skip_internal && nm == "internal") continue;
        names.push_back(nm);
    }
    return names;
}

// True if `obj` carries an attribute named `name`.
inline bool hasAttribute(H5::H5Object& obj, const std::string& name) {
    return obj.attrExists(name);
}

}  // namespace attr
}  // namespace mgcca
#endif  // MGCCA_HDF5ATTRIBUTES_HPP
