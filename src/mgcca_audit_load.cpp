// mgcca_audit_load.cpp — reads a stored availability audit in ONE file open.
//
// The mirror of mgcca_audit_save.cpp. The save was already one open, one
// owner, one close; the load was not: the manifest came back through mgcca's
// attribute reader and each of the sixteen tables through BigDataStatMeth's
// R-level matrix reader, so reading one audit opened and shut the same file
// seventeen times under two statically linked HDF5 instances taking turns.
// Here the manifest and every table are read on one handle and the file is
// closed once at the end.
//
// WHAT COMES BACK DOES NOT CHANGE. The dimname rule is the one
// BigDataStatMeth's own reader applies: a dimension whose names the file does
// not carry gets none, and a table with neither carries no dimnames at all.
// The caller still turns a NaN that came back into NA, as it always did.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "hdf5Attributes.h"
#include "mgcca_io.h"

using namespace Rcpp;

//' Read a stored availability audit from one HDF5 file in a single open
//'
//' Reads the attributes of \code{group} and every table of \code{tables} that
//' the file carries, through one open of \code{filename}.
//'
//' @param filename Path to the HDF5 file holding the audit.
//' @param group Group the audit was written to.
//' @param tables Names of the audit tables to read. A name the file does not
//'   carry comes back as \code{NULL}, so the caller reports what is missing
//'   rather than failing on the first gap.
//' @return A list with \code{attrs} (the group's attributes) and \code{tables}
//'   (one matrix per table, with its dimnames), named after \code{tables}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mgcca_read_audit_rcpp(std::string filename, std::string group,
                                 std::vector<std::string> tables) {
    H5::Exception::dontPrint();
    Rcpp::List attrs, stored(tables.size());
    try {
        // Read-write, although nothing here writes: the tables used to come
        // back through BigDataStatMeth's R-level reader, which opens the file
        // read-write for reading too, and HDF5 refuses a second open of a file
        // this process already holds in the other mode.
        std::unique_ptr<BigDataStatMeth::hdf5File> file(
            new BigDataStatMeth::hdf5File(filename, false));
        file->openFile("rw");
        mgcca::check_readable(file.get(), group, "", "mgcca_read_audit_rcpp");

        // The manifest first, so a group that holds no audit is reported from
        // the attributes rather than from a table that happens to be missing.
        {
            H5::Group g = file->getFileptr()->openGroup(group);
            try {
                std::vector<std::string> names =
                    mgcca::attr::listAttributeNames(g);
                for (const auto& n : names)
                    attrs[n] = mgcca::attr::readAttribute(g, n);
            } catch (...) {
                g.close();
                throw;
            }
            g.close();
        }

        for (std::size_t k = 0; k < tables.size(); ++k) {
            if (!BigDataStatMeth::exists_HDF5_element(file->getFileptr(),
                                                      group + "/" + tables[k]))
                continue;
            Rcpp::NumericMatrix m =
                mgcca::read_matrix(file.get(), group, tables[k]);
            Rcpp::List dn = mgcca::read_dimnames(file.get(), group, tables[k]);
            Rcpp::CharacterVector rn = dn["rownames"], cn = dn["colnames"];
            if (rn.size() > 0 || cn.size() > 0)
                m.attr("dimnames") = Rcpp::List::create(
                    rn.size() > 0 ? Rcpp::RObject(rn) : Rcpp::RObject(R_NilValue),
                    cn.size() > 0 ? Rcpp::RObject(cn) : Rcpp::RObject(R_NilValue));
            stored[k] = m;
        }

        file.reset();

    } catch (H5::Exception& ex) {
        Rf_error("mgcca_read_audit_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_read_audit_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_read_audit_rcpp: unknown error");
    }
    stored.names() = Rcpp::CharacterVector(tables.begin(), tables.end());
    return Rcpp::List::create(Rcpp::Named("attrs") = attrs,
                              Rcpp::Named("tables") = stored);
}
