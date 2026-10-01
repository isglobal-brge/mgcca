// mgcca_audit_save.cpp — saves an availability audit in ONE file open.
//
// The audit save is the only mgcca operation that writes many datasets into
// one file in a single go: the sixteen encoded tables plus the group manifest.
// Driven from R it alternated TWO statically linked HDF5 instances over the
// same file -- the tables through BigDataStatMeth's R-level matrix creator,
// which runs inside BigDataStatMeth's own HDF5 instance, and the manifest
// through mgcca's attribute writer, which runs inside mgcca's -- so one save
// opened the file seventeen times under two owners that cannot see each
// other's handles. BigDataStatMeth::hdf5File::openFile() states the hazard
// itself: a file this process already holds through a different statically
// linked instance is invisible to its isOpenInCurrentProcess() probe. On
// Windows the second owner's open is then refused outright and the save fails
// part-way through, at no particular table.
//
// Here the whole save is one open, one owner, one close: the file is opened
// once, every table is created, filled and named through BigDataStatMeth's
// dataset class bound to that one open file, the manifest is written as group
// attributes on the same handle, and the handle is closed once at the end.
//
// WHAT IS STORED DOES NOT CHANGE. Every call below is the call
// BigDataStatMeth's own matrix creator makes, in the same order and with the
// same arguments: the same dataset type, the same default compression, the
// same dimnames convention, and the same treatment of a matrix whose rownames
// or colnames are absent. A file written here is the file the R-level loop
// wrote, and mgcca_audit_load() reads either.
//
// [[Rcpp::depends(BH, RcppEigen, Rhdf5lib, BigDataStatMeth)]]
#include <BigDataStatMeth.hpp>
#include "hdf5Attributes.h"

using namespace Rcpp;

namespace {

// One table of the audit, checked against R before anything is opened: the
// matrix as R holds it, its R-view shape, and the two dimname vectors to be
// stored beside it.
struct av_table {
    std::string name;
    Rcpp::RObject values;
    Rcpp::CharacterVector rownames;
    Rcpp::CharacterVector colnames;
    int nrow = 0;
    int ncol = 0;
    bool with_dimnames = false;
};

// Refuse an identifier BigDataStatMeth cannot store, BEFORE the file is open.
//
// The dataset writer reports an over-long identifier with Rf_error(), and an R
// error is a longjmp, not a C++ exception: it does not unwind the C++ stack, so
// no destructor runs and the file handle this frame owns stays open for the
// rest of the session, invisible to every later attempt to release it. The
// length is a property of the R object alone, so it is checked here, with
// nothing open; a failure is then an ordinary error raised from a frame that
// holds no handle. Same discipline as mgcca_io.h's check_fail_closed(): never
// let that longjmp fire from under our own handles.
void check_storable(const Rcpp::CharacterVector& v, const char* which,
                    const std::string& table) {
    const std::size_t limit = (std::size_t)(MAXSTRING - 1);   // global constant
    for (int i = 0; i < v.size(); ++i) {
        std::string s = Rcpp::as<std::string>(v[i]);
        if (s.size() > limit) {
            std::string preview = s.substr(0, 40);
            if (s.size() > 40) preview += "...";
            throw std::runtime_error(
                std::string(which) + " " +
                std::to_string((long)(i + 1)) + " of table '" + table +
                "' is " + std::to_string(s.size()) + " bytes ('" + preview +
                "'), over the " + std::to_string(limit) +
                " bytes HDF5 can store for an identifier; identifiers are "
                "never truncated, so shorten it");
        }
    }
}

// Read one table out of the R list and settle what will be stored beside it.
//
// The dimname rule is the one BigDataStatMeth's matrix creator applies, kept
// identical so the stored layout is identical: a vector shorter than the
// dimension it names is replaced by a single NA, which the dataset writer then
// stores only if the dimension really is of length one. That is how a table
// with no rownames ends up with a rownames entry when it has a single row, and
// with none when it has more.
av_table prepare(const std::string& name, Rcpp::RObject obj) {
    av_table t;
    t.name = name;

    if (TYPEOF(obj) != REALSXP || !Rf_isMatrix(obj))
        throw std::runtime_error("table '" + name +
                                 "' must be a double matrix");
    Rcpp::NumericMatrix m(obj);
    t.values = obj;
    t.nrow = m.nrow();
    t.ncol = m.ncol();

    Rcpp::List dn(obj.attr("dimnames"));
    if (dn.size() > 0) {
        t.with_dimnames = true;
        if (dn.size() > 0 && !Rf_isNull(dn[0]))
            t.rownames = Rcpp::CharacterVector(dn[0]);
        if (dn.size() > 1 && !Rf_isNull(dn[1]))
            t.colnames = Rcpp::CharacterVector(dn[1]);

        if (t.rownames.size() < t.nrow) {
            t.rownames = Rcpp::CharacterVector(1);          // a single NA
        } else if (t.colnames.size() < t.ncol) {
            t.colnames = Rcpp::CharacterVector(1);
        }
        // Only the vector whose length matches its dimension is stored, so
        // only that one is worth checking.
        if (t.rownames.size() == t.nrow)
            check_storable(t.rownames, "rowname", name);
        if (t.colnames.size() == t.ncol)
            check_storable(t.colnames, "colname", name);
    }
    return t;
}

}  // namespace

//' Save an availability audit into one HDF5 file in a single open
//'
//' Creates the audit group if it is missing, writes every table of
//' \code{tables} as a double dataset with its dimnames, and writes
//' \code{attrs} as attributes of the group. The file is opened once for the
//' whole operation and closed once at the end, so no part of the save depends
//' on reopening a file this process already holds.
//'
//' @param filename Path to the HDF5 file; it is created if it does not exist,
//'   and never truncated, so the blocks the audit describes survive.
//' @param group Group to write the audit into.
//' @param tables Named list of double matrices with dimnames, one per audit
//'   table. An existing dataset of the same name is replaced, whatever its
//'   shape.
//' @param attrs Named list of manifest values, written as attributes of
//'   \code{group}.
//' @return Invisibly \code{NULL}; called for the side effect of writing the
//'   audit into the HDF5 file.
//' @keywords internal
// [[Rcpp::export]]
void mgcca_save_audit_rcpp(std::string filename, std::string group,
                           Rcpp::List tables, Rcpp::List attrs) {
    H5::Exception::dontPrint();
    try {
        // ---- check the R side first, with nothing open ---------------------
        Rcpp::RObject tnmObj = tables.names();
        if (tables.size() > 0 && tnmObj.isNULL())
            throw std::runtime_error("`tables` must be a named list");
        Rcpp::CharacterVector tnm(tnmObj);
        std::vector<av_table> prepared;
        prepared.reserve((std::size_t)tables.size());
        for (int i = 0; i < tables.size(); ++i)
            prepared.push_back(prepare(Rcpp::as<std::string>(tnm[i]),
                                       Rcpp::RObject(tables[i])));

        Rcpp::RObject anmObj = attrs.names();
        if (attrs.size() > 0 && anmObj.isNULL())
            throw std::runtime_error("`attrs` must be a named list");
        Rcpp::CharacterVector anm(anmObj);

        // ---- one open ------------------------------------------------------
        // createFile() creates the file when it is absent and reports
        // EXEC_WARNING when it already exists, which is the case that matters:
        // the audit is written into the file holding the blocks it describes,
        // so the file is opened read-write and never truncated.
        std::unique_ptr<BigDataStatMeth::hdf5File> file(
            new BigDataStatMeth::hdf5File(filename, false));
        int iRes = file->createFile();
        if (iRes == EXEC_WARNING)
            file->openFile("rw");
        else if (iRes != EXEC_OK)
            throw std::runtime_error("cannot create or open " + filename);
        if (file->getFileptr() == nullptr)
            throw std::runtime_error("cannot open " + filename);

        // The audit group itself, created on this same handle when the file
        // does not carry one yet, so both the tables and the manifest below
        // find it there.
        {
            std::unique_ptr<BigDataStatMeth::hdf5Group> grp(
                new BigDataStatMeth::hdf5Group(file.get(), group));
            (void) grp->getGroupName();
        }

        // ---- the tables, all on that one open file --------------------------
        for (std::size_t k = 0; k < prepared.size(); ++k) {
            const av_table& t = prepared[k];
            // Bound to the open file, not to the path: the dataset object
            // takes the file handle it is given and does not own it, so it
            // closes only itself and leaves the file open for the next table.
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds(
                new BigDataStatMeth::hdf5Dataset(file.get(), group, t.name,
                                                 true));
            // Compression is left at BigDataStatMeth's default, which is what
            // its matrix creator used for these datasets.
            ds->createDataset((std::size_t)t.nrow, (std::size_t)t.ncol,
                              "numeric");
            if (ds->getDatasetptr() == nullptr)
                throw std::runtime_error("cannot create '" + group + "/" +
                                         t.name + "' in " + filename);
            ds->writeDataset(t.values);
            if (t.with_dimnames)
                ds->writeDimnames(t.rownames, t.colnames);
        }   // the dataset is closed here; the file stays open

        // ---- the manifest, on the same open file ----------------------------
        if (attrs.size() > 0) {
            H5::Group g = file->getFileptr()->openGroup(group);
            try {
                for (int i = 0; i < attrs.size(); ++i)
                    mgcca::attr::writeAttribute(
                        g, Rcpp::as<std::string>(anm[i]), attrs[i]);
            } catch (...) {
                g.close();
                throw;
            }
            g.close();
        }

        // ---- one close -----------------------------------------------------
        // Explicit, so the file is released here rather than at some later
        // point of this frame; on any other path out the same destructor runs
        // while the exception unwinds, before the error reaches R.
        file.reset();

    // Every handle above is held by a unique_ptr or by an H5::Group closed on
    // the way out, so by the time any of these handlers runs the file is shut:
    // the error reaches R with nothing of ours left open.
    } catch (H5::Exception& ex) {
        Rf_error("mgcca_save_audit_rcpp HDF5 error: %s",
                 ex.getDetailMsg().c_str());
    } catch (std::exception& ex) {
        Rf_error("mgcca_save_audit_rcpp: %s", ex.what());
    } catch (...) {
        Rf_error("mgcca_save_audit_rcpp: unknown error");
    }
}
