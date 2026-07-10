// mgcca_phases.hpp — the mgcca pipeline as reusable inline stages, so both the
// per-stage [[Rcpp::export]] wrappers (for testing) and the single-call
// orchestrator (mgcca_rcpp) share one implementation. Each stage throws
// std::runtime_error on error; callers wrap in try/catch. Include after
// BigDataStatMeth.hpp and mgcca_io.hpp.
#ifndef MGCCA_PHASES_HPP
#define MGCCA_PHASES_HPP

#include <BigDataStatMeth.hpp>
#include "mgcca_io.hpp"
#include <algorithm>
#include <unordered_map>

namespace mgcca {

// ---- Phase 0 (optional): column-wise normalization, entirely in HDF5 --------
// Center+scale each variable (R-column) over its present individuals, streaming
// the table by blocks via BigDataStatMeth so no full copy is held in RAM. This
// reproduces base R scale(): the SD uses the (n-1) divisor (get_HDF5_mean_sd_*
// computes sqrt(sum((x-mu)^2)/(n-1)); the n/(n-1) RMS correction only applies
// when center=FALSE, which we never do here). MUST run on the RAW input (present
// individuals only), BEFORE getK pads missing individuals with zeros -- padding
// first would drag those zeros into the column means/SDs. Dimnames (rownames =
// individual IDs, required by getK; colnames = variable names, needed by
// scores/corsY) are copied from input to output.
inline void run_normalize(const std::string& filename, const std::string& in_group,
                          const std::string& out_group, const std::string& scale_group,
                          const std::vector<std::string>& datasets,
                          bool center, bool scale, Rcpp::Nullable<int> threads) {
    for (const std::string& ds : datasets) {
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, in_group, ds, false));
        dsA->openDataset();
        if (dsA->getDatasetptr() == nullptr)
            throw std::runtime_error("run_normalize: cannot open " + in_group + "/" + ds);

        // mean/sd per R-column (variable): datanormal is 2 x (HDF5 nrows = R ncols),
        // row 0 = mean, row 1 = sd. Column-wise == byrows=false.
        const hsize_t hnrows = dsA->nrows();
        Eigen::MatrixXd datanormal = Eigen::MatrixXd::Zero(2, hnrows);
        BigDataStatMeth::get_HDF5_mean_sd_by_column(dsA.get(), datanormal, true, true, R_NilValue);

        // Persist the per-variable centre/scale so predict() can transform new
        // data the same way: <scale_group>/<ds> is (p x 2) [center, scale].
        Rcpp::List dnA = dsA->readDimnames();
        Rcpp::CharacterVector varnames = dnA["colnames"];
        Eigen::MatrixXd cs(datanormal.cols(), 2);
        cs.col(0) = datanormal.row(0).transpose();     // center (mean)
        cs.col(1) = datanormal.row(1).transpose();     // scale  (sd)
        Rcpp::CharacterVector cscols = Rcpp::CharacterVector::create("center", "scale");
        write_full(filename, scale_group, ds, cs, 0, varnames, cscols);

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsOut(
            new BigDataStatMeth::hdf5Dataset(filename, out_group, ds, true));
        dsOut->setCompressionLevel(dsA->getCompressionLevel());
        dsOut->createDataset(dsA.get(), "real");
        if (dsOut->getDatasetptr() == nullptr)
            throw std::runtime_error("run_normalize: cannot create " + out_group + "/" + ds);

        BigDataStatMeth::RcppNormalizeHdf5(dsA.get(), dsOut.get(), datanormal, R_NilValue,
                                           center, scale, /*byrows*/false,
                                           /*bcorrected*/false, /*bparal*/false, threads);

        // Preserve dimnames (normalize does not copy them).
        BigDataStatMeth::hdf5Dims dims(dsOut.get());
        dims.writeDimnames(dnA["rownames"], dnA["colnames"]);
    }
}

// ---- Phase 1: getK — pad/reorder rows by individual ID, write X and K -------
// Returns the union of individual IDs (rn), in the order used for the rows.
inline Rcpp::CharacterVector run_getK(const std::string& filename,
                                      const std::string& in_group,
                                      const std::vector<std::string>& datasets,
                                      const std::string& out_group) {
    const std::size_t J = datasets.size();
    if (J < 2) throw std::runtime_error("need at least two tables");

    std::vector<std::vector<std::string>> ids(J);
    std::vector<std::size_t> ns(J);
    for (std::size_t j = 0; j < J; ++j) {
        Rcpp::CharacterVector rncv = read_rownames(filename, in_group, datasets[j]);
        ids[j].resize(rncv.size());
        for (R_xlen_t i = 0; i < rncv.size(); ++i) ids[j][i] = Rcpp::as<std::string>(rncv[i]);
        ns[j] = ids[j].size();
    }
    bool same_size = (*std::min_element(ns.begin(), ns.end()) ==
                      *std::max_element(ns.begin(), ns.end()));

    std::vector<std::string> rn;
    std::unordered_map<std::string, int> rn_pos;
    for (std::size_t j = 0; j < J; ++j)
        for (const auto& id : ids[j])
            if (rn_pos.emplace(id, (int)rn.size()).second) rn.push_back(id);
    if (!same_size) {
        std::sort(rn.begin(), rn.end());
        for (std::size_t k = 0; k < rn.size(); ++k) rn_pos[rn[k]] = (int)k;
    }
    const std::size_t m = rn.size();
    Rcpp::CharacterVector rn_cv(m);
    for (std::size_t k = 0; k < m; ++k) rn_cv[k] = rn[k];

    for (std::size_t j = 0; j < J; ++j) {
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsS(
            new BigDataStatMeth::hdf5Dataset(filename, in_group, datasets[j], false));
        dsS->openDataset();
        if (dsS->getDatasetptr() == nullptr)
            throw std::runtime_error("cannot open source " + datasets[j]);
        const std::size_t nj = dsS->nrows_r(), pj = dsS->ncols_r();
        Rcpp::List dnS = dsS->readDimnames();
        Rcpp::CharacterVector colnames = dnS["colnames"];

        // indicator diagonal K_j (depends only on presence, not on columns)
        Rcpp::NumericMatrix Kd((int)m, 1);
        for (std::size_t s = 0; s < nj; ++s) Kd(rn_pos[ids[j][s]], 0) = 1.0;

        // create padded X output (m x p_j) and stream the reorder COLUMN-BLOCK by
        // COLUMN-BLOCK so peak RAM is bounded to n_j x bcols, never the whole table.
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsX(
            new BigDataStatMeth::hdf5Dataset(filename, out_group + "/X", datasets[j], true));
        dsX->setCompressionLevel(dsS->getCompressionLevel());
        dsX->createDataset((int)m, (int)pj, "real");

        const std::size_t budget = 4000000;                 // ~32 MB of doubles
        const std::size_t bcols = std::max<std::size_t>(1, budget / std::max<std::size_t>(1, nj));
        for (std::size_t c0 = 0; c0 < pj; c0 += bcols) {
            const std::size_t bw = std::min(bcols, pj - c0);
            Rcpp::NumericMatrix S((int)nj, (int)bw);         // source block (n_j x bw)
            {   // read X_j[all rows, c0:c0+bw] (HDF5 axes swapped vs R)
                std::vector<hsize_t> off = {(hsize_t)c0, 0}, cnt = {(hsize_t)bw, (hsize_t)nj},
                                     st = {1, 1}, bl = {1, 1};
                dsS->readDatasetBlock(off, cnt, st, bl, S.begin());
            }
            Rcpp::NumericMatrix T((int)m, (int)bw);          // padded target block
            Eigen::Map<Eigen::MatrixXd> Sm(S.begin(), (Eigen::Index)nj, (Eigen::Index)bw);
            Eigen::Map<Eigen::MatrixXd> Tm(T.begin(), (Eigen::Index)m,  (Eigen::Index)bw);
            for (std::size_t s = 0; s < nj; ++s)
                Tm.row(rn_pos[ids[j][s]]) = Sm.row((Eigen::Index)s);
            {   // write the padded block into X output at cols c0:c0+bw
                std::vector<hsize_t> off = {(hsize_t)c0, 0}, cnt = {(hsize_t)bw, (hsize_t)m},
                                     st = {1, 1}, bl = {1, 1};
                dsX->writeDatasetBlock(Rcpp::wrap(T), off, cnt, st, bl, false);
            }
        }
        { BigDataStatMeth::hdf5Dims dims(dsX.get()); dims.writeDimnames(rn_cv, colnames); }

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsK(
            new BigDataStatMeth::hdf5Dataset(filename, out_group + "/K", datasets[j], true));
        dsK->setCompressionLevel(dsS->getCompressionLevel());
        dsK->createDataset((int)m, 1, "real");
        dsK->writeDataset(Rcpp::wrap(Kd));
    }
    return rn_cv;
}

// ---- Phase 1b: XKX — Mgram=X'X, inverse, Mi=X xkx X' ------------------------
// inv_method: 1=solve (SPD Cholesky), 2=penalized (+lambda I), 3=geninv/ginv (pseudoinv)
inline void run_XKX(const std::string& filename, const std::string& tmp_group,
                    const std::vector<std::string>& datasets, int inv_method,
                    const std::vector<double>& lam, Rcpp::Nullable<int> threads) {
    const std::size_t J = datasets.size();
    if (inv_method == 2 && lam.size() != J)
        throw std::runtime_error("penalized needs lambda of length #tables");
    const std::string gX = tmp_group + "/X";

    for (std::size_t j = 0; j < J; ++j) {
        const std::string& ds = datasets[j];
        {   // Mgram = X'X
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsX(
                new BigDataStatMeth::hdf5Dataset(filename, gX, ds, false)); dsX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsX2(
                new BigDataStatMeth::hdf5Dataset(filename, gX, ds, false)); dsX2->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsM(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/M", ds, true));
            dsM->setCompressionLevel(dsX->getCompressionLevel());
            int blk = BigDataStatMeth::getMaxBlockSize(
                dsX->nrows(), dsX->ncols(), dsX2->nrows(), dsX2->ncols(), 2, R_NilValue);
            BigDataStatMeth::crossprod(dsX.get(), dsX2.get(), dsM.get(),
                                       true, blk, 0, false, true, threads);
        }
        const std::string outXKXpath = tmp_group + "/XKX/" + ds;
        if (inv_method == 1 || inv_method == 2) {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsM(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/M", ds, false));
            dsM->openDataset();
            Rcpp::NumericVector d0;
            if (inv_method == 2) {
                d0 = BigDataStatMeth::getDiagonalfromMatrix(dsM.get());
                Rcpp::NumericVector d1 = Rcpp::clone(d0) + lam[j];
                BigDataStatMeth::setDiagonalMatrix(dsM.get(), d1);
            }
            std::unique_ptr<BigDataStatMeth::hdf5DatasetInternal> dsXKX(
                new BigDataStatMeth::hdf5DatasetInternal(filename, outXKXpath, true));
            dsXKX->setCompressionLevel(dsM->getCompressionLevel());
            dsXKX->createDataset((int)dsM->nrows(), (int)dsM->ncols(), "real");
            BigDataStatMeth::Rcpp_InvCholesky_hdf5(
                dsM.get(), dsXKX.get(), true, (long)MAXELEMSINBLOCK, threads);
            if (inv_method == 2) BigDataStatMeth::setDiagonalMatrix(dsM.get(), d0);
        } else if (inv_method == 3) {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsM(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/M", ds, false));
            dsM->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsXKX(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/XKX", ds, true));
            dsXKX->setCompressionLevel(dsM->getCompressionLevel());
            BigDataStatMeth::RcppPseudoinvHdf5(dsM.get(), dsXKX.get(), threads);
        } else {
            throw std::runtime_error("inv_method must be 1, 2 or 3");
        }
        {   // Mi = (X xkx) X'
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsX(
                new BigDataStatMeth::hdf5Dataset(filename, gX, ds, false)); dsX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsXKX(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/XKX", ds, false));
            dsXKX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsTmp(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/Xxkx", ds, true));
            dsTmp->setCompressionLevel(dsX->getCompressionLevel());
            BigDataStatMeth::multiplication(dsX.get(), dsXKX.get(), dsTmp.get(),
                                            false, false, R_NilValue, R_NilValue, threads);
        }
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsTmp(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/Xxkx", ds, false));
            dsTmp->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsX(
                new BigDataStatMeth::hdf5Dataset(filename, gX, ds, false)); dsX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsMi(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/Mi", ds, true));
            dsMi->setCompressionLevel(dsX->getCompressionLevel());
            BigDataStatMeth::multiplication(dsTmp.get(), dsX.get(), dsMi.get(),
                                            false, true, R_NilValue, R_NilValue, threads);
        }
    }
}

// ---- Phase 1b (DUAL): Mi from the thin SVD of X_j, without ever forming p x p.
// Small Gram matrix G = X_j X_j' (m x m, m = #individuals) via BigDataStatMeth's
// out-of-core tcrossprod, then its eigendecomposition (exact, in RAM since m is
// small). Left singular vectors of X_j are the eigenvectors U of G, with
// eigenvalues lambda_i = d_i^2. Then  M_j = U diag(w) U'  -- never forms p x p and
// never loads X_j fully (tcrossprod streams it). Weights:
//   solve/geninv (3): w_i = 1{lambda_i > tol}         (Moore-Penrose projection)
//   penalized  (2):   w_i = lambda_i / (lambda_i + lam_j)
// U and lambda are cached (tmp/U, tmp/lambda) for the scores stage. svd_method is
// accepted for API compatibility but the Gram route is exact & out-of-core.
inline void run_svd(const std::string& filename, const std::string& tmp_group,
                    const std::vector<std::string>& datasets, int inv_method,
                    const std::vector<double>& lam, const std::string& /*svd_method*/,
                    Rcpp::Nullable<int> threads) {
    const std::size_t J = datasets.size();
    if (inv_method == 2 && lam.size() != J)
        throw std::runtime_error("penalized needs lambda of length #tables");

    for (std::size_t j = 0; j < J; ++j) {
        const std::string& ds = datasets[j];
        // G = X X'  (m x m) via out-of-core tcrossprod (isSymmetric)
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/X", ds, false));
            dX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX2(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/X", ds, false));
            dX2->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dG(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/G", ds, true));
            dG->setCompressionLevel(dX->getCompressionLevel());
            int blk = BigDataStatMeth::getMaxBlockSize(
                dX->nrows(), dX->ncols(), dX2->nrows(), dX2->ncols(), 2, R_NilValue);
            BigDataStatMeth::tcrossprod(dX.get(), dX2.get(), dG.get(),
                                        true, blk, 0, false, true, threads);
        }
        Eigen::MatrixXd G = read_full(filename, tmp_group + "/G", ds);   // m x m (small)
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
        Eigen::VectorXd lambda = es.eigenvalues();          // ascending
        Eigen::MatrixXd U = es.eigenvectors();              // m x m
        const int m = (int)lambda.size();
        const double lmax = lambda.maxCoeff();
        const double tol = (lmax > 0 ? lmax : 1.0) * 1e-12;
        Eigen::VectorXd w(m);
        for (int i = 0; i < m; ++i) {
            const double li = lambda(i);
            w(i) = (inv_method == 2) ? li / (li + lam[j])
                                     : (li > tol ? 1.0 : 0.0);
        }
        Eigen::MatrixXd Mi = U * w.asDiagonal() * U.transpose();     // m x m (no p x p)
        write_full(filename, tmp_group + "/Mi", ds, Mi, 0);
        write_full(filename, tmp_group + "/U", ds, U, 0);            // cache for scores
        write_full(filename, tmp_group + "/lambda", ds,
                   Eigen::MatrixXd(lambda.transpose()), 0);
    }
}

// ---- Phase 2: Ksum, MKsum05, eigen, Y --------------------------------------
// Returns the eigenvalues.
inline Rcpp::NumericVector run_eigen(const std::string& filename,
                                     const std::string& tmp_group,
                                     const std::vector<std::string>& datasets,
                                     int nfac, const std::string& final_group,
                                     Rcpp::Nullable<int> threads) {
    const std::size_t J = datasets.size();
    Eigen::VectorXd Ksum;
    for (std::size_t j = 0; j < J; ++j) {
        Eigen::MatrixXd Kj = read_full(filename, tmp_group + "/K", datasets[j]);
        if (j == 0) Ksum = Kj.col(0); else Ksum += Kj.col(0);
    }
    Eigen::VectorXd Ksum05 = Ksum.array().pow(-0.5);

    Rcpp::CharacterVector rn = read_rownames(filename, tmp_group + "/X", datasets[0]);
    const std::size_t m = Ksum05.size();

    // Msum = sum_j M_j and MKsum05 = D Msum D built ROW-BLOCK by ROW-BLOCK, so
    // peak RAM is bounded to (brows x m) instead of the whole m x m matrix.
    std::vector<std::unique_ptr<BigDataStatMeth::hdf5Dataset>> dMi(J);
    for (std::size_t j = 0; j < J; ++j) {
        dMi[j].reset(new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/Mi", datasets[j], false));
        dMi[j]->openDataset();
    }
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> dMsum(
        new BigDataStatMeth::hdf5Dataset(filename, tmp_group, "Msum", true));
    dMsum->createDataset((int)m, (int)m, "real");
    std::unique_ptr<BigDataStatMeth::hdf5Dataset> dMK(
        new BigDataStatMeth::hdf5Dataset(filename, tmp_group, "MKsum05", true));
    dMK->createDataset((int)m, (int)m, "real");

    const std::size_t budget = 4000000;
    const std::size_t brows = std::max<std::size_t>(1, budget / std::max<std::size_t>(1, m));
    for (std::size_t r0 = 0; r0 < m; r0 += brows) {
        const std::size_t br = std::min(brows, m - r0);
        Eigen::MatrixXd buf = Eigen::MatrixXd::Zero((Eigen::Index)br, (Eigen::Index)m);
        for (std::size_t j = 0; j < J; ++j) {
            Rcpp::NumericMatrix Mb((int)br, (int)m);       // M_j[r0:r0+br, :]
            std::vector<hsize_t> off = {0, (hsize_t)r0}, cnt = {(hsize_t)m, (hsize_t)br},
                                 st = {1, 1}, bl = {1, 1};
            dMi[j]->readDatasetBlock(off, cnt, st, bl, Mb.begin());
            buf += Eigen::Map<Eigen::MatrixXd>(Mb.begin(), (Eigen::Index)br, (Eigen::Index)m);
        }
        std::vector<hsize_t> woff = {0, (hsize_t)r0}, wcnt = {(hsize_t)m, (hsize_t)br},
                             st = {1, 1}, bl = {1, 1};
        dMsum->writeDatasetBlock(Rcpp::wrap(buf), woff, wcnt, st, bl, false);
        // MKsum05[i,c] = Ksum05[r0+i] * buf(i,c) * Ksum05[c]
        Eigen::MatrixXd mk = Ksum05.segment((Eigen::Index)r0, (Eigen::Index)br).asDiagonal()
                             * buf * Ksum05.asDiagonal();
        dMK->writeDatasetBlock(Rcpp::wrap(mk), woff, wcnt, st, bl, false);
    }
    { BigDataStatMeth::hdf5Dims d(dMsum.get()); d.writeDimnames(rn, rn); }
    dMsum.reset(); dMK.reset();
    for (auto& d : dMi) d.reset();

    // The iterative (Spectra) eigensolver needs at least 2 requested eigenpairs
    // to build a valid Krylov subspace (it requires ncv > nev + 1); asking for a
    // single component (nfac == 1) makes it fail to converge. Request keig >= 2
    // and slice back to nfac, so nfac == 1 is supported without changing the
    // result for nfac >= 2 (the leading eigenpairs are identical).
    int keig = std::max(nfac, 2);
    BigDataStatMeth::RcppbdEigen_hdf5(filename, tmp_group, "MKsum05", keig, "LM", 0,
                                      false, false, 1e-10, 1000, true, true, threads);

    Eigen::MatrixXd V = read_full(filename, "EIGEN/MKsum05", "vectors");
    Eigen::MatrixXd Yast = V.leftCols(nfac);
    Eigen::MatrixXd Y = std::sqrt((double)J) * (Ksum05.asDiagonal() * Yast);
    Rcpp::CharacterVector comp(nfac);
    for (int c = 0; c < nfac; ++c) comp[c] = "comp" + std::to_string(c + 1);
    write_full(filename, final_group, "Y", Y, 0, rn, comp);

    Eigen::MatrixXd vals = read_full(filename, "EIGEN/MKsum05", "values");
    Eigen::Map<Eigen::VectorXd> vflat(vals.data(), vals.size());
    Eigen::VectorXd valsN = vflat.head(nfac);       // match the sliced Y columns
    return Rcpp::wrap(valsN);
}

// ---- Phase 3a: corsY (present individuals), p-values, AVE ------------------
inline void run_cor_ave(const std::string& filename, const std::string& tmp_group,
                        const std::vector<std::string>& datasets, int nfac,
                        const std::string& final_group) {
    const std::size_t J = datasets.size();
    Eigen::MatrixXd Y = read_full(filename, final_group, "Y");
    const std::size_t m = Y.rows();
    std::vector<int> p(J);
    Eigen::MatrixXd AVE_X(nfac, (int)J);

    auto standardise = [](Eigen::MatrixXd A) {
        for (int c = 0; c < A.cols(); ++c) {
            A.col(c).array() -= A.col(c).mean();
            double nrm = A.col(c).norm();
            if (nrm > 0) A.col(c) /= nrm;
        }
        return A;
    };

    const double df = (double)m - 2.0;
    Rcpp::CharacterVector comp(nfac);
    for (int c = 0; c < nfac; ++c) comp[c] = "comp" + std::to_string(c + 1);

    for (std::size_t j = 0; j < J; ++j) {
        const std::string& ds = datasets[j];
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsX(
            new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/X", ds, false));
        dsX->openDataset();
        const std::size_t pj = dsX->ncols_r();
        p[j] = (int)pj;
        Rcpp::CharacterVector varnames = dsX->readDimnames()["colnames"];

        Eigen::VectorXd Kd = read_full(filename, tmp_group + "/K", ds).col(0);  // m x 1, small
        std::vector<int> pres;
        for (std::size_t i = 0; i < m; ++i) if (Kd(i) != 0.0) pres.push_back((int)i);
        const int npr = (int)pres.size();
        Eigen::MatrixXd Yp(npr, nfac);                     // Y over present rows (small)
        for (int r = 0; r < npr; ++r) Yp.row(r) = Y.row(pres[r]);
        Eigen::MatrixXd YpN = standardise(Yp);

        // corsY / pval outputs (p_j x nfac); streamed column-block by column-block
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsCor(
            new BigDataStatMeth::hdf5Dataset(filename, final_group + "/corsY", ds, true));
        dsCor->createDataset((int)pj, nfac, "real");
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsPv(
            new BigDataStatMeth::hdf5Dataset(filename, final_group + "/pval", ds, true));
        dsPv->createDataset((int)pj, nfac, "real");

        for (int b = 0; b < nfac; ++b) AVE_X(b, (int)j) = 0.0;
        const std::size_t budget = 4000000;
        const std::size_t bcols = std::max<std::size_t>(1, budget / std::max<std::size_t>(1, m));
        for (std::size_t c0 = 0; c0 < pj; c0 += bcols) {
            const std::size_t bw = std::min(bcols, pj - c0);
            Rcpp::NumericMatrix Xb((int)m, (int)bw);        // block of variables (m x bw)
            {
                std::vector<hsize_t> off = {(hsize_t)c0, 0}, cnt = {(hsize_t)bw, (hsize_t)m},
                                     st = {1, 1}, bl = {1, 1};
                dsX->readDatasetBlock(off, cnt, st, bl, Xb.begin());
            }
            Eigen::Map<Eigen::MatrixXd> Xm(Xb.begin(), (Eigen::Index)m, (Eigen::Index)bw);
            Eigen::MatrixXd Xp(npr, (int)bw);               // present rows only
            for (int r = 0; r < npr; ++r) Xp.row(r) = Xm.row(pres[r]);
            Eigen::MatrixXd corsB = standardise(Xp).transpose() * YpN;   // bw x nfac
            Eigen::MatrixXd pvB((int)bw, nfac);
            for (int a = 0; a < (int)bw; ++a)
                for (int b = 0; b < nfac; ++b) {
                    double rr = corsB(a, b);
                    double t = rr * std::sqrt(df) / std::sqrt(1.0 - rr * rr);
                    pvB(a, b) = 2.0 * R::pt(std::fabs(t), df, 0, 0);
                    AVE_X(b, (int)j) += rr * rr;
                }
            std::vector<hsize_t> woff = {0, (hsize_t)c0}, wcnt = {(hsize_t)nfac, (hsize_t)bw},
                                 st = {1, 1}, bl = {1, 1};
            dsCor->writeDatasetBlock(Rcpp::wrap(corsB), woff, wcnt, st, bl, false);
            dsPv->writeDatasetBlock(Rcpp::wrap(pvB), woff, wcnt, st, bl, false);
        }
        for (int b = 0; b < nfac; ++b) AVE_X(b, (int)j) /= (double)pj;   // mean of corsY^2
        { BigDataStatMeth::hdf5Dims d1(dsCor.get()); d1.writeDimnames(varnames, comp); }
        { BigDataStatMeth::hdf5Dims d2(dsPv.get());  d2.writeDimnames(varnames, comp); }
    }
    double psum = 0; for (int j = 0; j < (int)J; ++j) psum += p[j];
    Eigen::VectorXd AVE_outer(nfac);
    for (int b = 0; b < nfac; ++b) {
        double s = 0; for (int j = 0; j < (int)J; ++j) s += p[j] * AVE_X(b, j);
        AVE_outer(b) = s / psum;
    }
    // The eigen stage may compute more eigenpairs than requested (>= 2, so the
    // solver converges even for nfac == 1); keep only the leading nfac.
    Eigen::MatrixXd AVE_inner_raw = read_full(filename, "EIGEN/MKsum05", "values");
    Eigen::Map<Eigen::VectorXd> aiflat(AVE_inner_raw.data(), AVE_inner_raw.size());
    Eigen::MatrixXd AVE_inner = aiflat.head(nfac);
    write_full(filename, final_group + "/AVE", "AVE_X", AVE_X, 0);
    write_full(filename, final_group + "/AVE", "AVE_outer", Eigen::MatrixXd(AVE_outer), 0);
    write_full(filename, final_group + "/AVE", "AVE_inner", AVE_inner, 0);
}

// ---- Phase 3b: scores (optional) ------------------------------------------
inline void run_scores(const std::string& filename, const std::string& tmp_group,
                       const std::vector<std::string>& datasets, int nfac,
                       const std::string& final_group) {
    const std::size_t J = datasets.size();
    Rcpp::CharacterVector comp(nfac);
    for (int c = 0; c < nfac; ++c) comp[c] = "comp" + std::to_string(c + 1);
    const std::string gX = tmp_group + "/X", gXKX = tmp_group + "/XKX",
                      Ypath = final_group + "/Y";

    // multiply helper (block-wise): C = op(A) %*% op(B)
    auto mult = [&](const std::string& gA, const std::string& dA,
                    const std::string& gB, const std::string& dB,
                    const std::string& gC, const std::string& dC,
                    bool tA, bool tB) {
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> a(
            new BigDataStatMeth::hdf5Dataset(filename, gA, dA, false)); a->openDataset();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> b(
            new BigDataStatMeth::hdf5Dataset(filename, gB, dB, false)); b->openDataset();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> c(
            new BigDataStatMeth::hdf5Dataset(filename, gC, dC, true));
        c->setCompressionLevel(a->getCompressionLevel());
        BigDataStatMeth::multiplication(a.get(), b.get(), c.get(), tA, tB,
                                        R_NilValue, R_NilValue, R_NilValue);
    };
    auto set_dimnames = [&](const std::string& g, const std::string& d,
                            Rcpp::CharacterVector rn, Rcpp::CharacterVector cn) {
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds(
            new BigDataStatMeth::hdf5Dataset(filename, g, d, false)); ds->openDataset();
        BigDataStatMeth::hdf5Dims dims(ds.get()); dims.writeDimnames(rn, cn);
    };

    for (std::size_t j = 0; j < J; ++j) {
        const std::string& ds = datasets[j];
        Rcpp::CharacterVector rn = read_rownames(filename, gX, ds);
        Rcpp::CharacterVector varnames;
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX(
                new BigDataStatMeth::hdf5Dataset(filename, gX, ds, false));
            dX->openDataset(); varnames = dX->readDimnames()["colnames"];
        }
        // XtY = X'Y (crossprod, block-wise)
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX(
                new BigDataStatMeth::hdf5Dataset(filename, gX, ds, false)); dX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dY(
                new BigDataStatMeth::hdf5Dataset(filename, final_group, "Y", false)); dY->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dXtY(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/XtY", ds, true));
            dXtY->setCompressionLevel(dX->getCompressionLevel());
            int blk = BigDataStatMeth::getMaxBlockSize(
                dX->nrows(), dX->ncols(), dY->nrows(), dY->ncols(), 2, R_NilValue);
            BigDataStatMeth::crossprod(dX.get(), dY.get(), dXtY.get(),
                                       false, blk, 0, false, true, R_NilValue);
        }
        mult(gXKX, ds, tmp_group + "/XtY", ds, final_group + "/A", ds, false, false); // A=xkx XtY
        mult(gX, ds, final_group + "/A", ds, tmp_group + "/KXA", ds, false, false);   // KXA=X A

        // small p x nfac / m x nfac matrices -> scale in RAM (nfac is tiny)
        Eigen::MatrixXd A   = read_full(filename, final_group + "/A", ds);
        Eigen::MatrixXd KXA = read_full(filename, tmp_group + "/KXA", ds);
        const double mm = (double)KXA.rows();
        Eigen::MatrixXd As = A;
        for (int b = 0; b < nfac; ++b) {
            double mean = KXA.col(b).mean();
            double sd = std::sqrt((KXA.col(b).array() - mean).square().sum() / (mm - 1.0));
            As.col(b) /= sd;
        }
        write_full(filename, final_group + "/weights", ds, As, 0, varnames, comp);
        mult(gX, ds, final_group + "/weights", ds, final_group + "/scores", ds, false, false);

        set_dimnames(final_group + "/A", ds, varnames, comp);
        set_dimnames(final_group + "/scores", ds, rn, comp);
    }
}

// ---- Phase 3b (DUAL): scores from the cached Gram eigen (tmp/U, tmp/lambda) -----
// A_j = X' (U diag(w_A) U' Y). The inner U diag(w_A) U' Y is m x nfac (tiny, RAM);
// A_j = crossprod(X, inner) is block-wise (p x nfac) -> no p x p, no full X_j in RAM.
// Then As/scores as usual (multiplication block-wise).
inline void run_scores_svd(const std::string& filename, const std::string& tmp_group,
                           const std::vector<std::string>& datasets, int inv_method,
                           const std::vector<double>& lam, int nfac,
                           const std::string& final_group) {
    const std::size_t J = datasets.size();
    Eigen::MatrixXd Y = read_full(filename, final_group, "Y");   // m x nfac (small)
    Rcpp::CharacterVector comp(nfac);
    for (int c = 0; c < nfac; ++c) comp[c] = "comp" + std::to_string(c + 1);

    auto mult = [&](const std::string& gA, const std::string& dA,
                    const std::string& gB, const std::string& dB,
                    const std::string& gC, const std::string& dC) {
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> a(
            new BigDataStatMeth::hdf5Dataset(filename, gA, dA, false)); a->openDataset();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> b(
            new BigDataStatMeth::hdf5Dataset(filename, gB, dB, false)); b->openDataset();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> c(
            new BigDataStatMeth::hdf5Dataset(filename, gC, dC, true));
        c->setCompressionLevel(a->getCompressionLevel());
        BigDataStatMeth::multiplication(a.get(), b.get(), c.get(), false, false,
                                        R_NilValue, R_NilValue, R_NilValue);
    };
    auto set_dimnames = [&](const std::string& g, const std::string& d,
                            Rcpp::CharacterVector rn, Rcpp::CharacterVector cn) {
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds(
            new BigDataStatMeth::hdf5Dataset(filename, g, d, false)); ds->openDataset();
        BigDataStatMeth::hdf5Dims dims(ds.get()); dims.writeDimnames(rn, cn);
    };

    for (std::size_t j = 0; j < J; ++j) {
        const std::string& ds = datasets[j];
        Rcpp::CharacterVector rn = read_rownames(filename, tmp_group + "/X", ds);
        Rcpp::CharacterVector varnames;
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/X", ds, false));
            dX->openDataset(); varnames = dX->readDimnames()["colnames"];
        }
        // inner = U diag(w_A) U' Y  (m x nfac, small) from the cached Gram eigen.
        // w_A: solve/geninv 1{lambda_i>tol}/lambda_i ; penalized 1/(lambda_i+lam_j).
        Eigen::MatrixXd U = read_full(filename, tmp_group + "/U", ds);       // m x m
        Eigen::MatrixXd L = read_full(filename, tmp_group + "/lambda", ds);  // 1 x m
        const int mm2 = (int)U.cols();
        double lmax = 0.0; for (int i = 0; i < mm2; ++i) lmax = std::max(lmax, L(0, i));
        const double tol = (lmax > 0 ? lmax : 1.0) * 1e-12;
        Eigen::VectorXd wA(mm2);
        for (int i = 0; i < mm2; ++i) {
            const double li = L(0, i);
            wA(i) = (inv_method == 2) ? 1.0 / (li + lam[j])
                                      : (li > tol ? 1.0 / li : 0.0);
        }
        Eigen::MatrixXd inner = U * wA.asDiagonal() * (U.transpose() * Y);   // m x nfac
        write_full(filename, tmp_group + "/inner", ds, inner, 0);           // small

        // A = X' inner = crossprod(X, inner) -> p x nfac (block-wise)
        {
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dX(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/X", ds, false));
            dX->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dI(
                new BigDataStatMeth::hdf5Dataset(filename, tmp_group + "/inner", ds, false));
            dI->openDataset();
            std::unique_ptr<BigDataStatMeth::hdf5Dataset> dA(
                new BigDataStatMeth::hdf5Dataset(filename, final_group + "/A", ds, true));
            dA->setCompressionLevel(dX->getCompressionLevel());
            int blk = BigDataStatMeth::getMaxBlockSize(
                dX->nrows(), dX->ncols(), dI->nrows(), dI->ncols(), 2, R_NilValue);
            BigDataStatMeth::crossprod(dX.get(), dI.get(), dA.get(),
                                       false, blk, 0, false, true, R_NilValue);
        }
        mult(tmp_group + "/X", ds, final_group + "/A", ds, tmp_group + "/KXA", ds); // KXA=X A (m x nfac)

        Eigen::MatrixXd A   = read_full(filename, final_group + "/A", ds);   // p x nfac (small)
        Eigen::MatrixXd KXA = read_full(filename, tmp_group + "/KXA", ds);   // m x nfac
        const double mm = (double)KXA.rows();
        Eigen::MatrixXd As = A;
        for (int b = 0; b < nfac; ++b) {
            double mean = KXA.col(b).mean();
            double sd = std::sqrt((KXA.col(b).array() - mean).square().sum() / (mm - 1.0));
            As.col(b) /= sd;
        }
        write_full(filename, final_group + "/weights", ds, As, 0, varnames, comp);
        mult(tmp_group + "/X", ds, final_group + "/weights", ds, final_group + "/scores", ds);

        set_dimnames(final_group + "/A", ds, varnames, comp);
        set_dimnames(final_group + "/scores", ds, rn, comp);
    }
}

}  // namespace mgcca
#endif  // MGCCA_PHASES_HPP
