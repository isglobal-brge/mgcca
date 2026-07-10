# mgcca 1.0.0

Complete, ground-up rewrite. The numerical core of Generalized Canonical
Correlation Analysis with missing individuals now runs **out-of-core in C++ on
HDF5** through the current [BigDataStatMeth](https://github.com/isglobal-brge/BigDataStatMeth)
API, scaling from small tables to full-omics data without loading a whole table
into memory or forming any variable-by-variable matrix.

## Highlights

- **Single entry point `mgcca()`** — imports the data into HDF5, scales in-file,
  runs the whole pipeline out-of-core, and returns a self-contained `"mgcca"`
  object (or a lightweight descriptor with `collect = FALSE`).
- **Scales to `p >> n`.** A per-table algebra route (`route = "auto"`) uses the
  Gram/dual formulation for wide tables, so integrating e.g. full methylation
  (~485k features) never forms a `p x p` matrix and completes in seconds.
- **Missing individuals** handled analytically via per-table indicator matrices,
  so every individual in the union is represented without imputation.
- **Inversion methods:** `solve`, `penalized` (per-table `lambda`),
  `geninv`/`ginv` (Moore–Penrose).

## New features

- `mgcca_load()` reopens a full result **from its HDF5 file alone**. Results now
  carry a self-describing provenance manifest (method, lambda, route, number of
  components, tables, eigenvalues, version) stored as native HDF5 attributes.
- `mgcca_results()` collects the numeric results into an `"mgcca"` object.
- `mgcca_associate()` (regress components on phenotypes; R²/p with BH
  adjustment) and `mgcca_permtest()` (permutation significance).
- `predict()` for out-of-sample projection using the stored weights and training
  scaling.
- Modern ggplot2 plots: `plotIndividuals()` (with confidence ellipses),
  `plotVariables()` (correlation circle), `plotLoadings()`, `plotScores()`,
  `plotBiplot()`, `plotAVE()`, `plotScree()`.
- S3 methods `print`, `summary`, `plot`, `predict` for `"mgcca"` objects.

## Under the hood

- Numerical core in C++ on the BigDataStatMeth C++/HDF5 API (crossproducts,
  eigendecomposition, pseudo-inverse, normalisation — all block-wise/out-of-core).
- A generic HDF5 attribute layer (scalar/vector, on datasets, groups and files)
  that extends BigDataStatMeth's `hdf5Dataset` without modifying it.
- `testthat` suite (152 tests) validated against the previous R implementation
  (0.9.9) as the reference oracle; a full-scale TCGA-ACC benchmark
  (methylation 485,577 × 80 + RNA-seq) runs in ~30 s.
- Documentation regenerated with roxygen2; pkgdown site; runnable vignette.

## Breaking changes

- The API is new. The 0.9.x functions `mgcca_bd()`, `mgcca_hdf5()`, `cv_gcca()`,
  `looe()` and the various in-memory helpers are removed or replaced; use
  `mgcca()` and the reader/analysis/plot functions above.
- Requires a current **BigDataStatMeth** and an HDF5 backend (`Rhdf5lib`); data
  is processed on HDF5 files rather than in memory.


# mgcca 0.9.9

Last release of the original in-memory / R-orchestrated implementation
(2018-09-03), archived for reference and reproducibility of earlier analyses.

- Core GCCA with missing individuals in R (`mgcca()`), plus `mgcca_bd()`
  (in-RAM) and `mgcca_hdf5()` (HDF5 orchestration in R against the earlier
  BigDataStatMeth API).
- Cross-validation (`cv_gcca()`, `looe()`), scores/weights helpers, import from
  `MultiAssayExperiment` (`getTables()`), base-R plots (`plotInds()`,
  `plotVars()`), significance helpers (`getSignif()`, `topVars()`).
- Requires each table to fit in memory; does not scale to full-omics tables.
  Superseded by 1.0.0.
