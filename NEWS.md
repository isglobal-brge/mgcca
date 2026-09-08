# mgcca 0.99.0

**Bioconductor submission.** The estimator is untouched: the C++/HDF5 core and
every number a fit produces are exactly those of the internal 1.4.0 release.
Besides packaging, metadata and documentation, this version carries one
addition that runs *beside* the estimator rather than inside it -- a pre-fit
audit of the availability design, which never calls `mgcca()` and cannot move
any quantity a fit computes (internal lineage 1.5.0, see *New* below).

## New

- **`mgcca_audit()` audits an availability design before anything is fitted.**
  It reads who is present in which block -- and nothing else -- and reports
  four things: how many individuals each block and each declared stratum holds;
  which pairs of blocks share enough individuals to be supported directly;
  which pairs are reachable only through other blocks, by which route, and how
  much support the weakest edge of that route carries; and how concentrated any
  declared availability weights are, per stratum, as a Kish effective sample
  size beside the raw counts.

  **It takes the object the fit will consume.** Besides a ready-made
  availability matrix, it accepts the same inputs `mgcca()` takes in either
  mode -- a list of blocks, a `MultiAssayExperiment`, a `SummarizedExperiment`,
  an `ExpressionSet`, an HDF5 file or its descriptor -- and derives the
  availability table itself through the package's own readers, with the same
  union-of-IDs rule the estimator uses. Only the individual IDs are read; no
  block data is materialised on either route.

  **It reports, it does not rule.** The pair codes (`DIRECT_SUPPORT`,
  `MODEL_TRANSFER_ELIGIBLE`, `UNSUPPORTED`, `STRUCTURAL_ZERO`) describe the
  design's support, not identification: a pair reachable only through a route
  is identified conditional on the declared rank-L model and its required
  structural conditions, which the audit does not check -- `$structural_checks`
  says `"not_evaluated"` and always will at this level, so that a consumer
  cannot mistake "not checked" for "checked and fine". The three design lines
  (`ess_min`, `availability_min`, `max_weight`) are declared lines with margin,
  configurable, and the value that raised a code is always reported beside it.
  Nothing is a verdict, and no design makes the audit refuse to return one.

  The edge connectivity it reports is **design-support fragility**: how many
  independent routes the design gives you between two blocks, and which single
  edge would cut them apart if it went. It is a statement about the design's
  topology, not about how well anything is estimated. Likewise, a joint count
  reaching `L` is descriptive support and no more -- it does not prove that the
  corresponding cross-block quantity has numerical rank `L` -- and a finite
  entry in the table evidences positivity for that cell without proving it, an
  empty cell being empty by chance and impossible by design looking exactly the
  same from here.

- **`mgcca_audit_load()` reads a stored audit back.** With
  `mgcca_audit(..., save_hdf5 = TRUE)` on an HDF5-backed design, the complete
  audit is written into an `availability_audit` group of that same file, so the
  certificate travels with the data. The encoding is structured and inspectable
  rather than a blob: numeric tables as datasets, block names, individual IDs
  and stratum labels as their dimnames, status codes as integers with the
  frozen vocabulary stored beside them as a legend, transfer routes as
  sequences of integer block indices. It goes in and comes out through
  BigDataStatMeth and the package's own attribute layer, and the loader checks
  what it rebuilds against what was stored rather than trusting either alone.

- **Six ways to look at an audit.** `print()` and `summary()` methods come with
  the class, and `plot()` defaults to `type = "overview"`: a curated four-panel
  composite with the pair classification as its anchor, realised availability,
  representation, and a small support graph, composed through base `grid`
  viewports and returning its panels invisibly. The five single types --
  `"matrix"`, `"graph"`, `"availability"`, `"representation"` and `"patterns"`
  -- come back as `ggplot` objects. `"patterns"` is an UpSet-like count of
  which combinations of blocks participants actually carry, aggregated at plot
  time from the availability table the object already holds.
  `plot(a, type = "graph", pair = c(...))` draws the stored route a
  transfer-eligible pair depends on, marks its bottleneck, and says in the
  subtitle that this holds under the declared rank-L assumption.

  **The grammar is fixed on purpose.** Colour encodes one variable per panel,
  and in the matrix that variable is the frozen classification -- never how much
  support a pair has, so a thinly supported direct pair is the same colour as a
  richly supported one with its count printed in the cell. There is no
  traffic-light scale and no red or green anywhere: a design has properties, not
  grades. Edge width in the graph is the raw count of jointly observed
  participants and nothing else, and says so in its legend and its caption; what
  the graph shows is design-support redundancy, which is a fact about topology,
  not about conditioning. The availability fill is pinned to `[0, 1]` rather
  than the observed range, representation is two aligned panels rather than one
  dual axis, and the declared design lines are drawn as neutral dashed lines
  labelled as declared. Every palette is Okabe-Ito and colour-blind-safe.

- **A vignette section walks through it** (*Auditing the availability design
  before fitting*, in `mgcca_example`): a four-block design where one block is
  reachable only through a single edge, the route and its bottleneck, the
  representation of the declared strata, the support graph, and a
  `save_hdf5`/`mgcca_audit_load()` round trip. All six plot types appear, each
  with a plain-language paragraph on what to look at first, what every visual
  element means, and what it decides for the reader.

## Unchanged

- No estimator code path was touched and no default changed. The audit is
  additive and estimator-independent: a fit made with this version is the fit
  1.4.0 made, which the 148-key identity bridge against the sealed 1.4.0
  reference confirms.
- No new dependency: the audit's `grid` and `ggplot2` are already declared
  `Imports`, and its optional `MultiAssayExperiment` / `Biobase` inputs are
  already declared `Suggests`.

## Version renumbering

- Bioconductor requires a new submission to carry version **0.99.0**, which the
  project auto-bumps to 1.0.0 at the first Bioconductor release. The package's
  internal lineage (1.0.0 -> 1.4.0, recorded below, plus the audit wave above,
  internally 1.5.0) is kept in this file as provenance: **0.99.0 is that
  lineage renumbered**, not a rollback. The full test suite passes unchanged
  against the 1.4.0 reference values, and the audit module above carries its
  own tests on top.

## Submission hygiene

- `DESCRIPTION`: added `biocViews`; `Depends` reduced to `R (>= 4.5.0)` and
  `BigDataStatMeth (>= 2.0.3)`, with `MASS` and `parallel` moved to `Imports`
  alongside the explicitly declared `grid`, `methods`, `stats` and `utils`;
  `LazyData` removed (the data set is loaded explicitly with
  `data(cardiovascular)`); minimum R version raised to 4.5.0.
- `LICENSE` rewritten as the standard MIT stub (`YEAR` / `COPYRIGHT HOLDER`).
- `cor.test.p()` is now exported as a plain function. It was registered as an S3
  method for `cor.test`, which is not a generic, so the namespace could not be
  loaded with only the base namespace attached. Its computation is unchanged.
- `print.listMAE()` now has the S3-consistent signature `(x, ...)`. The printed
  output is unchanged.
- `useDynLib(mgcca, .registration = TRUE)`.
- The low-level orchestrator `mgcca_rcpp()` is no longer exported; it remains
  available as an internal entry point and is still reached through `mgcca()`.
- The three data objects `X1`, `X2` and `X3` shipped in `data(cardiovascular)`
  are now documented.
- Every man page documents its return value, and the user-facing pages carry
  runnable examples.
- `src/`: the private headers were renamed `*.hpp` -> `*.h` (with their
  `#include`s) so `R CMD check` recognises them as source files; no code change.
  `src/Makevars.win` no longer overrides the compiler's optimisation flags and
  now mirrors `src/Makevars`.
- Both vignettes gained an Installation section.
- The heaviest integration tests are marked `skip_on_bioc()` so the Bioconductor
  builders stay inside the check time limit; they still run everywhere else.
- Style sweeps requested by `BiocCheck` (`TRUE`/`FALSE` instead of `T`/`F`,
  `inherits()` instead of `class(x) ==`, `vapply()`, `seq_len()`/`seq_along()`,
  `<-` for assignment), all behaviour-neutral.

# mgcca 1.4.0

Three additions that widen how a fit is chosen, collected and run, and **no
change to any number the estimator produces on the HDF5 path**. The C++ core is
untouched, the file-backed route returns exactly what 1.3.0 returned, and the
two new call paths (a penalty selector and an in-memory backend) sit beside the
existing one without displacing it. See the regression note at the end.

## New

- **`mgcca_select_lambda()` chooses the ridge penalty for
  `method = "penalized"` by held-out cross-block agreement.** It runs K-fold over
  individuals, projects the held-out individuals through each block's weights and
  scores how well the per-block scores agree across blocks (mean pairwise RV),
  returning the chosen `lambda`, the agreement diagnostic and the whole candidate
  grid so the choice is on the record instead of buried in a fit. Choosing
  `lambda` by minimising predicted instability is the intuitive alternative and
  it is wrong: it lands in the over-regularised regime, where each block has
  fallen onto its own non-shared axis and nothing shared is left to move.
  Cross-block agreement is a concordance criterion, and a concordance guard warns
  and returns the boundary candidate — a `mgcca_lambda_boundary` warning — when
  the chosen `lambda` sits at an end of the grid, so a grid that does not bracket
  the optimum is visible rather than silent.

  **It is available, never automatic.** `mgcca(method = "penalized")` still
  requires an explicit `lambda` and still stops without one. No default changes.

- **`mgcca()` and `mgcca_results()` take `outputs =`,** naming the subset of
  components to collect (`"Y"`, `"corsY"`, `"scores"`, `"pval"`, `"weights"`,
  `"scaling"`, `"AVE"`, `"eigen"`, `"overlap"`); the datasets behind the rest are
  never read. The choice is **purely subtractive**: with the default (`NULL`) the
  collected object is bit-for-bit what earlier versions returned, and with a
  subset every component asked for is `identical()` to its full-collection
  counterpart. Selecting outputs changes what is *read*, never what was
  *computed*. `print()` and `summary()` now print only the blocks the object
  actually carries, so a partially collected fit prints instead of failing; a
  fully collected fit prints exactly what it printed before.

- **In-memory estimation with `filename = NULL`.** `mgcca()` now defaults
  `filename` to `NULL`; passing a file is unchanged in every respect, and
  omitting it fits **entirely in RAM** for tables that fit in RAM, with no HDF5
  file created. The in-memory path accepts a list of matrices, a MultiAssayExperiment,
  a SummarizedExperiment or an ExpressionSet; a path to an existing HDF5 file
  still requires `filename` and is refused clearly without one. It covers the
  four methods, the `auto`/`cov`/`dual` routes, `scores`, the overlap report,
  the spectral diagnostics and `outputs =`, and returns an `"mgcca"` object of
  the same shape whose descriptor reports `backend = "memory"` and
  `filename = NA`. `print`, `summary`, `plot*`, `predict`, `mgcca_associate` and
  `mgcca_permtest` work on it.

  The reliability layer still requires a file-backed fit. `mgcca_sensitivity()`
  and `mgcca_stability()` on an in-memory fit stop with a message pointing to the
  HDF5-backed source blocks, and `mgcca_load()` and `collect = FALSE` do not
  apply without a file and say so. The in-memory `print()` marks the fit as not
  reloadable from a file.

## Parity

- **The in-memory route reproduces the pure-R mathematical reference to machine
  precision** (measured `2e-16` against the covariance oracle). Its cross-backend
  equivalence with the HDF5 route is certified at the package's own established
  `1e-5` standard (the one already used by `test-route.R` and the `eqs`/`eqm`
  helpers). One observation is documented rather than smoothed over: on
  ill-conditioned data the **HDF5 route** departs from the pure-R reference by up
  to `~2e-7` in `weights` (and `~6e-10` in `pval`), so the RAM-to-HDF5 difference
  on those fixtures is larger than the difference of the RAM route from the
  reference. That gap is a pre-existing property of the sealed binary's
  block-Cholesky inverse and eigen solver, **not** of this version; the in-memory
  route is the more accurate of the two, matching the reference exactly. The full
  measured per-quantity table is archived with the release pin.

## Regression

- **No estimator, kernel or computed quantity changes on the HDF5 path: `src/`
  is byte-for-byte the 1.2.0 source.** Every output the file-backed route
  produced is returned unchanged, bridged `identical()` against a 1.3.0 build on
  the shared fixtures (`paper_PLOS/pin_bridge_130_140.R`): the default fits across
  the seven fixture x method combinations, the components, the descriptor and the
  sensitivity/stability surface all match bit-for-bit, `mgcca_results(outputs =)`
  subsets are `identical()` to their full counterparts, and `mgcca_select_lambda`
  is present in 1.4.0 and absent in 1.3.0 — the change is additive only. The
  bridge reports `bridge_ok = TRUE` with `max|Delta| = 0` over all pre-existing
  quantities.

# mgcca 1.3.0

One quantity the pipeline already computed but discarded, now exposed. **No
estimator, kernel or computed quantity changes**: `src/` is byte-for-byte the
1.2.0 source, and every output 1.2.0 returned is returned unchanged (see the
regression note at the end).

- **Per-participant sensitivity contributions: `by_individual`.** The sealed C++
  kernel already returned per-participant accumulators on every
  `mgcca_sensitivity()` call and the R layer threw them away. They are now
  reported as a long-format data frame with columns `query`, `id`, `block`,
  `S_total`, `share_total` and — only when `group` is supplied — `S_between` and
  `share_between`. Writing `s_ij` for a participant's contribution and `S_j` for
  the block quantity already in `by_block`, the identity is `sum_i s_ij = S_j`:
  it is the same sensitivity cut over participants instead of over blocks, not a
  new estimand.

  The rows are exactly the ordered participant set used by the grouped
  calculation, so the contributions are **conditioned on that set** — a grouping
  covering fewer participants gives different numbers, which is a different
  question rather than an instability. **Key the table by `query`, `id` and
  `block`; row order is not part of the contract.** Shares are taken within each
  query and block (`share_total = s_ij / S_j`), so they sum to one there, and a
  denominator that is not finite and strictly positive gives `NA` rather than a
  silent zero. A participant absent from a block contributes a numerical zero to
  it (of order 1e-34), not an exact zero, and is deliberately not thresholded.

  Note that `S_between` is built from group means and is therefore constant
  across the members of a group: it is the group's between-group contribution
  carried in each member's row. `S_total` is the per-participant quantity.

  **This decomposition is exploratory.** Its block-level parent has empirical
  calibration evidence with a MIXTO verdict; **no individual-level calibration
  has been established**, and that remains future work. The values are
  sensitivity contributions — not variance, not uncertainty, not causal
  influence, and not guarantees of correctness or reliability. The absolute
  contributions inherit the non-identified scale of the parent quantity; the
  shares remove that common multiplicative scale and are the preferred
  dimensionless representation, but they stay conditional on the fitted model,
  the ridge, the query, the block and the grouped participant set. No invariant
  ranking of participants is implied.

- **`by_individual_check`.** A new element reporting the decomposition identity
  actually achieved: `sum_rel_total`, and `sum_rel_between` when grouped, are the
  largest `|sum_i s_ij - S_j| / max(1, |S_j|)` over queries and blocks. The
  identity is definitional but not bit-for-bit, because a trace and a sum over
  rows accumulate in different orders. It is deliberately **not** placed inside
  `validity`, which stays exactly what earlier versions produced.

- **New view `plot(s, type = "individuals", top_n = 15)`.** One panel per query,
  one horizontal bar per participant, stacked by block, ranked by contribution
  summed over blocks, with the participant's group in the axis label when a
  grouping was supplied. Everyone outside `top_n` is aggregated into a single
  **"all others" bar** rather than dropped, so each panel accounts for the whole
  of the query's sensitivity: a top-N chart without the remainder silently
  overstates concentration, and the remainder is what lets a reader see how much
  the named bars leave out.

- **New vignette section**, "Which participants is it coming from?", in
  *How much of a shared subspace is carried by group structure?*. It extracts the
  table, shows `summary()` and the new plot on the `miniACC` grouping, and states
  the exploratory framing and the conditioning in prose rather than only in a
  comment.

## Bug fix

- **`print()`, `summary()` and `plot()` looked up `by_individual` with `$`,
  which partial-matches.** On a result carrying `by_individual_check` but no
  `by_individual`, `x$by_individual` returned the check element instead of
  `NULL`, so the guard for older objects never fired and the failure surfaced
  later as an unreadable type error. All four lookups now use exact `[[`
  indexing. Results written by released versions carry neither element and were
  never affected.

## Regression

- **Every pre-existing output is unchanged.** `overall`, `by_block`,
  `validity`, `settings`, `resamples` and `summary_metrics` are bridged
  `identical()` against a 1.2.0 build on a shared fixture
  (`paper_PLOS/pin_bridge_120_130.R`), and the 1.2.0 object is confirmed to lack
  `by_individual` and `by_individual_check` while the 1.3.0 object carries them —
  the change is additive only.

# mgcca 1.2.0

Two diagnostics that **report on a fit without changing it**. No computed
quantity moves: the estimator, the C++ core and every number a 1.1.0 fit
returned are untouched, and both diagnostics are read back from datasets the
pipeline had already written.

- **Breakdown warning for `method = "solve"`.** The eigen stage already computes
  the external spectral gap of the mgcca operator, `mu_L - mu_{L+1}`. That gap
  is what separates *the solver converged* from *the leading subspace is
  identified*, and the two come apart in a regime that occurs in practice: when
  a table carries about as many variables as it has observed rows, the top
  eigenvalues coincide, the gap falls to machine zero, and the directions the
  solver returns are arbitrary — while the solver still reports success. A fit
  in that state now raises a warning of class `"mgcca_degenerate_subspace"`
  naming the mechanism and recommending `method = "penalized"`, which restores
  the separation. **It is a warning and nothing else**: the fit is returned
  exactly as computed, no result is altered and nothing is aborted. The screen
  is `gap_rel < 1e-8` (with an absolute safety net at `1e-10`). The two regimes
  observed so far are separated by more than nine orders of magnitude of empty
  space — machine zero on one side, `2.3e-4` at the very lowest on the other —
  so the cut has about five orders of clearance above the broken side and four
  below the healthy side and does not depend on where exactly it is placed.

- **The spectral gap is now visible.** `fit$eigen` (and, for
  `collect = FALSE`, `attr(fit, "desc")$eigen`) carries the absolute gap, the
  relative gap and the eigen-residual, and `print()` shows the relative gap on
  its own line, marked `NOT IDENTIFIED` when it is degenerate.

- **Pairwise overlap report in `summary()`.** `fit$overlap` gives, for every
  pair of tables, the number of individuals they share (`n_jk`) and its
  normalised version `alpha_jk = n_jk / n` over the union, together with the
  per-table counts and the dispersion of `alpha` across pairs. `summary()`
  prints the table, and adds a one-line note when the overlap is heterogeneous
  (the regime in which the treatment of missing individuals matters most) or
  when some pair shares very few individuals. Heterogeneity is screened on the
  coefficient of variation of `alpha` rather than its standard deviation,
  because a reweighting that is the same for every pair cancels out of the
  eigenproblem: what bites is the spread of `alpha` relative to its level. Both
  notes are display heuristics with no inferential status, and the quantity each
  is based on is printed next to it.

- `mgcca_results()` gains a `tmp_group` argument (default `"MGCCA_TMP"`) naming
  the group that holds the per-table presence masks `K`, the source of the
  overlap report. A file without that group simply gets `overlap = NULL`.

# mgcca 1.1.0

Three fixes to how individuals that are **absent from a table** are treated.
mgcca aligns the tables on the union of individual identifiers and pads the
rows a table does not hold with zeros, so that `K_j X_j = X_j`. Those padded
rows are an algebraic device, not data, and until now three per-table results
let them through. The shared components `Y`, the eigenvalues, `corsY`, the AVE
and the scaling parameters are **unchanged** — they never depended on the
padding.

- **Per-table scores are now standardised over the individuals present in that
  table.** The standard deviation used to normalise a table's scores was taken
  over the whole union, zero-padded rows included, which deflated it and made
  the returned scores and weights too large by a factor
  `sqrt((n_union - 1) / (n_present - 1))` — 12% for a table holding 24 of 30
  individuals. `scores[[j]]` and `weights[[j]]` change in scale for any table
  with missing individuals; their direction, and everything derived from the
  shared components, does not.

- **Correlation p-values use each table's own sample size.** `pval[[j]]` is the
  significance of `corsY[[j]]`, a correlation computed on the individuals
  present in table *j*, but it carried `n_union - 2` degrees of freedom. It now
  carries `n_present - 2`. The old p-values were anticonservative wherever a
  table had missing individuals.

- **Individuals absent from a table get `NA` in that table's scores.** They used
  to get exactly `0`, indistinguishable from a genuinely average measured
  individual, and `plotScores()` piled them at the origin. They are now stored
  as `NaN` in the HDF5 file and collected as `NA`; `plotScores()` leaves them
  out and reports how many. The consensus `Y` still scores every individual of
  the union — it is built by K-weighting the observed tables and is unaffected.

Code that reads `scores[[j]]` should expect `NA` for individuals outside table
*j* (e.g. `complete.cases()` before a downstream model), and absolute score
values from earlier runs are comparable to 1.1.0 values only after the scale
factor above.

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
