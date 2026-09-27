# MetaQCAT 0.0.3.0000

## MetaPALM-based release

- Renamed the package, shared-library registration, generated Rcpp bindings,
  project file, and references to MetaQCAT. MetaPALM's original directory is
  retained outside this release as the reference implementation.
- Unified exported argument names: counts, X, X.index, cluster.id, taxonomy,
  method, min.prev, save.coef, and n.threads. Zero-indicator APIs use zeros.
  Old aliases and partial argument names are not accepted. Six low-level
  Score_test functions remain available with explicit intercept conventions.
- Added a public seed argument for rarefaction and resampling.

## Correctness and diagnostics

- Internal component dispatch uses named arguments, preventing accidental
  binding of taxon names or cluster IDs to prevalence parameters.
- Restored random-effect method labels in the retained Poisson backend.
- Aligned taxa and named study lists before analysis; checked sample alignment,
  covariates, taxonomy, tested indices, and cluster IDs before calling C++.
- Corrected empty-study detection after depth filtering and no-available-study
  handling in the composition component.
- Kept failed lineage columns as NA with explicit failure reasons rather than
  silently dropping them. Added lineage.pval and lineage.qval matrices and
  retained component results and reference taxa.
- Required an explicit clustered permutation scheme for FE as well as RE;
  checked exposure constancy for between-cluster and nontrivial within-cluster
  permutations. Removed unused factor levels before entering C++.
- Capped adaptive batches at the requested budget so remaining draws are not
  omitted. Rejected non-finite permutation statistics and degenerate RE
  covariance rather than returning a misleading numeric p-value.
- Replaced read-expansion rarefaction with sequential hypergeometric sampling:
  the sampling distribution is unchanged, but exact seeded draws may differ.
- Removed machine-specific compiler include paths and a deprecated finite check.

## Documentation and validation

- Added a styled HTML vignette with runnable CRC and clustered IBD examples,
  all six methods, result interpretation, failure diagnostics, and migration.
- Integrated current CRC/IBD dataset documentation into R/data_doc.R.
- Added automated regression tests for interfaces, all six independent and
  correlated resampling paths, dimensions, alignment, failures, and reproducibility.

## Scope of validation

The zero resampling interface in this MetaPALM-based source is connected; the
missing-interface issue in the older MetaQCAT source does not apply here.
Execution tests are not proofs of type-I error control. Cross-study subject
overlap, covariate-dependent exchangeability, sparse-data calibration, and
adaptive-resampling precision still require scientific sensitivity assessment.
