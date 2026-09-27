#' Meta-analysis of microbiome composition associations
#'
#' Combine study-specific Quasi-Conditional Association Test statistics using
#' fixed-effect, random-effect, or omnibus methods. This release is based on
#' the MetaPALM implementation.
#'
#' @param counts A list of nonnegative integer-valued, sample-by-taxon count
#'   matrices, one per study. Taxon column names are required. Different study
#'   columns are aligned by name; absent taxa are treated as zero counts.
#' @param X A list of numeric covariate matrices. Rows must match `counts` in
#'   number and order. Do not include an intercept; one is added internally.
#'   Adjustment covariates may differ between studies. Named study lists are
#'   aligned by study name; sample rows are never silently reordered.
#' @param X.index Column indices in each `X` matrix for the covariate(s) to test,
#'   before adding an intercept. If omitted, all columns of the first matrix
#'   are tested. Specify this explicitly when adjustment covariates are present.
#' @param cluster.id A list of subject/cluster ID vectors, one per study, or
#'   `NULL` for independent samples. Entries may be `NULL` for individual
#'   independent studies. Each supplied vector must have one nonmissing ID
#'   per sample. IDs are converted to factors with unused levels removed.
#' @param taxonomy A taxon-by-rank matrix with unique row names covering all
#'   count-matrix taxa. Columns must be named `Rank1`, `Rank2`, etc., from
#'   higher to lower rank. Extra rows are ignored. With `NULL`, perform one
#'   composition-wide test; otherwise test parents with multiple child taxa.
#' @param method One of `"FE-MetaQCAT"`, `"RE-MetaQCAT"`, `"FE-MetaQCAT*"`,
#'   `"RE-MetaQCAT*"`, `"FE-MetaQCAT-O"`, or `"RE-MetaQCAT-O"`.
#'   Unstarred methods use composition information, starred methods use zero
#'   indicators, and `-O` combines component p-values using a Cauchy combination.
#' @param min.depth Remove samples with total count less than or equal to this
#'   threshold, within each study. Default: 0 (remove zero-depth samples).
#' @param min.prev Prevalence threshold between 0 and 1. For composition tests, keep
#'   taxa with nonzero proportion at least this value and some positive counts.
#'   For zero tests, require the zero proportion strictly between `min.prev`
#'   and `1 - min.prev`. Default: 0.1. Omnibus tests may use two values, for
#'   composition and zero components respectively. This is not a maximum
#'   zero-proportion threshold.
#' @param n.perm Positive integer resampling budget, or `NULL` for asymptotic
#'   fixed-effect tests only. Random-effect methods require resampling.
#'   Adaptive resampling may stop early once ten exceedances have accumulated.
#' @param permute.strata For clustered resampling, `"within"` shuffles tested
#'   covariates within clusters and `"between"` shuffles cluster-level values
#'   between clusters. Between-cluster tested covariates must be constant
#'   within each cluster. Specify a scheme whenever clusters and `n.perm`
#'   are supplied. Choose according to the design and exchangeability assumptions.
#' @param fdr.alpha Benjamini-Hochberg FDR threshold for lineage tests.
#'   Default: 0.05. Use resampling p-values when requested, otherwise asymptotic
#'   p-values. Ineligible or failed lineage tests are not declared significant.
#' @param save.coef Whether to return descendant-taxon coefficient estimates
#'   and variances. Default: `FALSE`. The composition reference taxon is always
#'   reported when its test succeeds.
#' @param n.threads Maximum number of forked workers for lineage tests.
#'   Default: 1. Windows falls back to one worker with a warning.
#' @param seed Integer random seed for rarefaction and resampling. Default:
#'   1234. Every lineage uses the same seed, independently of worker scheduling.
#'
#' @details
#' The zero component is based on rarefied counts. After depth filtering, each
#' study is rarefied once to its own minimum sample depth, before aggregation
#' into lineages. This follows MetaPALM, not the older MetaQCAT implementation.
#' Rarefaction uses sequential hypergeometric sampling without expanding reads;
#' its distribution is unchanged, but exact seeded draws can differ from MetaPALM.
#'
#' Taxa absent from a study's count matrix are treated as zero, not as missing
#' measurements. Do not use this convention for unmeasured taxa. Studies must
#' be independent of one another; `cluster.id` only accounts for dependence
#' within each study. Failure messages are retained rather than silently
#' removing lineage results. Finite p-values do not establish statistical
#' calibration; see the vignette for design and simulation considerations.
#'
#' @return Without `taxonomy`, component fields `pvalue.comp`, `pvalue.pa`,
#'   and/or `pvalue.comb`, named `Asy` and/or `Perm`; score statistics; and
#'   requested coefficient estimates. With `taxonomy`, a list containing:
#'   * `lineage.pval`: calibration-by-lineage p-value matrix;
#'   * `lineage.qval`: corresponding BH-adjusted p-values;
#'   * `lineage.coef`: detailed component results, including optional coefficients;
#'   * `sig.lineage`: BH-significant lineage names;
#'   * `failed.lineage`: named failure reasons. Failed tests retain `NA` columns.
#' @seealso [CRC_data], [IBD_data], [Score_test_QCAT_meta],
#'   `vignette("MetaQCAT", package = "MetaQCAT")`
#' @references
#' Tang ZZ, Chen G, Alekseyenko AV, Li H (2017). A general framework for
#' association analysis of microbial communities on a taxonomic tree.
#' Bioinformatics. \doi{10.1093/bioinformatics/btw804}.
#'
#' Lee S, Teslovich TM, Boehnke M, Lin X (2013). General framework for
#' meta-analysis of rare variants in sequencing association studies.
#' American Journal of Human Genetics. \doi{10.1016/j.ajhg.2013.05.010}.
#' @examples
#' data("CRC_data")
#' names(CRC_data$counts)
#' # See the vignette for runnable CRC and clustered IBD analyses.
#' @export
#' @import MASS
#' @import brglm2
#' @import stats
#' @importFrom parallelly availableCores
MetaQCAT <- function(counts, X, X.index, cluster.id = NULL, taxonomy = NULL,
                     method = "FE-MetaQCAT", min.depth = 0, min.prev = 0.1,
                     n.perm = NULL, permute.strata = NULL, fdr.alpha = 0.05,
                     save.coef = FALSE, n.threads = 1, seed = 1234) {
  .check_arg_names(sys.call(), names(formals(MetaQCAT)))
  .check_scalar(min.depth, "min.depth")
  .check_scalar(fdr.alpha, "fdr.alpha", 0, 1)
  .check_scalar(n.threads, "n.threads", 1, .Machine$integer.max, integer = TRUE)
  .check_options(method, min.prev, n.perm, permute.strata, cluster.id, save.coef, seed)
  if (missing(X.index)) {
    if (!is.list(X) || !length(X)) stop("X must be a nonempty list.", call. = FALSE)
    X.index <- seq_len(ncol(as.matrix(X[[1L]])))
  }
  d <- .prepare_inputs(counts, X, X.index, cluster.id, taxonomy,
                        permute.strata = permute.strata, n.perm = n.perm)
  for (i in seq_along(d$counts)) {
    keep <- rowSums(d$counts[[i]]) > min.depth
    d$counts[[i]] <- d$counts[[i]][keep, , drop = FALSE]
    d$X[[i]] <- d$X[[i]][keep, , drop = FALSE]
    if (!is.null(d$cluster.id) && !is.null(d$cluster.id[[i]]))
      d$cluster.id[[i]] <- factor(d$cluster.id[[i]][keep])
  }
  keep.study <- vapply(d$counts, nrow, integer(1)) > 0L
  if (!any(keep.study)) stop("No studies remain after depth filtering.", call. = FALSE)
  if (any(!keep.study)) warning("Excluded empty studies after depth filtering: ",
    paste(names(d$counts)[!keep.study], collapse = ", "), call. = FALSE)
  counts <- d$counts[keep.study]; X <- d$X[keep.study]
  cluster.id <- if (is.null(d$cluster.id)) NULL else d$cluster.id[keep.study]
  X <- lapply(X, function(x) cbind(`(Intercept)` = 1, x))
  X.index <- d$X.index + 1L
  zeros <- if (grepl("\\*$|-O$", method))
    lapply(counts, .rarefy_counts, seed = seed) else NULL
  if (grepl("\\*$", method)) counts <- NULL
  if (is.null(taxonomy)) {
    if (!is.null(zeros)) zeros <- lapply(zeros, function(y) (y == 0) * 1)
    res <- .run_meta_test(counts, zeros, X, X.index, cluster.id, method,
                          min.prev, n.perm, permute.strata, save.coef, seed)
    if (!is.null(res$error)) warning("Test failed; returning NA: ", res$error, call. = FALSE)
    return(res)
  }
  .scan_lineages(counts, zeros, X, X.index, cluster.id, d$taxonomy, method,
                 min.prev, n.perm, permute.strata, fdr.alpha, save.coef, n.threads, seed)
}
