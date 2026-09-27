.score_input <- function(counts, X, X.index, min.prev, cluster.id, method,
                          n.perm = NULL, permute.strata = NULL,
                          save.coef = FALSE, seed = 1234, binary = FALSE,
                          stat.only = FALSE) {
  .check_options(method, min.prev,
                 if (stat.only && grepl("^RE-", method)) 1L else n.perm,
                 permute.strata, if (stat.only) NULL else cluster.id,
                 save.coef, seed, omnibus = FALSE)
  .prepare_inputs(counts, X, X.index, cluster.id = cluster.id, binary = binary,
                   intercept = TRUE, permute.strata = permute.strata, n.perm = n.perm)
}

#' Low-level composition summary statistics
#'
#' Compute study and meta-analysis score summaries. These advanced interfaces
#' do not add intercepts or rarefy counts. Prefer [MetaQCAT()] for routine use.
#' `Score_test_stat()` uses the retained Poisson score implementation;
#' `Score_test_QCAT_stat()` uses the multinomial composition implementation
#' that underlies `MetaQCAT()`.
#'
#' @inheritParams MetaQCAT
#' @param X A list of numeric design matrices including an intercept of ones
#'   in column 1. Sample rows align with `counts`.
#' @param X.index Tested design columns, including the intercept's position
#'   in the indexing. Do not test column 1. For disease in column 2 use 2.
#' @param method `"FE-MetaQCAT"` or `"RE-MetaQCAT"`. Random-effect summary
#'   statistics can be computed without a permutation budget.
#' @param min.prev Minimum nonzero proportion. The Poisson variant retains
#'   taxa strictly above this threshold; the multinomial variant retains taxa
#'   at or above it. Entirely zero taxa are always excluded.
#' @return Backend score summaries, coefficients, covariance matrices, and
#'   filtered inputs used for resampling. The multinomial variant also returns
#'   `reference.taxon`. These low-level fields are intended for method development.
#' @export
Score_test_stat <- function(counts, X, X.index, min.prev = 0.1,
                            cluster.id = NULL, method = "FE-MetaQCAT") {
  .check_arg_names(sys.call(), names(formals(Score_test_stat)))
  d <- .score_input(counts, X, X.index, min.prev, cluster.id, method, stat.only = TRUE)
  .Score_test_stat(Y.list = d$counts, X.list = d$X, X.par.index = d$X.index,
                   prev.filter = min.prev, cluster.list = d$cluster.id, meta.method = method)
}

#' @rdname Score_test_stat
#' @export
Score_test_QCAT_stat <- function(counts, X, X.index, min.prev = 0.1,
                                 cluster.id = NULL, method = "FE-MetaQCAT") {
  .check_arg_names(sys.call(), names(formals(Score_test_QCAT_stat)))
  d <- .score_input(counts, X, X.index, min.prev, cluster.id, method, stat.only = TRUE)
  .Score_test_QCAT_stat(Y.list = d$counts, X.list = d$X, X.par.index = d$X.index,
                        prev.filter = min.prev, cluster.list = d$cluster.id, meta.method = method)
}

#' Low-level composition meta-analysis
#'
#' Compute asymptotic and/or resampling p-values without adding an intercept
#' or carrying out taxonomy scans. `Score_test_meta()` uses the retained
#' Poisson score variant; `Score_test_QCAT_meta()` uses the multinomial
#' composition implementation underlying [MetaQCAT()].
#' @inheritParams Score_test_stat
#' @inheritParams MetaQCAT
#' @param method `"FE-MetaQCAT"` or `"RE-MetaQCAT"`; random-effect p-values
#'   require a positive `n.perm`.
#' @return A list with `score.stat.comp`, `pvalue.comp` (named `Asy` and/or
#'   `Perm`), and optional coefficient/variance estimates. The QCAT variant
#'   also reports the reference taxon.
#' @export
Score_test_meta <- function(counts, X, X.index, min.prev = 0.1,
                            cluster.id = NULL, n.perm = NULL,
                            permute.strata = NULL, method = "FE-MetaQCAT",
                            save.coef = FALSE, seed = 1234) {
  .check_arg_names(sys.call(), names(formals(Score_test_meta)))
  d <- .score_input(counts, X, X.index, min.prev, cluster.id, method, n.perm,
                     permute.strata, save.coef, seed)
  .Score_test_meta(Y.list = d$counts, X.list = d$X, X.par.index = d$X.index,
                   prev.filter = min.prev, cluster.list = d$cluster.id,
                   n.perm = n.perm, permute.strata = permute.strata,
                   meta.method = method, save.desc = save.coef, seed = seed)
}

#' @rdname Score_test_meta
#' @export
Score_test_QCAT_meta <- function(counts, X, X.index, min.prev = 0.1,
                                 cluster.id = NULL, n.perm = NULL,
                                 permute.strata = NULL, method = "FE-MetaQCAT",
                                 save.coef = FALSE, seed = 1234) {
  .check_arg_names(sys.call(), names(formals(Score_test_QCAT_meta)))
  d <- .score_input(counts, X, X.index, min.prev, cluster.id, method, n.perm,
                     permute.strata, save.coef, seed)
  .Score_test_QCAT_meta(Y.list = d$counts, X.list = d$X, X.par.index = d$X.index,
                        prev.filter = min.prev, cluster.list = d$cluster.id,
                        n.perm = n.perm, permute.strata = permute.strata,
                        meta.method = method, save.desc = save.coef, seed = seed)
}

#' Low-level zero-pattern summary statistics
#'
#' Compute summaries from binary zero indicators (1 = zero, 0 = present).
#' This function does not rarefy counts or add an intercept.
#' @inheritParams Score_test_stat
#' @param zeros A list of sample-by-taxon binary matrices, with 1 denoting
#'   an observed zero and 0 denoting presence. Taxon column names are required.
#' @param min.prev Keep taxa with zero proportions strictly between this
#'   threshold and `1 - min.prev`, within each study. Default: 0.1.
#' @param method `"FE-MetaQCAT*"` or `"RE-MetaQCAT*"`.
#' @return Backend score summaries, coefficients, covariance matrices, and
#'   filtered inputs used for zero-pattern resampling.
#' @export
Score_test_zero_stat <- function(zeros, X, X.index, min.prev = 0.1,
                                 cluster.id = NULL, method = "FE-MetaQCAT*") {
  .check_arg_names(sys.call(), names(formals(Score_test_zero_stat)))
  d <- .score_input(zeros, X, X.index, min.prev, cluster.id, method,
                     binary = TRUE, stat.only = TRUE)
  .Score_test_zero_stat(Y.bin.list = d$counts, Z.list = d$X, Z.par.index = d$X.index,
                        prev.filter = min.prev, cluster.list = d$cluster.id, meta.method = method)
}

#' Low-level zero-pattern meta-analysis
#'
#' Compute asymptotic and/or resampling p-values from binary zero indicators.
#' This is the C++-backed zero component used by [MetaQCAT()].
#' @inheritParams Score_test_zero_stat
#' @inheritParams MetaQCAT
#' @param method `"FE-MetaQCAT*"` or `"RE-MetaQCAT*"`; random-effect p-values
#'   require a positive `n.perm`.
#' @return A list with `score.stat.pa`, `pvalue.pa` (named `Asy` and/or `Perm`),
#'   and optional coefficient/variance estimates. Coefficients refer to the
#'   probability of a zero, not to the probability of presence.
#' @export
Score_test_zero_meta <- function(zeros, X, X.index, min.prev = 0.1,
                                 cluster.id = NULL, n.perm = NULL,
                                 permute.strata = NULL, method = "FE-MetaQCAT*",
                                 save.coef = FALSE, seed = 1234) {
  .check_arg_names(sys.call(), names(formals(Score_test_zero_meta)))
  d <- .score_input(zeros, X, X.index, min.prev, cluster.id, method, n.perm,
                     permute.strata, save.coef, seed, binary = TRUE)
  .Score_test_zero_meta(Y.bin.list = d$counts, Z.list = d$X, Z.par.index = d$X.index,
                        prev.filter = min.prev, cluster.list = d$cluster.id,
                        n.perm = n.perm, permute.strata = permute.strata,
                        meta.method = method, save.desc = save.coef, seed = seed)
}
