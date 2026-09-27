test_that("the release uses consistent names without old aliases", {
  expect_identical(names(formals(MetaQCAT)), c("counts", "X", "X.index",
    "cluster.id", "taxonomy", "method", "min.depth", "min.prev", "n.perm",
    "permute.strata", "fdr.alpha", "save.coef", "n.threads", "seed"))
  expect_identical(formals(MetaQCAT)$min.prev, 0.1)
  expect_identical(formals(MetaQCAT)$seed, 1234)
  for (old in c("OTU", "Tax", "cluster", "meta.method", "prev.filter",
                "save.desc", "nthreads")) {
    a <- list(counts = list(), X = list(), X.index = 1)
    a[[old]] <- TRUE
    expect_error(do.call(MetaQCAT, a), "unused argument")
  }
  for (f in list(Score_test_stat, Score_test_meta, Score_test_QCAT_stat,
                 Score_test_QCAT_meta, Score_test_zero_stat, Score_test_zero_meta)) {
    expect_true(all(c("X", "X.index", "min.prev", "cluster.id", "method") %in%
                      names(formals(f))))
    expect_false(any(c("Y.list", "X.list", "Z.list", "prev.filter",
                        "cluster.list", "meta.method") %in% names(formals(f))))
  }
})

test_that("internal component arguments are forwarded by name", {
  d <- make_input()
  local_mocked_bindings(.Score_test_QCAT_meta = function(Y.list, X.list,
    X.par.index, prev.filter, taxa.names, cluster.list, n.perm,
    permute.strata, meta.method, save.desc, seed) {
      expect_identical(prev.filter, .2)
      expect_null(taxa.names)
      expect_identical(meta.method, "FE-MetaQCAT")
      expect_identical(save.desc, TRUE)
      expect_identical(seed, 71)
      list(pvalue.comp = c(Asy = .5))
    })
  r <- MetaQCAT(d$counts, d$X, 1, min.prev = .2, save.coef = TRUE, seed = 71)
  expect_identical(r$pvalue.comp, c(Asy = .5))
})

test_that("input and cluster errors are explicit", {
  d <- make_input()
  ids <- lapply(d$X, function(x) factor(rep(1:20, each = 2)))
  expect_error(MetaQCAT(d$counts, d$X, 1, cluster.id = ids, n.perm = 9),
               "Choose permute.strata")
  expect_error(MetaQCAT(d$counts, d$X, 1, cluster.id = ids, n.perm = 9,
                        permute.strata = "within"), "constant in every cluster")
  bad.ids <- lapply(d$X, function(x) factor(rep(1:20, 2)))
  expect_error(MetaQCAT(d$counts, d$X, 1, cluster.id = bad.ids, n.perm = 9,
                        permute.strata = "between"), "constant within")
  bad.x <- d$X; bad.x[[1]] <- bad.x[[1]][40:1, , drop = FALSE]
  expect_error(MetaQCAT(d$counts, bad.x, 1), "row names")
  bad.x <- lapply(d$X, function(x) cbind(1, x))
  expect_error(MetaQCAT(d$counts, bad.x, 2), "rank deficient")
  bad.y <- d$counts; bad.y[[1]][1, 1] <- NA
  expect_error(MetaQCAT(bad.y, d$X, 1), "Invalid")
  expect_error(MetaQCAT(d$counts, d$X, 1, method = "RE-MetaQCAT"), "positive n.perm")
  expect_error(MetaQCAT(d$counts, d$X, 1, n.perm = 0), "n.perm")
  expect_error(MetaQCAT(d$counts, d$X, 1, min.prev = c(.1, .2)), "two values")
})
