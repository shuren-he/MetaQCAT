test_that("all six methods work globally and by lineage", {
  d <- make_input()
  methods <- c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O",
                "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")
  for (m in methods) {
    fit <- MetaQCAT(d$counts, d$X, 1, method = m, taxonomy = d$taxonomy,
                    n.perm = 9, save.coef = TRUE)
    expect_length(fit$failed.lineage, 0L)
    expect_identical(colnames(fit$lineage.pval), c("FamilyA", "FamilyB"))
    expect_true(all(is.finite(fit$lineage.pval)))
    expect_true(all(fit$lineage.pval >= 0 & fit$lineage.pval <= 1))
    expect_identical(dim(fit$lineage.qval), dim(fit$lineage.pval))
    global <- MetaQCAT(d$counts, d$X, 1, method = m, n.perm = 9)
    expect_null(global$error)
    expect_true(all(is.finite(.extract_pvalues(global, m))))
  }
})

test_that("all six clustered resampling paths work", {
  d <- make_input()
  ids <- lapply(d$X, function(x) factor(rep(1:20, each = 2), levels = 1:25))
  for (m in c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O",
               "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")) {
    fit <- MetaQCAT(d$counts, d$X, 1, method = m, taxonomy = d$taxonomy,
                    cluster.id = ids, permute.strata = "between", n.perm = 9)
    expect_length(fit$failed.lineage, 0L)
    expect_true(all(is.finite(fit$lineage.pval)))
  }
  within.x <- lapply(d$X, function(x) {x[, 1] <- rep(0:1, 20); x})
  fit <- MetaQCAT(d$counts, within.x, 1, cluster.id = ids,
                  permute.strata = "within", n.perm = 9,
                  method = "RE-MetaQCAT-O", taxonomy = d$taxonomy)
  expect_length(fit$failed.lineage, 0L)
  expect_true(all(is.finite(fit$lineage.pval)))
})

test_that("alignment and single-lineage dimensions are preserved", {
  d <- make_input()
  ref <- MetaQCAT(d$counts, d$X, 1, taxonomy = d$taxonomy)
  y <- d$counts; y[[2]] <- y[[2]][, 6:1]
  fit <- MetaQCAT(y, rev(d$X), 1, taxonomy = d$taxonomy[6:1, ])
  expect_equal(fit$lineage.pval, ref$lineage.pval)
  tax <- d$taxonomy[1:3, , drop = FALSE]
  y <- lapply(d$counts, function(y) y[, 1:3])
  fit <- MetaQCAT(y, d$X, 1, taxonomy = tax)
  expect_identical(dim(fit$lineage.pval), c(1L, 1L))
  expect_identical(colnames(fit$lineage.pval), "FamilyA")
  empty.tax <- cbind(Rank1 = LETTERS[1:6], Rank2 = LETTERS[1:6])
  rownames(empty.tax) <- LETTERS[1:6]
  empty <- MetaQCAT(d$counts, d$X, 1, taxonomy = empty.tax)
  expect_identical(dim(empty$lineage.pval), c(1L, 0L))
})

test_that("failed lineages remain as NA with diagnostics", {
  d <- make_input()
  y <- lapply(d$counts, function(y) {y[, 1:3] <- 0; y})
  expect_warning(fit <- MetaQCAT(y, d$X, 1, taxonomy = d$taxonomy), "NA columns retained")
  expect_identical(colnames(fit$lineage.pval), c("FamilyA", "FamilyB"))
  expect_true(all(is.na(fit$lineage.pval[, "FamilyA"])))
  expect_true(all(is.finite(fit$lineage.pval[, "FamilyB"])))
  expect_identical(names(fit$failed.lineage), "FamilyA")
  expect_false("FamilyA" %in% fit$sig.lineage)
})

test_that("depth filtering excludes only truly empty studies", {
  d <- make_input()
  y <- d$counts; y[[1]][1:20, ] <- 0
  expect_silent(fit <- MetaQCAT(y, d$X, 1))
  expect_null(fit$error)
  y[[1]][, ] <- 0
  expect_warning(fit <- MetaQCAT(y, d$X, 1), "Excluded empty studies")
  expect_null(fit$error)
  y[[2]][, ] <- 0
  expect_error(MetaQCAT(y, d$X, 1), "No studies remain")
})

test_that("rarefaction preserves dimensions, bounds, and depth", {
  d <- make_input()
  y <- d$counts[[1]]
  r <- .rarefy_counts(y, 77)
  expect_identical(dimnames(r), dimnames(y))
  expect_true(all(rowSums(r) == min(rowSums(y))))
  expect_true(all(r >= 0 & r <= y))
  expect_identical(r, .rarefy_counts(y, 77))
  one <- y[1, , drop = FALSE]
  expect_identical(dim(.rarefy_counts(one, 77)), dim(one))
})

test_that("advanced interfaces have explicit intercept and index rules", {
  d <- make_input()
  x <- lapply(d$X, function(x) cbind(1, x))
  for (f in list(Score_test_meta, Score_test_QCAT_meta)) {
    fit <- f(d$counts, x, 2, n.perm = 9, save.coef = TRUE)
    expect_true(all(is.finite(fit$pvalue.comp)))
  }
  z <- lapply(d$counts, function(y) (y == 0) * 1)
  fit <- Score_test_zero_meta(z, x, 2, n.perm = 9, save.coef = TRUE)
  expect_true(all(is.finite(fit$pvalue.pa)))
  for (f in list(Score_test_stat, Score_test_QCAT_stat))
    expect_true(is.finite(f(d$counts, x, 2)$score.stat.meta))
  expect_true(is.finite(Score_test_zero_stat(z, x, 2)$score.stat.meta))
  expect_error(Score_test_QCAT_meta(d$counts, x, 1), "intercept")
})

test_that("results are reproducible across worker schedules", {
  d <- make_input()
  a <- MetaQCAT(d$counts, d$X, 1, taxonomy = d$taxonomy, n.perm = 9, seed = 7,
                method = "RE-MetaQCAT-O")
  b <- MetaQCAT(d$counts, d$X, 1, taxonomy = d$taxonomy, n.perm = 9, seed = 7,
                method = "RE-MetaQCAT-O", n.threads = 2)
  expect_identical(a$lineage.pval, b$lineage.pval)
})

test_that("packaged data keep their distinct designs", {
  data("CRC_data", package = "MetaQCAT", envir = environment())
  data("IBD_data", package = "MetaQCAT", envir = environment())
  expect_false("design" %in% names(CRC_data))
  expect_null(CRC_data$subject_ids)
  expect_true("design" %in% names(IBD_data))
  expect_identical(sum(vapply(CRC_data$counts, nrow, integer(1))), 573L)
  expect_identical(sum(vapply(IBD_data$counts, nrow, integer(1))), 3472L)
})

test_that("adaptive resampling completes a non-round budget when needed", {
  d <- make_input()
  x <- lapply(d$X, function(x) cbind(1, x))
  q <- Score_test_QCAT_stat(d$counts, x, 2)
  p <- resample_QCAT_pvalue(q$X.list, q$cluster.list, q$col.index.list,
                            q$Y.list, q$coef.list, 2, 1e300,
                            q$n.par.interest.beta, 123L, "FE-MetaQCAT", NULL)
  expect_equal(p, 1 / 124)
  z <- lapply(d$counts, function(y) (y == 0) * 1)
  a <- Score_test_zero_stat(z, x, 2)
  p <- resample_zero_pvalue(a$Z.list, a$cluster.list, a$col.index.list,
                            a$Y.bin.list, a$coef.list, 2, 1e300,
                            a$n.par.interest.alpha, 123L, "FE-MetaQCAT*", NULL)
  expect_equal(p, 1 / 124)
  b <- Score_test_stat(d$counts, x, 2)
  p <- resample_pvalue(b$X.list, b$cluster.list, b$col.index.list,
                       b$Y.R.list, b$Y.I.list, 2, 1e300,
                       b$n.par.interest.beta, 123L, "FE-MetaQCAT", NULL)
  expect_equal(p, 1 / 124)
  re <- Score_test_stat(d$counts, x, 2, method = "RE-MetaQCAT")
  expect_gt(re$score.stat.meta, b$score.stat.meta)
})

test_that("multiple tested covariates retain correct coefficient dimensions", {
  d <- make_input()
  x <- lapply(d$X, function(x) {
    cbind(x, trend = seq_len(nrow(x)) %% 3)
  })
  for (m in c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O")) {
    fit <- MetaQCAT(d$counts, x, 1:2, taxonomy = d$taxonomy,
                    method = m, n.perm = 9, save.coef = TRUE)
    expect_length(fit$failed.lineage, 0L)
    expect_true(all(is.finite(fit$lineage.pval)))
    r <- fit$lineage.coef[[1L]]
    if (!is.null(r$est.beta.comp)) expect_identical(nrow(r$est.beta.comp), 2L)
    if (!is.null(r$est.beta.pa)) expect_identical(nrow(r$est.beta.pa), 2L)
  }
})
