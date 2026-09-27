# Rarefy once per study before any lineage aggregation, as in MetaPALM.
.rarefy_counts <- function(counts, seed) {
  set.seed(seed)
  depth <- min(rowSums(counts))
  draw <- function(x) {
    out <- numeric(length(x)); remaining <- depth; population <- sum(x)
    if (length(x) > 1L) {
      for (j in seq_len(length(x) - 1L)) {
        out[j] <- rhyper(1L, x[j], population - x[j], remaining)
        remaining <- remaining - out[j]; population <- population - x[j]
      }
    }
    out[length(x)] <- remaining
    out
  }
  out <- t(vapply(seq_len(nrow(counts)), function(i) draw(counts[i, ]),
                  numeric(ncol(counts))))
  dimnames(out) <- dimnames(counts)
  out
}

.run_meta_test <- function(counts, zeros, X, X.index, cluster.id, method,
                            min.prev, n.perm, permute.strata, save.coef, seed) {
  args <- list(X.list = X, X.par.index = X.index, prev.filter = min.prev[1L],
               taxa.names = NULL, cluster.list = cluster.id, n.perm = n.perm,
               permute.strata = permute.strata, meta.method = sub("-O$", "", method),
               save.desc = save.coef, seed = seed)
  composition <- function() do.call(.Score_test_QCAT_meta,
                                    c(list(Y.list = counts), args))
  zero <- function() {
    do.call(.Score_test_zero_meta, list(
      Y.bin.list = zeros, Z.list = X, Z.par.index = X.index,
      prev.filter = min.prev[length(min.prev)], taxa.names = NULL,
      cluster.list = cluster.id, n.perm = n.perm,
      permute.strata = permute.strata,
      meta.method = if (grepl("-O$", method)) paste0(sub("-O$", "", method), "*") else method,
      save.desc = save.coef, seed = seed))
  }
  tryCatch({
    if (grepl("\\*$", method)) {
      res <- zero()
    } else if (grepl("-O$", method)) {
      res <- composition(); z <- zero()
      if (!identical(names(res$pvalue.comp), names(z$pvalue.pa)))
        stop("Component calibrations do not align.")
      combined <- setNames(vapply(seq_along(res$pvalue.comp), function(i)
        ACAT(c(res$pvalue.comp[i], z$pvalue.pa[i])), numeric(1)), names(res$pvalue.comp))
      res <- c(res, z, list(pvalue.comb = combined))
    } else {
      res <- composition()
    }
    p <- .extract_pvalues(res, method)
    if (!length(p) || any(!is.finite(p)) || any(p < 0 | p > 1))
      stop("The backend returned invalid p-values.")
    res
  }, error = function(e) {
    calibration <- if (grepl("^RE-", method)) "Perm" else
      if (is.null(n.perm)) "Asy" else c("Asy", "Perm")
    field <- if (grepl("-O$", method)) "pvalue.comb" else
      if (grepl("\\*$", method)) "pvalue.pa" else "pvalue.comp"
    res <- list(error = conditionMessage(e))
    res[[field]] <- setNames(rep(NA_real_, length(calibration)), calibration)
    res
  })
}

.extract_pvalues <- function(result, method) {
  field <- if (grepl("-O$", method)) "pvalue.comb" else
    if (grepl("\\*$", method)) "pvalue.pa" else "pvalue.comp"
  result[[field]]
}

.scan_lineages <- function(counts, zeros, X, X.index, cluster.id, taxonomy,
                            method, min.prev, n.perm, permute.strata, fdr.alpha,
                            save.coef, n.threads, seed) {
  tasks <- list()
  for (rank in rev(seq_len(ncol(taxonomy) - 1L))) {
    for (parent in unique(taxonomy[, rank])) {
      rows <- which(taxonomy[, rank] == parent)
      children <- taxonomy[rows, rank + 1L]
      if (length(unique(children)) > 1L)
        tasks[[length(tasks) + 1L]] <- list(parent = parent, rank = rank,
                                          rows = rows, children = children)
    }
  }
  keys <- vapply(tasks, function(t) t$parent, character(1))
  duplicate <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  keys[duplicate] <- vapply(tasks[duplicate], function(t)
    paste0("Rank", t$rank, ":", t$parent), character(1))
  keys <- make.unique(keys)
  aggregate <- function(y, task) {
    if (is.null(y)) return(NULL)
    lapply(y, function(m) {
      a <- t(rowsum(t(m[, task$rows, drop = FALSE]), task$children, reorder = FALSE))
      rownames(a) <- rownames(m)
      a
    })
  }
  work <- function(i) {
    task <- tasks[[i]]
    y <- aggregate(counts, task)
    z <- aggregate(zeros, task)
    if (!is.null(z)) z <- lapply(z, function(m) (m == 0) * 1)
    .run_meta_test(counts = y, zeros = z, X = X, X.index = X.index,
                   cluster.id = cluster.id, method = method, min.prev = min.prev,
                   n.perm = n.perm, permute.strata = permute.strata,
                   save.coef = save.coef, seed = seed)
  }
  if (.Platform$OS.type == "windows" && n.threads > 1L) {
    warning("Fork-based parallelism is unavailable on Windows; using one worker.",
            call. = FALSE)
    n.threads <- 1L
  }
  results <- if (n.threads == 1L || length(tasks) < 2L) lapply(seq_along(tasks), work) else
    parallel::mclapply(seq_along(tasks), work,
                       mc.cores = min(n.threads, parallelly::availableCores()))
  names(results) <- keys
  calibration <- if (grepl("^RE-", method)) "Perm" else
    if (is.null(n.perm)) "Asy" else c("Asy", "Perm")
  p <- matrix(NA_real_, length(calibration), length(tasks),
              dimnames = list(paste0(ifelse(calibration == "Asy", "Asymptotic-", "Resampling-"), method), keys))
  for (i in seq_along(results)) p[, i] <- .extract_pvalues(results[[i]], method)[calibration]
  failed <- vapply(results, function(r) if (is.null(r$error)) NA_character_ else r$error,
                   character(1))
  failed <- failed[!is.na(failed)]
  if (length(failed)) warning(length(failed),
    " lineage test(s) failed; NA columns retained. See failed.lineage for reasons.", call. = FALSE)
  q <- p
  for (i in seq_len(nrow(q))) q[i, ] <- p.adjust(p[i, ], method = "BH")
  final.q <- q[nrow(q), ]
  list(lineage.pval = p, lineage.qval = q, lineage.coef = results,
       sig.lineage = keys[!is.na(final.q) & final.q <= fdr.alpha],
       failed.lineage = failed)
}
