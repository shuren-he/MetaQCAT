make_input <- function() {
  counts <- X <- vector("list", 2L)
  names(counts) <- names(X) <- c("Study1", "Study2")
  for (i in 1:2) {
    set.seed(100 + i)
    y <- t(rmultinom(40, 120, c(.2, .15, .15, .2, .15, .15)))
    for (j in seq_len(ncol(y))) y[sample(40, 10), j] <- 0
    dimnames(y) <- list(paste0("S", i, "_", 1:40), LETTERS[1:6])
    counts[[i]] <- y
    X[[i]] <- matrix(rep(0:1, each = 20), ncol = 1,
                      dimnames = list(rownames(y), "disease"))
  }
  taxonomy <- cbind(Rank1 = rep(c("FamilyA", "FamilyB"), each = 3),
                    Rank2 = LETTERS[1:6])
  rownames(taxonomy) <- LETTERS[1:6]
  list(counts = counts, X = X, taxonomy = taxonomy)
}
