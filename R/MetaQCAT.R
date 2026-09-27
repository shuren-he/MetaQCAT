
########################################
#                                      #
#        Positive Part Model           #
#                                      #
########################################

.Score_test_stat <- function(Y.list, X.list, X.par.index, prev.filter = 0.1, cluster.list = NULL, meta.method = "FE-MetaQCAT") {
  # generate the summary statistics for later meta analysis
  n_studies = length(Y.list); n_avail = 0
  n.par.interest = length(X.par.index);
  study.excluded = c()
  col.index.list = Y.R.list = Y.I.list = list()
  has.cluster = !is.null(cluster.list)
  if(is.null(colnames(Y.list[[1]]))){
    taxa.names <- paste0(rep('V',ncol(Y.list[[1]])),c(1:ncol(Y.list[[1]])))
  }else{
    taxa.names <- colnames(Y.list[[1]])
  }
  taxa.index.list <- vector("list", n_studies)
  for(j in 1:n_studies){
    Y = Y.list[[j]]; X = X.list[[j]]
    if(prev.filter > 0){taxa.index = which(colMeans(Y!=0) > prev.filter)}
    else{taxa.index = which(colSums(Y) > 0)}
    if(length(taxa.index) == 0){study.excluded = append(study.excluded, j); next}
    Y = Y[,taxa.index, drop = FALSE]
    remove_rows = which(rowSums(Y) == 0)
    if(length(remove_rows) == nrow(Y)){
      study.excluded = append(study.excluded, j)
      next
    }
    if(length(remove_rows)>0){
      Y <- Y[-remove_rows, ,drop=FALSE]
      X <- X[-remove_rows, , drop=FALSE]
      if(has.cluster){
        if(!is.null(cluster.list[[j]])){
          cluster.list[[j]] <- factor(cluster.list[[j]][-remove_rows])
        }
      }
    }
    col_nonzero <- colSums(Y != 0) > 0 
    Y <- Y[, col_nonzero, drop = FALSE]
    taxa.index <- taxa.index[col_nonzero]
    if(ncol(Y) <= 1){
      study.excluded <- append(study.excluded, j)
      next
    }
    p <- ncol(X); nY <- rowSums(Y)
    n <- nrow(Y); m <- ncol(Y)
    n.beta <- m * p;
    if(is.null(colnames(Y))){ colnames(Y) <- taxa.names[taxa.index] }
    if(is.null(rownames(Y))){ rownames(Y) <- paste0(rep('S',n),c(1:n)) }
    X.reduce <- X[, -X.par.index, drop = FALSE]
    p.reduce <- p - n.par.interest
    # the index of parameters of interest
    # par.interest.index.beta <- kronecker(((0:(m - 1)) * p), rep(1, n.par.interest)) + X.par.index
    # n.par.interest.beta <- length(par.interest.index.beta)
    # beta.ini.reduce <- rep(0, p.reduce * m) # initialize the those elements of beta which we are interested in
    est.reduce.beta <- rep(0, n.beta); remove.taxa.id = c()
    Y_b <- matrix(0, nrow = nrow(Y), ncol = ncol(Y), dimnames = list(rownames(Y), colnames(Y)))
    for(i in 1:m){
      data.reduce.beta <- list(Y = Y[,i], X = X.reduce)
      glm.out.tmp <- suppressWarnings(
        try(
          glm(Y ~ X - 1,
              data   = data.reduce.beta,
              family = poisson(link = "log"),
              offset = log(nY),
              method = brglm2::brglmFit,
              type   = "AS_mean"),
          silent = TRUE
        )
      )
      if(class(glm.out.tmp)[1] != "try-error"){
        names(glm.out.tmp$coefficients) <- paste0(":", names(glm.out.tmp$coefficients))
        par.est.index.beta = (((i-1)*p + 1):(i*p))[-X.par.index]
        est.reduce.beta[par.est.index.beta] <- c(glm.out.tmp$coefficients)
        if(glm.out.tmp$converged){
          Qmat <- qr.Q(glm.out.tmp$qr); Y_b[,i] <- 0.5 * rowSums(Qmat * Qmat)
        }else{
          warning("Cannot converge for feature ", colnames(Y)[i],
                  ", remove this taxa in nul model.\n")
          remove.taxa.id = c(remove.taxa.id, i)
        }
      }else{
        warning("Cannot generate score estimate for feature ", colnames(Y)[i],
                ", remove this taxa in null model.\n")
        remove.taxa.id = c(remove.taxa.id, i)
      }
    }
    if(length(remove.taxa.id) == m){study.excluded = append(study.excluded, j); next}
    n_avail = n_avail + 1
    if(length(remove.taxa.id)>0){
      taxa.index = taxa.index[-remove.taxa.id]
      est.reduce.beta = est.reduce.beta[-(rep((remove.taxa.id-1)*p, each = p) + (1:p))]
      Y = Y[,-remove.taxa.id, drop=FALSE]
      Y_b = Y_b[,-remove.taxa.id, drop=FALSE]
      m = m - length(remove.taxa.id)
    }
    est.reduce.beta.mat = matrix(est.reduce.beta, ncol = m)
    est.reduce.beta.mat = est.reduce.beta.mat[-X.par.index, ,drop=FALSE]
    dd = crossprod(t(X.reduce), est.reduce.beta.mat)
    Y.I.list[[n_avail]] = exp(dd) * nY
    Y.R.list[[n_avail]] = Y - Y.I.list[[n_avail]] + Y_b
    taxa.index.list[[n_avail]] = taxa.index
    X.list[[j]] = X
    # idx = kronecker((col.index-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest)
    # col.zero.index.list[[n_avail]] = idx  # the index of alpha parameter of interest in this study
    # col.index.list[[n_avail]] = kronecker((taxa.index-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest)
  }
  taxa.index.list <- taxa.index.list[!vapply(taxa.index.list, is.null, logical(1))]
  taxa.include = sort(Reduce(union, taxa.index.list, init = integer(0)))
  col.index.list = lapply(taxa.index.list, function(taxa.index){
    taxa.pos =  match(taxa.index, taxa.include)
    return(kronecker((taxa.pos-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest))
  })
  n.par.interest.beta = n.par.interest * length(taxa.include)
  if(length(study.excluded) > 0){
    X.list = X.list[-study.excluded]
    if(has.cluster){cluster.list = cluster.list[-study.excluded]}
  }
  if(length(study.excluded) == n_studies){ stop("Error: No studies available for meta-analysis")}
  res <- tryCatch(
    score_test_stat_meta(X.list, cluster.list, col.index.list,
                         Y.R.list, Y.I.list, X.par.index,
                         n.par.interest.beta, meta.method),
    error = function(e) {
      stop("score_test_stat_meta failed: ", conditionMessage(e))
    }
  )
  res$Y.R.list = Y.R.list; res$Y.I.list = Y.I.list
  res$X.list = X.list; res$cluster.list = cluster.list
  res$col.index.list = col.index.list
  res$n.par.interest.beta = n.par.interest.beta
  res$beta.hat = matrix(res$beta.hat, nrow = n.par.interest, dimnames = list(NULL, taxa.names[taxa.include]))
  res$var.beta.hat = matrix(res$var.beta.hat, nrow = n.par.interest, dimnames = list(NULL, taxa.names[taxa.include]))
  return(res)
}



.Score_test_meta <- function(Y.list, X.list, X.par.index, prev.filter = 0.1, taxa.names = NULL, cluster.list = NULL,
                            n.perm = NULL, permute.strata = NULL, meta.method = "FE-MetaQCAT", save.desc = FALSE, seed = 1234){
  n_studies = length(X.list)
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if(! meta.method %in% c("FE-MetaQCAT", "RE-MetaQCAT")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # the taxa of each study should be the same
  if (length(unique(sapply(1:n_studies,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  # initialize the test statistics
  if(is.null(X.par.index)){
    stop("Error: Please provide the index(es) of covariate(s) of interest")
  }
  if (ncol(Y.list[[1]])<=1){ stop("Error: The number of taxa in each study should be greater than 1")  }
  # initialize for later use
  beta.meta.results <- tryCatch(
    .Score_test_stat(Y.list, X.list, X.par.index, prev.filter = prev.filter, cluster.list = cluster.list, meta.method = meta.method),
    error = function(e) stop(".Score_test_stat failed: ", conditionMessage(e))
  )
  score.stat.meta = beta.meta.results$score.stat.meta
  if(meta.method == "FE-MetaQCAT"){
    score.pvalue <- setNames(beta.meta.results$score.pvalue, "Asy")
    res = list(score.stat.comp = score.stat.meta, pvalue.comp = score.pvalue)
  }else{
    res = list(score.stat.comp = score.stat.meta)
  }
  if(save.desc){ res$est.beta.comp =  beta.meta.results$beta.hat; res$var.beta.comp = beta.meta.results$var.beta.hat }
  # adaptive resampling test
  if(!is.null(n.perm) && n.perm > 0){
    set.seed(seed)
    X.list = beta.meta.results$X.list
    Y.list = beta.meta.results$Y.list
    cluster.list = beta.meta.results$cluster.list
    n.par.interest.beta = beta.meta.results$n.par.interest.beta
    score.Rpvalue = resample_pvalue(X.list, cluster.list, beta.meta.results$col.index.list, beta.meta.results$Y.R.list,
                                    beta.meta.results$Y.I.list, X.par.index, score.stat.meta, n.par.interest.beta,
                                    n.perm, meta.method, permute.strata)
    if(meta.method == "FE-MetaQCAT"){ res$pvalue.comp = setNames(c(res$pvalue.comp, score.Rpvalue), c("Asy", "Perm")) }
    else{ res$pvalue.comp = setNames(c(score.Rpvalue), c("Perm")) }
  }
  return(res)
}


.Score_test_QCAT_stat <- function(Y.list, X.list, X.par.index, prev.filter = 0.1, cluster.list = NULL, meta.method = "FE-MetaQCAT") {
  # generate the summary statistics for later meta analysis
  n_studies = length(Y.list); n_avail = 0
  n.par.interest = length(X.par.index)
  study.excluded = c()
  col.index.list = list()
  has.cluster = !is.null(cluster.list)
  if(is.null(colnames(Y.list[[1]]))){
    taxa.names <- paste0(rep('V',ncol(Y.list[[1]])),c(1:ncol(Y.list[[1]])))
  }else{
    taxa.names <- colnames(Y.list[[1]])
  }
  taxa.index.list <- vector("list", n_studies)
  taxa.zero.ratio.mat <- matrix(0, nrow = n_studies, ncol = ncol(Y.list[[1]]))
  study.included = c()
  for(j in 1:n_studies){
    Y = Y.list[[j]]; X = X.list[[j]]
    taxa.zero.ratio = colMeans(Y == 0)
    if(prev.filter > 0){taxa.index = which(taxa.zero.ratio <= (1-prev.filter))}
    else{taxa.index = which(taxa.zero.ratio != 1)}
    if(length(taxa.index) <= 1){study.excluded = append(study.excluded, j); next}
    Y = Y[,taxa.index, drop = FALSE]
    remove_rows = which(rowSums(Y) == 0)
    if(length(remove_rows) == nrow(Y)){
      study.excluded = append(study.excluded, j)
      next
    }
    if(length(remove_rows)>0){
      Y <- Y[-remove_rows, ,drop=FALSE]
      X <- X[-remove_rows, , drop=FALSE]
      if(has.cluster){
        if(!is.null(cluster.list[[j]])){
          cluster.list[[j]] <- factor(cluster.list[[j]][-remove_rows])
        }
      }
    }
    col_nonzero <- colSums(Y != 0) > 0 
    Y <- Y[, col_nonzero, drop = FALSE]
    taxa.index <- taxa.index[col_nonzero]
    if(ncol(Y) <= 1){
      study.excluded <- append(study.excluded, j)
      next
    }
    n <- nrow(Y)
    if(is.null(colnames(Y))){ colnames(Y) <- taxa.names[taxa.index] }
    if(is.null(rownames(Y))){ rownames(Y) <- paste0(rep('S',n),c(1:n)) }
    n_avail = n_avail + 1
    study.included = c(study.included, j)
    taxa.index.list[[n_avail]] = taxa.index
    taxa.zero.ratio.mat[n_avail, ] = taxa.zero.ratio
    X.list[[j]] = X; Y.list[[j]] = Y
  }
  if (n_avail == 0) stop("No studies available for composition meta-analysis.")
  taxa.index.list <- taxa.index.list[!vapply(taxa.index.list, is.null, logical(1))]
  taxa.zero.ratio.mat <- taxa.zero.ratio.mat[1:n_avail, , drop = FALSE]
  taxa.reference.candidate = Reduce(intersect, taxa.index.list)
  if(length(taxa.reference.candidate) == 0){
    stop("Error: No common taxon across studies to serve as reference for compositional analysis.
        Please consider using a less stringent prevalence filter (i.e., decrease `prev.filter`).")
  }
  taxa.reference = taxa.reference.candidate[which.min(colMeans(taxa.zero.ratio.mat[, taxa.reference.candidate, drop = FALSE]))]
  if(length(study.excluded) == n_studies){ stop("Error: No studies available for meta-analysis") }
  if(length(study.excluded) > 0){
    X.list = X.list[-study.excluded]
    Y.list = Y.list[-study.excluded]
    if(has.cluster){cluster.list = cluster.list[-study.excluded]}
  }
  coef.list = vector("list", n_avail); study.excluded = c()
  for(j in 1:n_avail){
    taxa.ref =  match(taxa.reference, taxa.index.list[[j]])
    re.order <- c(c(1:ncol(Y.list[[j]]))[-taxa.ref], taxa.ref)
    Y.list[[j]] = Y.list[[j]][, re.order, drop = FALSE]
    taxa.index.list[[j]] = taxa.index.list[[j]][-taxa.ref]
    X.reduce <- X.list[[j]][, -X.par.index, drop = FALSE]
    data.reduce <- list(Y = Y.list[[j]], X = X.reduce)
    m = ncol(Y.list[[j]])
    glm.out.tmp <- suppressWarnings(
      try(brmultinom(Y ~ X - 1, data = data.reduce, type = "AS_mean", ref = m), silent = TRUE )
    )
    if(class(glm.out.tmp)[1] != "try-error"){ 
      if(glm.out.tmp$converged){
        coef.list[[j]] = t(coef(glm.out.tmp)) 
      }else{
        warning("Cannot converge for study ", study.included[j], ", remove this study in meta-analysis.\n")
        study.excluded = append(study.excluded, j)
      }
    }else{
      warning("Cannot generate score estimate for study ", study.included[j], ", remove this study in meta-analysis.\n")
      study.excluded = append(study.excluded, j)
    }
  }
  if(length(study.excluded) == n_avail){ stop("Error: No studies available for meta-analysis") }
  if(length(study.excluded) > 0){
    X.list = X.list[-study.excluded]
    Y.list = Y.list[-study.excluded]
    coef.list = coef.list[-study.excluded]
    if(has.cluster){cluster.list = cluster.list[-study.excluded]}
    taxa.index.list = taxa.index.list[-study.excluded]
  }
  taxa.include = sort(Reduce(union, taxa.index.list, init = integer(0)))
  col.index.list = lapply(taxa.index.list, function(taxa.index){
    taxa.pos = match(taxa.index, taxa.include)
    return(kronecker((taxa.pos-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest))
  })
  n.par.interest.beta = n.par.interest * length(taxa.include)
  if(length(study.excluded) == n_studies){ stop("Error: No studies available for meta-analysis") }
  res <- tryCatch(
    score_test_stat_QCAT_meta(X.list, cluster.list, col.index.list, Y.list, coef.list, X.par.index,
                              n.par.interest.beta, meta.method),
    error = function(e) {
      stop("score_test_stat_meta failed: ", conditionMessage(e))
    }
  )
  res$Y.list = Y.list
  res$X.list = X.list; res$cluster.list = cluster.list
  res$col.index.list = col.index.list
  res$coef.list = coef.list
  res$reference.taxon = taxa.names[taxa.reference]
  res$n.par.interest.beta = n.par.interest.beta
  res$beta.hat = matrix(res$beta.hat, nrow = n.par.interest, dimnames = list(NULL, taxa.names[taxa.include]))
  res$var.beta.hat = matrix(res$var.beta.hat, nrow = n.par.interest, dimnames = list(NULL, taxa.names[taxa.include]))
  return(res)
}



.Score_test_QCAT_meta <- function(Y.list, X.list, X.par.index, prev.filter = 0.1, taxa.names = NULL, cluster.list = NULL,
                            n.perm = NULL, permute.strata = NULL, meta.method = "FE-MetaQCAT", save.desc = FALSE, seed = 1234){
  n_studies = length(X.list)
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if(! meta.method %in% c("FE-MetaQCAT", "RE-MetaQCAT")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # the taxa of each study should be the same
  if (length(unique(sapply(1:n_studies,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  # initialize the test statistics
  if(is.null(X.par.index)){
    stop("Error: Please provide the index(es) of covariate(s) of interest")
  }
  if (ncol(Y.list[[1]])<=1){ stop("Error: The number of taxa in each study should be greater than 1")  }
  # initialize for later use
  beta.meta.results <- tryCatch(
    .Score_test_QCAT_stat(Y.list, X.list, X.par.index, prev.filter = prev.filter, cluster.list = cluster.list, meta.method = meta.method),
    error = function(e) stop(".Score_test_stat failed: ", conditionMessage(e))
  )
  score.stat.meta = beta.meta.results$score.stat.meta
  if(meta.method == "FE-MetaQCAT"){
    score.pvalue <- setNames(beta.meta.results$score.pvalue, "Asy")
    res = list(score.stat.comp = score.stat.meta, pvalue.comp = score.pvalue)
  }else{
    res = list(score.stat.comp = score.stat.meta)
  }
  res$reference.taxon = beta.meta.results$reference.taxon
  if(save.desc){ res$est.beta.comp = beta.meta.results$beta.hat; res$var.beta.comp = beta.meta.results$var.beta.hat }
  # adaptive resampling test
  if(!is.null(n.perm) && n.perm > 0){
    set.seed(seed)
    X.list = beta.meta.results$X.list
    Y.list = beta.meta.results$Y.list
    cluster.list = beta.meta.results$cluster.list
    n.par.interest.beta = beta.meta.results$n.par.interest.beta
    score.Rpvalue = resample_QCAT_pvalue(X.list, cluster.list, beta.meta.results$col.index.list, Y.list,
                                         beta.meta.results$coef.list, X.par.index, score.stat.meta, n.par.interest.beta,
                                         n.perm, meta.method, permute.strata)
    if(meta.method == "FE-MetaQCAT"){ res$pvalue.comp = setNames(c(res$pvalue.comp, score.Rpvalue), c("Asy", "Perm")) }
    else{ res$pvalue.comp = setNames(c(score.Rpvalue), c("Perm")) }
  }
  return(res)
}




.Score_test_zero_stat <- function(Y.bin.list, Z.list, Z.par.index, prev.filter = 0.1, cluster.list = NULL, meta.method = "FE-MetaQCAT*") {
  # generate the summary statistics for later meta analysis
  n_studies = length(Y.bin.list); n_avail = 0
  n.par.interest = length(Z.par.index);
  study.excluded = c()
  col.index.list = coef.list = list()
  has.cluster = !is.null(cluster.list)
  if(sum(Z.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if(! meta.method %in% c("FE-MetaQCAT*", "RE-MetaQCAT*")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # the taxa of each study should be the same
  if (length(unique(sapply(1:n_studies,function(j) ncol(Y.bin.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  # initialize the test statistics
  if(is.null(Z.par.index)){
    stop("Error: Please provide the index(es) of covariate(s) of interest")
  }
  if(is.null(colnames(Y.bin.list[[1]]))){
    taxa.names <- paste0(rep('V',ncol(Y.bin.list[[1]])),c(1:ncol(Y.bin.list[[1]])))
  }else{
    taxa.names <- colnames(Y.bin.list[[1]])
  }
  taxa.index.list <- lapply(seq_len(n_studies), function(j) {
    p <- colMeans(Y.bin.list[[j]])
    filter = 0
    if(prev.filter > 0){filter = prev.filter}
    which(p < 1 - filter & p > filter)
  })
  taxa.include <- sort(Reduce(union, taxa.index.list, init = integer(0)))
  n.par.interest.alpha = n.par.interest * length(taxa.include)
  taxa.names <- taxa.names[taxa.include]
  n_avail <- 0
  for (j in seq_len(n_studies)) {
    taxa.index <- taxa.index.list[[j]]
    if(length(taxa.index) == 0) { study.excluded <- c(study.excluded, j); next }
    taxa.pos <- match(taxa.index, taxa.include)
    Y.bin <- Y.bin.list[[j]][, taxa.index, drop = FALSE]
    Z  <- Z.list[[j]]
    n <- nrow(Y.bin); m <- ncol(Y.bin); p <- ncol(Z)
    if (is.null(colnames(Y.bin))) colnames(Y.bin) <- taxa.names[taxa.pos]
    if (is.null(rownames(Y.bin))) rownames(Y.bin) <- paste0("S", seq_len(n))
    Z.reduce <- Z[, -Z.par.index, drop = FALSE]
    outcome <- as.vector(t(Y.bin))
    cova.reduce <- kronecker(diag(m), Z.reduce)
    perm <- as.vector(t(matrix(1:(n*m), nrow = n, ncol = m)))  
    cova.reduce <- cova.reduce[perm, ]
    data.reduce <- data.frame(outcome = outcome, cova.reduce, row.names = NULL)
    gee.reduce <- suppressWarnings(
      try(glm(outcome ~ . - 1, family = binomial("logit"), data = data.reduce))
    )
    if(class(gee.reduce)[1] != "try-error"){ 
      if(gee.reduce$converged){
        n_avail <- n_avail + 1
        coef.list[[n_avail]] <- coef(gee.reduce)
        col.index.list[[n_avail]] <- kronecker((taxa.pos - 1) * n.par.interest,
                                              rep(1, n.par.interest)) + seq_len(n.par.interest)
      }else{
        warning("Cannot converge for study ", j, ", remove this study in meta-analysis.\n")
        study.excluded = append(study.excluded, j)
      }
    }else{
      warning("Cannot generate score estimate for study ", j, ", remove this study in meta-analysis.\n")
      study.excluded = append(study.excluded, j)
    }
    Y.bin.list[[j]] <- Y.bin
  }
  if(length(study.excluded) > 0){
    Z.list = Z.list[-study.excluded]
    Y.bin.list = Y.bin.list[-study.excluded]
    if(has.cluster){cluster.list = cluster.list[-study.excluded]}
  }
  if(length(study.excluded) == n_studies){ stop("Error: No studies available for meta-analysis") }
  res <- tryCatch(
    score_test_stat_zero_meta(Z.list, cluster.list, col.index.list, Y.bin.list,
                              coef.list, Z.par.index, n.par.interest.alpha, meta.method),
    error = function(e) { stop("score_test_stat_zero_meta failed: ", conditionMessage(e)) }
  )
  res$coef.list = coef.list
  res$col.index.list = col.index.list
  res$Z.list = Z.list
  res$Y.bin.list = Y.bin.list
  res$cluster.list = cluster.list
  res$n.par.interest.alpha = n.par.interest.alpha
  res$alpha.hat = matrix(res$alpha.hat, nrow = n.par.interest, dimnames = list(NULL, taxa.names))
  res$var.alpha.hat = matrix(res$var.alpha.hat, nrow = n.par.interest, dimnames = list(NULL, taxa.names))
  return(res)
}


.Score_test_zero_meta <- function(Y.bin.list, Z.list, Z.par.index, prev.filter = 0.1, taxa.names = NULL, cluster.list = NULL,
                                 n.perm=NULL, permute.strata = NULL, meta.method = "FE-MetaQCAT*", save.desc = FALSE, seed = 1234){
  n_studies = length(Z.list)
  for(i in 1:n_studies)
  {
    if((!is.null(taxa.names)) && is.null(colnames(Y.bin.list[[i]]))){
      colnames(Y.bin.list[[i]]) = taxa.names
    }
  }
  if(sum(Z.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Alpha part)")
  }
  if (!meta.method %in% c("FE-MetaQCAT*", "RE-MetaQCAT*")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if (length(unique(sapply(1:n_studies,function(j) ncol(Y.bin.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  if (ncol(Y.bin.list[[1]])<=1){ stop("Error: The number of taxa in each study should be greater than 1")  }
  if (is.null(Z.par.index)){ stop("Error: Please provide the index(es) of covariate(s) of interest")  }

  # initialize for later use
  alpha.meta.results <- tryCatch(
    .Score_test_zero_stat(Y.bin.list, Z.list, Z.par.index, prev.filter = prev.filter,
                         cluster.list = cluster.list, meta.method = meta.method),
    error = function(e) stop(".Score_test_zero_stat failed: ", conditionMessage(e))
  )
  score.stat.meta = alpha.meta.results$score.stat.meta
  if(meta.method == "FE-MetaQCAT*"){
    score.pvalue <- setNames(alpha.meta.results$score.pvalue, "Asy")
    res = list(score.stat.pa = score.stat.meta, pvalue.pa = score.pvalue)
  }else{
    res = list(score.stat.pa = score.stat.meta)
  }
  if(save.desc){ res$est.beta.pa = alpha.meta.results$alpha.hat; res$var.beta.pa = alpha.meta.results$var.alpha.hat }
  if(!is.null(n.perm) && n.perm > 0){
    set.seed(seed)
    Y.bin.list = alpha.meta.results$Y.bin.list
    Z.list = alpha.meta.results$Z.list
    cluster.list = alpha.meta.results$cluster.list
    n.par.interest.alpha = alpha.meta.results$n.par.interest.alpha
    score.Rpvalue = resample_zero_pvalue(Z.list, cluster.list, alpha.meta.results$col.index.list, Y.bin.list,
                                         alpha.meta.results$coef.list, Z.par.index, score.stat.meta, n.par.interest.alpha, n.perm,
                                         meta.method, permute.strata)
    if(meta.method == "FE-MetaQCAT*"){  res$pvalue.pa = setNames(c(res$pvalue.pa, score.Rpvalue), c("Asy", "Perm"))  }
    else{ res$pvalue.pa = setNames(c(score.Rpvalue), c("Perm")) }
  }
  return(res)
}

## Cauchy Combine
ACAT<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    return(1)
  }

  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }


  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}
