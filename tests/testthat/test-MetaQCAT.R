library(testthat)
library(PALMcov)

data("CRC_data", package = "PALMcov")
CRC_abd <- CRC_data$CRC_abd
CRC_meta <- CRC_data$CRC_meta

rel.abd <- list()
covariate.interest <- list()
cluster.id.list <- list()
for(d in unique(CRC_meta$Study)){
  rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
  disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
  names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
  covariate.interest[[d]] <- matrix(c(disease), ncol = 1)
  covariate.interest[[d]] <- cbind(1, covariate.interest[[d]], rbinom(length(disease), 1, 0.5))
  covariate.interest[[d]] <- cbind(1, covariate.interest[[d]])
  cluster.id.list[[d]] <- sample(1:10, size = length(disease), replace = T)
}

Score.test.stat <- function(Y.list, X.list, X.par.index, prev.filter = 0.1, cluster.id.list = NULL, meta.method = "FE-MetaQCAT") {
  # generate the summary statistics for later meta analysis
  total.num = length(Y.list); ava.cnt = 0
  n.par.interest = length(X.par.index); 
  n.par.interest.beta = n.par.interest * ncol(Y.list[[1]]);
  remove.study = c()
  col.index.list = Y.R.list = Y.I.list = list()
  flag = !is.null(cluster.id.list)
  for(j in 1:total.num){
    Y = Y.list[[j]]; X = X.list[[j]]
    col.index = which(colMeans(Y!=0) > prev.filter)
    Y = Y[,col.index, drop = FALSE]
    p <- ncol(X); nY <- rowSums(Y)
    n <- nrow(Y); m <- ncol(Y)
    n.beta <- m * p; 
    if(is.null(colnames(Y))){
      colnames(Y) <- paste0(rep('V',m),c(1:m))
    }
    if(is.null(rownames(Y))){
      rownames(Y) <- paste0(rep('S',n),c(1:n))
    }
    taxa.names = colnames(Y);
    if (sum(X.par.index == 1)) {
      stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
    }
    if (is.null(X.par.index) || n == 0) {
      stop("Error: Not proper data")
    }else {
      X.reduce <- X[, -X.par.index, drop = FALSE]
      p.reduce <- p - n.par.interest
      # the index of parameters of interest
      # par.interest.index.beta <- kronecker(((0:(m - 1)) * p), rep(1, n.par.interest)) + X.par.index
      # n.par.interest.beta <- length(par.interest.index.beta)
      # beta.ini.reduce <- rep(0, p.reduce * m) # initialize the those elements of beta which we are interested in
      est.reduce.beta <- rep(0, n.beta); remove.tax.id = c()
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
            warning("Cannot converge for feature ", tax.names[i],
                    ", remove this taxa in nul model.\n")
            remove.tax.id = c(remove.tax.id, i)
          }
        }else{
          warning("Cannot generate score estimate for feature ", tax.names[i],
                  ", remove this taxa in null model.\n")
          remove.tax.id = c(remove.tax.id, i)
        }
      }
    }
    if(length(remove.tax.id) == p){ remove.study = append(remove.study, j); next}
    ava.cnt = ava.cnt + 1
    if(length(remove.tax.id)>0){
      col.index= col.index[remove.tax.id]
      est.reduce.beta = est.reduce.beta[-(rep((remove.tax.id-1)*p, each = p) + (1:p))]
      Y = Y[,-remove.tax.id, drop=FALSE]
      Y_b = Y_b[,-remove.tax.id, drop=FALSE]
      m = m - length(remove.tax.id)
    }
    est.reduce.beta.mat = matrix(est.reduce.beta, ncol = m)
    est.reduce.beta.mat = est.reduce.beta.mat[-X.par.index, ,drop=FALSE]
    dd = crossprod(t(X.reduce), est.reduce.beta.mat)
    Y.I.list[[ava.cnt]] = exp(dd) * nY
    Y.R.list[[ava.cnt]] = Y - Y.I.list[[ava.cnt]] + Y_b

    # idx = kronecker((col.index-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest)
    # col.zero.index.list[[ava.cnt]] = idx  # the index of alpha parameter of interest in this study
    col.index.list[[ava.cnt]] = kronecker((col.index-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest)
  }
  if(length(remove.study) > 0){
    X.list = X.list[-remove.study]
    if(flag){cluster.id.list = cluster.id.list[-remove.study]}
  }
  if(length(remove.study) == total.num){ stop("Error: Not proper data") }
  res <- try(score_test_stat_meta(X.list, cluster.id.list, col.index.list,
            Y.R.list, Y.I.list, X.par.index, n.par.interest.beta, meta.method))
  res$Y.R.list = Y.R.list; res$Y.I.list = Y.I.list
  res$col.index.list = col.index.list
  return(res)
}

Score.test.meta <- function(Y.list, X.list, X.par.index, tax.names = NULL, cluster.id.list = NULL, permute.strata = NULL, 
                            seed=11, resample=FALSE, n.replicates=NULL, meta.method = "FE-MetaQCAT"){
  set.seed(seed)
  stu.num = length(X.list)
  n.par.interest = length(X.par.index)
  m = ncol(Y.list[[1]])
  # p = ncol(X.list[[1]])
  n.par.interest.beta = m*length(X.par.index)
  cls.flg = !is.null(cluster.id.list)
  remove.study = NULL
  for(i in 1:stu.num){
    remove.rows.idx = which(rowSums(Y.list[[i]]) == 0)
    if(length(remove.rows.idx) == nrow(Y.list[[i]])){
      remove.study = append(remove.study, i)
      next
    }
    if(length(remove.rows.idx)>0){
      Y.list[[i]] <- Y.list[[i]][-remove.rows.idx, ,drop=FALSE]
      X.list[[i]] <- X.list[[i]][-remove.rows.idx, , drop=FALSE]
      if(cls.flg){
        if(!is.null(cluster.id.list[[i]])){
          cluster.id.list[[i]] <- factor(cluster.id.list[[i]][-remove.rows.idx])
        }
      }
    }
  }
  if(!is.null(remove.study)){
    if(length(remove.study) == stu.num){
      stop("Error: Not proper data")
    }else{
      Y.list <- Y.list[-remove.study]
      X.list <- X.list[-remove.study]
      if(cls.flg){
        cluster.id.list <- cluster.id.list[-remove.study]
      }
      stu.num = length(X.list)
    }
  }

  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if(! meta.method %in% c("FE-MetaQCAT", "RE-MetaQCAT")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # the taxa of each study should be the same
  if (length(unique(sapply(1:stu.num,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  # initialize the test statistics
  if(is.null(X.par.index)){
    stop("Error: Please provide the index(es) of covariate(s) of interest")
  }

  # initialize for later use
  beta.meta.results = try(Score.test.stat(Y.list, X.list, X.par.index, prev.filter = 0.1, cluster.id.list = cluster.id.list, meta.method = meta.method))
  if(class(beta.meta.results)[1] != "try-error"){
    beta.hat = matrix(beta.meta.results$beta.hat, nrow = n.par.interest)
    var.beta.hat = matrix(beta.meta.results$var.beta.hat, nrow = n.par.interest)
    colnames(beta.hat) <- tax.names; colnames(var.beta.hat) <- tax.names;
    score.stat.meta = beta.meta.results$score.stat.meta
  }else{ score.pvalue = NA; score.stat.meta = NA }
  if(meta.method == "FE-MetaQCAT"){
    score.pvalue = beta.meta.results$score.pvalue
    res = list(beta.hat = beta.hat, var.beta.hat = var.beta.hat, score.stat.meta = score.stat.meta, score.pvalue = score.pvalue)
  }else{
    res = list(beta.hat = beta.hat, var.beta.hat = var.beta.hat, score.stat.meta = score.stat.meta)
  }
  # if resample = TRUE then will apply permutation method to get permuted p-value
  # adaptive resampling test
  if(resample){
    score.Rpvalue = resample_pvalue(X.list, cluster.id.list, beta.meta.results$col.index.list, beta.meta.results$Y.R.list,
                                    beta.meta.results$Y.I.list, X.par.index, score.stat.meta, n.par.interest.beta,
                                    n.replicates, meta.method, permute.strata)
    res$score.Rpvalue = score.Rpvalue
  }
  return(res)
}

# test_res <- Score.test.stat(rel.abd, covariate.interest, c(2,3), prev.filter = 0.1, cluster.id.list = class.id, meta.method = "FE-MetaQCAT")
test_res <-  Score.test.meta(rel.abd, covariate.interest, c(2), tax.names = NULL, cluster.id.list = cluster.id.list, permute.strata = "within",  
                              seed=1234, resample = TRUE, n.replicates = 500, meta.method = "FE-MetaQCAT")

save(test_res, file = "/Users/shurenhe/Downloads/code/Meta_test/test_res_1.RData")