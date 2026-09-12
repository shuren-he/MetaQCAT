F.test <- function(x) {
  # Fisher's p-value combination
  x.stat <- -2 * sum(log(x))
  return(1 - pchisq(x.stat, df = 2 * length(x)))
}

diag2 <- function(x){
  # transform the numeric into diag matrix
  if(length(x)>1){  return(diag(x)) }
  else{ return(as.matrix(x)) }
}

########################################
#                                      #
#        Positive Part Model           #
#                                      #
########################################

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
    if(prev.filter > 0){col.index = which(colMeans(Y!=0) > prev.filter)}
    else{col.index = which(colSums(Y) > 0)}
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

Score.test.meta <- function(Y.list, X.list, X.par.index, prev.filter = 0.1, tax.names = NULL, cluster.id.list = NULL, permute.strata = NULL, 
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
  beta.meta.results = try(Score.test.stat(Y.list, X.list, X.par.index, prev.filter = prev.filter, cluster.id.list = cluster.id.list, meta.method = meta.method))
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


########################################
#                                      #
#          Zero Part Model             #
#                                      #
########################################

Pi.alpha<-function(m, p, alpha, X.i){
  # calculate the exponential of alpha times X
  Pi.out = rep(NA,m)

  for(j in 1:m){

    tmp = exp(crossprod(alpha[((j-1)*p+1):(j*p)], X.i))
    if(is.infinite(tmp)){
      Pi.out[j] = 1
    }else{
      Pi.out[j] = tmp/(tmp + 1)
    }

  }

  # no need for base in GEE method
  return (Pi.out)
}

fun.score.i.alpha <- function(alpha, data, save.list=FALSE){
  # return  the negative log likelihood of estimated pi for each taxa
  Y = data$Y; Z = data$Z;

  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)

  vA.list = list()
  Vinv.list = list()
  VY.list = list()

  n.alpha = m*p
  # check the dimension of alpha
  if(length(alpha)!=n.alpha){

    warning("Dim of initial alpha does not match the dim of covariates")

  }else{

    Score.alpha.i = matrix(0, n, n.alpha)
    nY = rowSums(Y)

    for(i in 1:n){

      Pi.i = Pi.alpha(m, p, alpha, Z[i,])
      vA.tmp = Pi.i*(1-Pi.i)
      A.i = diag2(vA.tmp) # transform into diagonal matrices
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure

      tmp.V.i = ginv(V.i)
      tmp.VY = crossprod(t(tmp.V.i), (Y[i,] - Pi.i))
      Score.alpha.i[i,] = crossprod(t(t.D.i), tmp.VY) # score value

      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }


  }
  # if save.list = TRUE, this will be used for resampling test
  if(save.list){

    return ( list(Score.alpha=Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list=VY.list) )

  }else{

    return (Score.alpha.i)
  }

}

#fun.hessian.alpha(est.reduce.alpha, data.alpha)
fun.hessian.alpha <- function(alpha, data){

  Y = data$Y; Z = data$Z

  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  n.alpha = m*p

  if(length(alpha)!=n.alpha){
    print("Waring: dim of alpha is not the same as alpha\n")

  }else{

    Hessian.alpha = matrix(0, nrow=n.alpha, ncol=n.alpha)
    nY = rowSums(Y)


    for(i in 1:n){

      Pi.i = Pi.alpha(m, p, alpha, Z[i,])
      tmp = Pi.i*(1-Pi.i)
      A.i = diag2(tmp)
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure

      # the Hessian matrix for GEE model
      Hessian.alpha = Hessian.alpha + crossprod(t(t.D.i), tcrossprod(ginv(V.i), t.D.i))


    }

    return (Hessian.alpha)


  }


}

# Get the permutation statistics by R
Score.test.stat.zero <- function(Y0, Z, Z.par.index, class.id = NULL){

  # Z.reduce preserve covariates that are not interested in
  Z.reduce = Z[,-Z.par.index,drop=FALSE]
  n = nrow(Y0)
  m = ncol(Y0)
  p = ncol(Z)

  if(is.null(colnames(Y0))){
    colnames(Y0) <- paste0(rep('V',m),c(1:m))
  }
  p.reduce = ncol(Z.reduce)
  outcome = NULL
  id = NULL
  cova.reduce = NULL
  for(i in 1:n){

    outcome = c(outcome, Y0[i,])
    index.start = 1
    index.end = p

    index.start.reduce = 1
    index.end.reduce = p.reduce

    for(j in 1:m){

      tmp = rep(0, m*p.reduce)
      tmp[index.start.reduce:index.end.reduce] = Z.reduce[i,]
      cova.reduce = rbind(cova.reduce, tmp)
      index.start.reduce = index.start.reduce + p.reduce
      index.end.reduce = index.end.reduce + p.reduce

    }
    # id = c(id, rep(i, m))
  }

  # data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  # gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  # # use geeglm to estimate the parameter not interested in
  # data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  # gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  # wald.test = anova(gee.full, gee.reduce)

  data.reduce = data.frame(outcome=outcome, cova.reduce, row.names=NULL)
  gee.reduce = glm(formula = outcome ~ . - 1, family = binomial("logit"),  data = data.reduce)
  # gee.reduce = brglmFit(x = cova.reduce, y = outcome, family = binomial("logit"), control = list(type=c("AS_mean")), intercept = F)
  # gee.reduce = update(gee.reduce, method="brglmFit", type = "AS_mean")

  ########### perform score test
  n.alpha = m * p
  par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha)
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0 # set the est.reduce.aplha corresponding to parameters which are interested in to 0
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce

  data.alpha = list(Y=Y0, Z=Z)
  # estimate the Score statistics for parameter of interest
  tmp = fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list

  Hess.reduce.alpha =  fun.hessian.alpha(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix according to the index of par.interest.index.alpha
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ),
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))


  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]

  B1 <- ginv(Hess.reduce.reorg[(1:n.par.interest.alpha), (1:n.par.interest.alpha), drop = F] - crossprod(t(Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha), drop = F]),
                                                                                                         crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), ((n.par.interest.alpha + 1):n.alpha), drop = F])), Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), (1:n.par.interest.alpha), drop = F])))

  alpha.hat <- crossprod(t(B1), A)


  U <- Score.reduce.reorg[ ,1:n.par.interest.alpha] - crossprod(t(Score.reduce.reorg[ ,((n.par.interest.alpha + 1):n.alpha), drop = F]),tcrossprod(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha),
                                                                                                                                                                          ((n.par.interest.alpha + 1):n.alpha), drop = F]), Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha), drop = F] ))

  B2 <- matrix(0, n.par.interest.alpha, n.par.interest.alpha)

  if(is.null(class.id)){
    for (i in 1:n) {
      B2 <- B2 + tcrossprod(U[i, ])
    }
  }else{
    init_vec = rep(0, n.par.interest.alpha)
    grp.lvls = levels(class.id)
    grp.num <- length(grp.lvls)
    vec.list <- rep(list(init_vec),grp.num)
    for(i in 1:n){
      n.lvl = which(class.id[i] == grp.lvls)
      vec.list[[n.lvl]] = vec.list[[n.lvl]] + U[i, ]
    }
    vec.list <- lapply(vec.list, tcrossprod)
    B2 <- Reduce('+', vec.list)
  }

  cov.alpha = crossprod(t(B1), crossprod(t(B2), B1))

  # save these outcomes for later resampling test
  return(list(score.alpha = alpha.hat, est.cov.zero=cov.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list )   )

}

# Get the permutation statistics by cpp
Score.test.stat.zero.meta.4Gresampling <- function(Z.perm.list, Z.par.index, cluster.id.list.c, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, meta.method = "FE-MetaQCAT"){
  tmp = score_test_stat_zero_meta_resampling(Z.perm.list, cluster.id.list.c, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Z.par.index, n.par.interest.alpha)
  est.cov.meta = tmp$est_cov_meta
  score.alpha.meta = tmp$score_alpha_meta
  est.cov.zero = tmp$est_cov
  score.alpha = tmp$score_alpha
  # save the index of those elements that have values greater than zero in score.aplha.meta vector
  save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
  score.alpha.meta = score.alpha.meta[save.index.zero]
  est.cov.meta = est.cov.meta[save.index.zero,save.index.zero]
  cov.alpha.hat = ginv(est.cov.meta)
  stu.num = length(Z.perm.list)
  if(meta.method == "FE-MetaQCAT*"){
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(cov.alpha.hat), score.alpha.meta))
  }
  # if(meta.method == "FE-VC*"){
  #   score.stat.alpha.perm = crossprod(score.alpha.meta)
  # }
  if(meta.method == "RE-MetaQCAT*"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov.zero[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(cov.alpha.hat), score.alpha.meta)) + U.theta^2/V.theta
  }
  # if(meta.method == "RE-VC*"){
  #   U.theta = 0
  #   V.theta = 0
  #   U.tau = 0
  #   est.inv.sum = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha)
  #   for( i in 1:stu.num){
  #     idx = col.zero.index.list[[i]]
  #     est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
  #     U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
  #     # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
  #     V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
  #     est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
  #   }
  #   V.tau = 1/2 * tr(crossprod(est.inv.sum))
  #   U.tau = 1/2 * crossprod(score.alpha.meta) - 1/2 * tr(est.inv.sum)
  #   score.stat.alpha.perm = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
  # }
  return(as.numeric(score.stat.alpha.perm))

}

# add 10/22/2022 for adaptive resampling
resample.work.zero.meta <- function(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, cluster.id.list, permute.strata, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, meta.method = "FE-MetaQCAT"){
  stu.num = length(Z.list)
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  cls.flg = !is.null(cluster.id.list)
  # cls.num.vec = NULL
  cluster.id.list.c = NULL
  for(k in start.nperm:end.nperm){
    Z.perm.list = Z.list
    if(!cls.flg){
      # cls.num.vec = NULL
      for(p in 1:length(Z.list)){
        idx = sample(1:nrow(Z.list[[p]])) # sampling the index and reset the design matrix based on the new index
        Z.perm.list[[p]][,Z.par.index] = Z.list[[p]][idx,Z.par.index]
      }
    }else{
      cluster.id.list.c = lapply(cluster.id.list, function(x){if(is.null(x)){NULL}else{as.numeric(x)}})
      # cls.num.vec = sapply(cluster.id.list.c, function(x){length(unique(x))})
      if(permute.strata == "within"){
        for(p in 1:length(Z.list)){
          class.id = cluster.id.list[[p]]
          if(!is.null(class.id)){
            perm.idx = idx = seq(nrow(Z.list[[p]]))
            cls.lvls = levels(class.id)
            perm.idx.lst = tapply(idx, class.id, function(x){if(length(x)==1){x}else{sample(x)}})
            for(i in 1:length(cls.lvls)){
              Z.perm.list[[p]][class.id == cls.lvls[i], Z.par.index] = Z.list[[p]][perm.idx.lst[[cls.lvls[i]]], Z.par.index]
              #group.perm.index[group.perm.index == cls.lvls[i]] = perm.cls.lvls[i]
            }
          }else{
            idx = sample(1:nrow(Z.list[[p]])) # sampling the index and reset the design matrix based on the new index
            Z.perm.list[[p]][,Z.par.index] = Z.list[[p]][idx,Z.par.index]
          }
        }
      }else if(permute.strata == "between"){
        for(p in 1:length(Z.list)){
          # group.perm.index = cluster.id.list.c[[p]]
          # Z.perm.tmp = Z.tmp = Z.list[[p]][, Z.par.index]
          class.id = cluster.id.list[[p]]
          if(!is.null(class.id)){
            cls.lvls = levels(class.id)
            cls.cova.lst = lapply(cls.lvls, function(x){Z.list[[p]][class.id == x, Z.par.index, drop = F][1, , drop = F]})
            names(cls.cova.lst) <- cls.lvls
            perm.cls.lvls = if(length(cls.lvls) == 1){cls.lvls}else{sample(cls.lvls)}
            for(i in 1:length(cls.lvls)){
              if(sum(class.id == cls.lvls[i]) == 0){
                print('stop here')
              }
              if(is.null(cls.cova.lst[[perm.cls.lvls[i]]])){
                print('stop here')
              }
              Z.perm.list[[p]][class.id == cls.lvls[i], Z.par.index] = cls.cova.lst[[perm.cls.lvls[i]]][rep(1,sum(class.id == cls.lvls[i])), ]
              #group.perm.index[group.perm.index == cls.lvls[i]] = perm.cls.lvls[i]
            }
          }else{
            idx = sample(1:nrow(Z.list[[p]])) # sampling the index and reset the design matrix based on the new index
            Z.perm.list[[p]][,Z.par.index] = Z.list[[p]][idx,Z.par.index]
          }
          # cluster.id.list.c = as.numeric(group.perm.index) - 1
        }
      }
    }
    # get the permutation test statistics
    score.stat.alpha.perm <- try(Score.test.stat.zero.meta.4Gresampling(Z.perm.list, Z.par.index, cluster.id.list.c, n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, meta.method = meta.method))
    if (!("try-error" %in% class(score.stat.alpha.perm))) {# if score.stat.alpha.perm exists, then n.one.new + 1
      # if the permutation test statistics greater than original test statistics, cnt + 1
      n.zero.new = n.zero.new + 1
      if(score.stat.alpha.perm >= score.stat.zero.meta){
        zero.acc.new = zero.acc.new + 1 # if score.stat.alpha.perm >= score.stat.zero.meta one.acc.new + 1

      }
    }
  }

  # adaptive adjust the total number of iterations according to the number of resampling statistics which are more extreme than the original one in each loop
  # if the total number of permutation results which are greater than original results too small, enlarge the
  # number of total iterations
  if(zero.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;

  }else if(zero.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;
  }

  return(list(n.zero.new=n.zero.new, zero.acc.new=zero.acc.new,
              flag=flag, next.end.nperm=next.end.nperm))

}

# Z.list: a list of covariates for zero part: first column is always intercept
# Z.par.index: index for the parameter of interest for the Z part

Score.test.zero.meta <- function(Y.list, Z.list, Z.par.index, tax.names = NULL, cluster.id.list = NULL, permute.strata = NULL, seed=11, resample=FALSE, n.replicates=NULL, meta.method = "FE-MetaQCAT*"){
  stu.num = length(Z.list)
  n.par.interest = length(Z.par.index)
  m = ncol(Y.list[[1]])
  # p.zero = ncol(Z.list[[1]])
  Y0.list = Y.list
  # the total number of parameter of interest
  n.par.interest.alpha = m*length(Z.par.index)
  remove.study = NULL
  cls.flg = !is.null(cluster.id.list)
  for(i in 1:stu.num)
  {
    # for each study set those value which are zero to 1
    Y0.list[[i]][Y.list[[i]]==0] = 1
    # for each study set those value which are non-zero to 0
    Y0.list[[i]][Y.list[[i]]>0] = 0
    remove.rows.idx = which(rowSums(Y0.list[[i]]) == 0)
    if(length(remove.rows.idx) == nrow(Y0.list[[i]])){
      remove.study = append(remove.study, i)
      next
    }
  }
  if(!is.null(remove.study)){
    if(length(remove.study) == stu.num){
      stop("Error: Not proper data")
    }else{
      Y0.list <- Y0.list[-remove.study]
      Z.list <- Z.list[-remove.study]
      if(cls.flg){
        cluster.id.list <- cluster.id.list[-remove.study]
      }
      stu.num = length(Z.list)
    }
  }
  if(! meta.method %in% c("FE-MetaQCAT*", "RE-MetaQCAT*")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if (length(unique(sapply(1:stu.num,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  if(m<=1){
    stop("Error: Improper dimension for OTU table")
  }
  if(is.null(Z.par.index)){
    stop("Error: Testing parameters for the intercept is not informative. (Alpha part)")
  }
  ava.cnt = 0
  # col.pos.index.lst = lapply(1:stu.num, function(j) )
  ############################# Asymptotic: zero part


  # if all 0 in one group across across taxa, then output NA
  #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){

  # initialize for later use
  score.pvalue = score.stat.meta = df = NA
  var.alpha.hat = alpha.hat = save.index.zero = NA
  score.stat.alpha = NULL
  remove.index = NULL
  score.alpha = list()
  est.cov.zero = list()
  zero.vA.list.meta = list()
  zero.Vinv.list.meta = list()
  zero.VY.list.meta = list()
  col.zero.index.list = list()
  # initialize the score statistics and estimate covariance matrix for meta analysis
  score.alpha.meta  = rep(0,n.par.interest.alpha) ## A
  est.cov.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
  for(i in 1:stu.num){
    Y0 = Y0.list[[i]]
    Z = Z.list[[i]]
    zero.cnt.ave = colMeans(Y0)
    col.zero.index = which((zero.cnt.ave<0.95)&(zero.cnt.ave>0.05)) # keep those taxa which have both zero and one (both positive and zero observations)
    Y0 = Y0[, col.zero.index , drop=FALSE] # only save those columns which have both 0 and 1 values
    if(cls.flg){
      class.id = cluster.id.list[[i]]
    }else{
      class.id = NULL
    }
    if(length(col.zero.index)<1)
    {
      remove.index = append(remove.index,i)
      next
    }else{
      # nY0 = rowSums(Y0)
      # ## remove 03/28/2016
      # index.subj.zero = which(nY0>0)
      # if(length(index.subj.zero)== 0){
      #     next
      #}
      #
      # Y0 = Y0[index.subj.zero, , drop=FALSE]
      # Z = Z[index.subj.zero, , drop=FALSE]
      #
      tmp.zero = try( Score.test.stat.zero(Y0, Z, Z.par.index, class.id) )
    }
    if("try-error" %in% class(tmp.zero)){
      remove.index= append(remove.index,i)
      next

    }else{
      # number of study which can get score statistics
      ava.cnt = ava.cnt + 1
      # get the index for parameter of interest after score statistics and estimate covariance matrix are reorginzed
      # different across studies because of different column numbers
      idx = kronecker((col.zero.index-1)*n.par.interest, rep(1,n.par.interest)) + c(1:n.par.interest)
      col.zero.index.list[[ava.cnt]] = idx  # the index of alpha parameter of interest in this study
      # save these outcome in list form for later meta-analysis as well as resampling test
      score.stat.alpha = append(score.stat.alpha, tmp.zero$score.stat.alpha)
      score.alpha[[ava.cnt]] = tmp.zero$score.alpha
      est.cov.zero[[ava.cnt]] = tmp.zero$est.cov.zero
      zero.vA.list.meta[[ava.cnt]] = tmp.zero$vA.list
      zero.Vinv.list.meta[[ava.cnt]] = tmp.zero$Vinv.list
      zero.VY.list.meta[[ava.cnt]] = tmp.zero$VY.list
      score.alpha.meta[idx] =  score.alpha.meta[idx] + crossprod(ginv(tmp.zero$est.cov.zero), tmp.zero$score.alpha) # add according to the index of parameter of interest for each study
      est.cov.meta[idx, idx] =  est.cov.meta[idx, idx] + ginv(tmp.zero$est.cov.zero)
    }
  }

  ############################# Asymptotic: combined
  if(ava.cnt>0){
    if(length(remove.index) != 0){
      Z.list = Z.list[-remove.index]
      Y0.list = Y0.list[-remove.index]
      if(cls.flg){
        cluster.id.list = cluster.id.list[-remove.index]
      }
    }
    # save the index of those elements that have values greater than zero in score.aplha.meta vector
    save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
    n.par.save.alpha = length(save.index.zero)
    score.alpha.meta = score.alpha.meta[save.index.zero]
    est.cov.meta = est.cov.meta[save.index.zero,save.index.zero]
    cov.alpha.hat = ginv(est.cov.meta)
    alpha.hat = as.numeric(crossprod(cov.alpha.hat, score.alpha.meta))
    if(!is.null(tax.names)){
      names(alpha.hat) <- tax.names[save.index.zero]
    }
    if(meta.method == "FE-MetaQCAT*"){
      score.stat.meta = crossprod( score.alpha.meta,crossprod(t(cov.alpha.hat), score.alpha.meta))
      score.pvalue = 1- pchisq(score.stat.meta,df = n.par.save.alpha)
    }
    # if(meta.method == "FE-VC*"){
    #   weight.cov.meta = eigen(est.cov.meta)$values
    #   score.stat.meta = crossprod(score.alpha.meta)
    #   score.pvalue = davies(score.stat.meta, weight.cov.meta, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
    #   score.pvalue = ifelse(score.pvalue>0,score.pvalue,0)
    #   df = n.par.save.alpha
    # }
    if(meta.method == "RE-MetaQCAT*"){
      U.theta = 0
      V.theta = 0
      for( i in 1:ava.cnt){
        est.inv = ginv(est.cov.zero[[i]])
        U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      }
      score.stat.meta = crossprod( score.alpha.meta,crossprod(t(cov.alpha.hat), score.alpha.meta)) + U.theta^2/V.theta
    }
    # if(meta.method == "RE-VC*"){
    #   U.theta = 0
    #   V.theta = 0
    #   U.tau = 0
    #   est.inv.sum = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha)
    #   for( i in 1:ava.cnt){
    #     idx = col.zero.index.list[[i]]
    #     est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
    #     U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
    #     # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
    #     V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    #     est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
    #   }
    #   V.tau = 1/2 * tr(crossprod(est.inv.sum))
    #   U.tau = 1/2 * crossprod(score.alpha.meta) - 1/2 * tr(est.inv.sum)
    #   score.stat.meta = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
    # }
    var.alpha.hat = as.numeric(diag2(cov.alpha.hat))
    if(!is.null(tax.names)){
      names(var.alpha.hat) <- tax.names[save.index.zero]
    }
  }
  zero.results = list(score.stat = score.stat.meta, score.pvalue = score.pvalue, alpha.hat = alpha.hat, var.alpha.hat = var.alpha.hat, save.index.zero = save.index.zero)
  # if resample = TRUE then will apply permutation method to get permuted p-value
  # adaptive resampling test
  if(resample){

    #print("simulated stat:")
    set.seed(seed)
    if(!is.na(score.stat.meta)){

      n.zero = 0
      zero.acc = 0

      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      while(flag & end.nperm <= n.replicates){
        results = resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, cluster.id.list = cluster.id.list, permute.strata = permute.strata, score.stat.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, meta.method = meta.method)
        n.zero = results$n.zero.new
        zero.acc = results$zero.acc.new
        flag = results$flag
        next.end.nperm = results$next.end.nperm

        if(flag){
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;

        }

        if(start.nperm < n.replicates & end.nperm > n.replicates){
          #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))

          results = resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, cluster.id.list = cluster.id.list, permute.strata = permute.strata, score.stat.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, meta.method = meta.method)

          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new

        }

      }

      tmp = (zero.acc+1)/(n.zero+1) # to avoid n.one be zero # resampling p value


    }else{

      tmp = NA
    }

    zero.results = c(zero.results, score.Rpvalue = tmp)

  }
  return(zero.results)
}

lineage.score.test.meta <- function(tax.info, W.data.list = NULL,  W.data.rare.list = NULL, otucols, X, X.index, class.id, permute.strata = NULL, n.replicates=NULL, meta.method = "FE-MetaQCAT"){
  Rank.low = unname(tax.info['Rank.low'])
  Rank.high = unname(tax.info['Rank.high'])
  level.uni = unname(tax.info['level.uni'])
  n.OTU = length(X)
  X.rare = X
  class.id.rare = class.id
  # partition and merge OTU table according to taxonomy information
  if(!is.null(W.data.list)){
    tt = lapply(1:n.OTU, function(j) W.data.list[[j]][, lapply(.SD , sum, na.rm=TRUE), .SDcols=as.vector(unlist(otucols[j])), by=list( get(Rank.low), get(Rank.high) )])
    tt = lapply(1:n.OTU,function(j) setnames(tt[[j]], 1:2, c(Rank.low, Rank.high)))
    W.count = lapply(1:n.OTU,function(j) tt[[j]][, otucols[[j]], with=FALSE])
    W.tax = as.vector(unlist(tt[[1]][, Rank.low, with=FALSE]))
    tax.names <- as.vector(unlist(tt[[1]][, Rank.high, with=FALSE]))[which(W.tax == level.uni)]
    Y = lapply(1:n.OTU, function(i) t(W.count[[i]][which(W.tax == level.uni), , drop=FALSE]))
    if(is.null(Y)){
      return(NA)
    }
    if(ncol(Y[[1]])==1){
      return(NA)
    }
    zero.ave.lst <- lapply(Y, function(x){colMeans(x==0)})
    zero.cnt = as.matrix(Reduce('+', zero.ave.lst), nrow = 1)
    if(length(zero.ave.lst) > 1){
      zero.cnt.mat <- Reduce(rbind, zero.ave.lst)
    }else{
      zero.cnt.mat <- as.matrix(zero.ave.lst[[1]], nrow = 1)
    }
    non.base.id <- apply(zero.cnt.mat, 2, function(x){any(x == 1)})
    if(sum(non.base.id) == ncol(zero.cnt.mat)){
      remove.study.cnt <- apply(zero.cnt.mat, 2, function(x){sum(x == 1)})
      reference.column.candidates = which(remove.study.cnt == min(remove.study.cnt))
      reference.column = reference.column.candidates[rev(which(zero.cnt[reference.column.candidates] == min(zero.cnt[reference.column.candidates])))[1]]
      remove.study.id = which(zero.cnt.mat[,reference.column] == 1)
      Y = Y[-remove.study.id]
      X = X[-remove.study.id]
      if(!is.null(class.id)){
        class.id = class.id[-remove.study.id]
      }
      # return(NA) # no proper base for positive part
    }else if(sum(non.base.id)>0){
      zero.cnt[non.base.id] <- Inf
      reference.column = rev(which(zero.cnt == min(zero.cnt)))[1]
    }else{
      reference.column = rev(which(zero.cnt == min(zero.cnt)))[1]
    }

    re.order <- c(c(1:ncol(Y[[1]]))[-reference.column], reference.column)
    Y =  lapply(Y, function(x){x[, re.order]})# perform one test using all OTUs
    # Y =  lapply(Y, function(x){x[, order(-zero.cnt)]})
    tax.names.reorder = tax.names[re.order]
    reference.taxon <- tax.names[reference.column]
  }
  if(!is.null(W.data.rare.list)){
    # partition and merge OTU table according to taxonomy information
    tt.rare = lapply(1:n.OTU, function(j) W.data.rare.list[[j]][, lapply(.SD , sum, na.rm=TRUE), .SDcols=as.vector(unlist(otucols[j])), by=list( get(Rank.low), get(Rank.high) )])
    tt.rare = lapply(1:n.OTU,function(j) setnames(tt.rare[[j]], 1:2, c(Rank.low, Rank.high)))
    W.rare.count = lapply(1:n.OTU,function(j) tt.rare[[j]][, otucols[[j]], with=FALSE])
    W.tax = as.vector(unlist(tt.rare[[1]][, Rank.low, with=FALSE]))
    Y.rare = lapply(1:n.OTU, function(i) t(W.rare.count[[i]][which(W.tax == level.uni), , drop=FALSE]))
    if(ncol(Y.rare[[1]])==1){
      return(NA)
    }
    tax.names.rare <- as.vector(unlist(tt.rare[[1]][, Rank.high, with=FALSE]))[which(W.tax == level.uni)]
  }
  # subtree = c(subtree, level.uni[j])
  #print(level.uni[j])
  if(is.null(n.replicates)){ # asymptotic test only
    # (Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, meta.method = "FE-MetaQCAT", Weight=NULL )
    #  run test for each lineage
    if(meta.method == "FE-MetaQCAT"){
      tmp = try(Score.test.meta(Y, X, X.index, tax.names.reorder, class.id, permute.strata, meta.method = meta.method))
      if(!("try-error" %in% class(tmp))){
        pval = tmp$score.pvalue
        # re.map <- order(match(names(tmp$beta.hat), tax.names))
        pos.coef <- list(est.par = tmp$beta.hat, est.par.var = tmp$var.beta.hat, reference.taxon = reference.taxon)
        res <- list(pval = pval, pos.coef = pos.coef)
        return(res)
      }else{
        return(NA)
      }
    }else if(meta.method == "FE-MetaQCAT*"){
      tmp = try(Score.test.zero.meta(Y.rare, X.rare, X.index, tax.names.rare, class.id.rare, permute.strata, meta.method = meta.method))
      if(!("try-error" %in% class(tmp))){
        pval = tmp$score.pvalue
        zero.coef <- list(est.par = tmp$alpha.hat, est.par.var = tmp$var.alpha.hat)
        res <- list(pval = pval, zero.coef = zero.coef)
        return(res)
      }else{
        return(NA)
      }
    }else if(meta.method == 'FE-MetaQCAT-O'){
      Pval_c = c()
      pos.coef = zero.coef = NULL
      tmp = try(Score.test.meta(Y, X, X.index, tax.names.reorder, class.id, permute.strata, meta.method = "FE-MetaQCAT"))
      if(!("try-error" %in% class(tmp))){
        pos.coef <- list(est.par = tmp$beta.hat, est.par.var = tmp$var.beta.hat, reference.taxon = reference.taxon)
        Pval_c = append(Pval_c, tmp$score.pvalue)
      }
      tmp = try(Score.test.zero.meta(Y.rare, X.rare, X.index, tax.names.rare, class.id.rare, permute.strata, meta.method = "FE-MetaQCAT*"))
      if(!("try-error" %in% class(tmp))){
        zero.coef <- list(est.par = tmp$alpha.hat, est.par.var = tmp$var.alpha.hat)
        Pval_c = append(Pval_c, tmp$score.pvalue)
      }
      Pval_c = na.omit(Pval_c)
      if(length(Pval_c) == 0){
        return(NA)
      }else{
        res <- list(pval = ACAT(Pval_c), pos.coef = pos.coef, zero.coef = zero.coef)
        return(res)
      }
    }
  }
  else{
    # if n.perm in not null, select the significant lineage according to the resampling pvalue
    #  run test for each lineage
    if(meta.method %in% c("FE-MetaQCAT", "RE-MetaQCAT")){
      tmp = try(Score.test.meta(Y, X, X.index, tax.names.reorder, class.id, permute.strata, resample=TRUE, n.replicates=n.replicates, meta.method = meta.method))
      if(!("try-error" %in% class(tmp))){
        if(meta.method == "RE-MetaQCAT"){
          pval = tmp$score.Rpvalue
          pos.coef <- list(est.par = tmp$beta.hat, est.par.var = tmp$var.beta.hat, reference.taxon = reference.taxon)
          res <- list(pval = pval, pos.coef = pos.coef)
          return(res)
        }else{
          pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
          if(all(is.na(pval))){
            pval <- NA
          }
          pos.coef <- list(est.par = tmp$beta.hat, est.par.var = tmp$var.beta.hat, reference.taxon = reference.taxon)
          res <- list(pval = pval, pos.coef = pos.coef)
          return(res)
        }
      }else{
        return(NA)
      }
    }else if(meta.method %in% c("FE-MetaQCAT*", "RE-MetaQCAT*")){
      tmp = try(Score.test.zero.meta(Y.rare, X.rare, X.index, tax.names.rare, class.id.rare, permute.strata, seed=11, resample=TRUE, n.replicates=n.replicates, meta.method = meta.method))
      if(!("try-error" %in% class(tmp))){
        if(meta.method == "RE-MetaQCAT*"){
          pval = tmp$score.Rpvalue
          zero.coef <- list(est.par = tmp$alpha.hat, est.par.var = tmp$var.alpha.hat)
          res <- list(pval = pval, zero.coef = zero.coef)
          return(res)
        }else{
          pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
          zero.coef <- list(est.par = tmp$alpha.hat, est.par.var = tmp$var.alpha.hat)
          res <- list(pval = pval, zero.coef = zero.coef)
          return(res)
        }
      }else{
        return(NA)
      }
    }else if(meta.method == 'FE-MetaQCAT-O'){
      Pval_asy_c = c()
      Pval_res_c = c()
      pos.coef = zero.coef = NULL
      tmp = try(Score.test.meta(Y, X, X.index, tax.names.reorder, class.id, permute.strata, resample=TRUE, n.replicates=n.replicates, meta.method = "FE-MetaQCAT"))
      if(!("try-error" %in% class(tmp))){
        Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
        Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
        pos.coef <- list(est.par = tmp$beta.hat, est.par.var = tmp$var.beta.hat, reference.taxon = reference.taxon)
      }
      tmp = try(Score.test.zero.meta(Y.rare, X.rare, X.index, tax.names.rare, class.id.rare, permute.strata, resample=TRUE, n.replicates=n.replicates, meta.method = "FE-MetaQCAT*"))
      if(!("try-error" %in% class(tmp))){
        Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
        Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
        zero.coef <- list(est.par = tmp$alpha.hat, est.par.var = tmp$var.alpha.hat)
      }
      Pval_asy_c = na.omit(Pval_asy_c)
      Pval_res_c = na.omit(Pval_res_c)
      if((length(Pval_asy_c) == 0)|(length(Pval_res_c) == 0)){
        return(NA)
      }else{
        pval = c(ACAT(Pval_asy_c),ACAT(Pval_res_c))
        res <- list(pval = pval, pos.coef = pos.coef, zero.coef = zero.coef)
        return(res)
      }
    }else if(meta.method == "RE-MetaQCAT-O"){
      Pval_res_c = c()
      pos.coef = zero.coef = NULL
      tmp = try(Score.test.meta(Y, X, X.index, tax.names.reorder, class.id, permute.strata, resample=TRUE, n.replicates=n.replicates, meta.method = "RE-MetaQCAT"))
      if(!("try-error" %in% class(tmp))){
        Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
        pos.coef <- list(est.par = tmp$beta.hat, est.par.var = tmp$var.beta.hat, reference.taxon = reference.taxon)
      }
      tmp = try(Score.test.zero.meta(Y.rare, X.rare, X.index, tax.names.rare, class.id.rare, permute.strata, seed=11, resample=TRUE, n.replicates=n.replicates, meta.method = "RE-MetaQCAT*"))
      if(!("try-error" %in% class(tmp))){
        Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
        zero.coef <- list(est.par = tmp$alpha.hat, est.par.var = tmp$var.alpha.hat)
      }
      Pval_res_c = na.omit(Pval_res_c)
      if(length(Pval_res_c) == 0){
        return(NA)
      }else{
        pval = c(ACAT(Pval_res_c))
        res <- list(pval = pval, pos.coef = pos.coef, zero.coef = zero.coef)
        return(res)
      }
    }
  }

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
    warning("There are p-values that are exactly 1!")
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


Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab)), seed = 1234){
  # Rarefaction function: downsample to equal depth
  #
  # Args:
  #		otu.tab: OTU count table, row - n sample, column - q OTU
  #		depth: required sequencing depth
  #
  # Returns:
  # 	otu.tab.rff: Rarefied OTU table
  #		discard: labels of discarded samples
  #
  otu.tab <- as.matrix(otu.tab)
  ind <- (rowSums(otu.tab) < depth)
  sam.discard <- rownames(otu.tab)[ind]
  otu.tab <- otu.tab[!ind, ]
  set.seed(seed)
  rarefy <- function(x, depth){
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    z
  }
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
  rownames(otu.tab.rff) <- rownames(otu.tab)
  colnames(otu.tab.rff) <- colnames(otu.tab)
  return(list(otu.tab.rff=otu.tab.rff, discard=sam.discard))
}



#' meta-analysis version of Quasi-Conditional Association Tests for microbiome data
#'
#' @param OTU a list of matrices contains OTU counts with each row corresponding to a sample and each column corresponding to an OTU or taxa. Column names, denoting the species names of each taxon, are mandatory.
#' @param X a list of matrices contains covariates, with each column pertaining to one variable. The number of elements in X and OTU must match.
#' @param X.index a vector indicates the columns in X for the covariate(s) of interest.
#' @param meta.method The meta-analysis method to be used. Including composition data based fixed effect test "FE-MetaQCAT", presence-absence data based fixed effect test "FE-MetaQCAT*" and the Cauchy-combined p-value of the two fix effect tests denoted as "FE-MetaQCAT-O",
#'                  composition data based random effect test "RE-MetaQCAT", presence-absence data based random effect test "RE-MetaQCAT*"  and the Cauchy-combined p-value of all the random effect tests denoted as "RE-MetaQCAT-O".
#' @param class.id a list of factors represents the group index to which each sample belongs if samples are correlated. The default is null, but this must be provided if samples are correlated.
#' @param Tax a matrix defines the taxonomy ranks , with each row corresponding to an OTU or taxa and each column to a rank (starting from the higher taxonomic level). Row names are mandatory and should be consistent with the column names of the OTU table. Column names should be formatted as "Rank1", "Rank2", etc.
#'        If provided, tests will be performed for lineages based on the taxonomic rank. The output includes P-values for all lineages; a list of significant lineages controlling the false discovery rate (based on the resampling p-value if a resampling test is performed).
#'        If not provided, a single test will be performed with all OTUs, and one p-value will be output.
#' @param min.depth keep samples with depths >= min.depth.
#' @param max.zero.prop  keep taxa with percentage of zero counts <= max.zero.prop.
#' @param n.perm perform asymptotic test if n.perm is null, otherwise, perform permutation tests using the specified number of resamplings. If using random effect meta-analysis methods, then a resampling test must be performed
#' @param permute.strata a specification for how to perform resampling test for correlated samples. The default is null. The shuffle will be held within each group and among groups if the value is "within" and "between", respectively.
#' @param fdr.alpha false discovery rate for multiple tests on the lineages.
#' @param n.threads number of cores for parallel computation, with default being one.
#'
#' @return If Tax is null, only the global p-value for all OTUs will be returned. If Tax is provided, the output will include:
#'    \item{lineage.pval}{p-values for all lineages. By default ( meta.method = "FE-MetaQCAT", n.perm = NULL ), only the asymptotic test will be performed.}
#'    \item{sig.lineage}{a vector of significant lineages.}
#'    \item{lineage.coef}{the estimated coefficent for each taxon, including the estimated parameter and the its variance. If provided method is composition data based, the taxon regarded as the reference will also be returned.}
#'
#' @export
#'
#' @import MASS
#' @import data.table
#' @import CompQuadForm
#' @import brglm2
#' @import psych
#' @import parallel
#' @import stats
#' @importFrom dplyr bind_rows
#' @importFrom parallelly availableCores
#' @references
#' Tang ZZ, Chen G, Alekseyenko AV, Li H. (2017) A general framework for association analysis of microbial communities on a taxonomic tree.
#' \emph{Bioinformatics}
#' \doi{10.1093/bioinformatics/btw804}.
#' @references
#' Lee S, Teslovich TM, Boehnke M, Lin X. (2013) General framework for meta-analysis of rare variants in sequencing association studies. Am J Hum Genet.
#' \emph{Am J Hum Genet}
#' \doi{10.1016/j.ajhg.2013.05.010}.
#' @references
#' Benjamini, Yoav, and Yosef Hochberg.(1995) Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
#' \emph{Journal of the Royal Statistical Society. Series B}
MetaQCAT <- function(OTU, X, X.index,  meta.method = "FE-MetaQCAT", class.id = NULL, Tax=NULL, min.depth=0, max.zero.prop = 0.9, n.perm=NULL, permute.strata = NULL, fdr.alpha=0.05, n.threads = 1){
  n.OTU = length(OTU)
  n.X = length(X)
  # class.id
  # permuate.strata
  group.flag = !is.null(class.id)
  # drop.col = NULL
  if(meta.method %in% c("RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")){
    if(is.null(n.perm)){
      stop("The p-value for random effect meta-analysis method must be got by resampling test")
    }
    if(group.flag){
      if(! permute.strata %in% c("within", "between")){
        stop("Please choose a proper permutation method for cluster data")
      }
    }
  }

  if(group.flag){
    n.id = length(class.id)
    if(n.id != n.OTU)
    {
      stop("The length of OTU table and class.id should be the same once class.id is provided")
    }
  }

  if(n.OTU != n.X)
  {
    stop("The length of OTU table and Covariate should be the same")
  }
  remove.study = NULL
  # remove.col.lst = vector(mode = 'list', length = n.OTU)
  for(i in 1:n.OTU)
  {
    if(!is.matrix(OTU[[i]])){
      OTU[[i]] = as.matrix(OTU[[i]])
    }

    if(!is.matrix(X[[i]])){
      X[[i]] = as.matrix(X[[i]])
    }
    if(group.flag){
      if(!is.null(class.id[[i]])){
        if(!is.factor(class.id[[i]])){
          class.id[[i]] = as.factor(class.id[[i]])
        }
      }
    }

    if(nrow(OTU[[i]])!=nrow(X[[i]])){
      stop(paste0("Number of samples in the OTU table and the covariate table of study ", i,
                  " should be the same"))
    }
    remove.subject = which(rowSums(OTU[[i]])<=min.depth)
    if(length(remove.subject)>0){
      cat(paste("Remove",length(remove.subject), "samples with read depth less or equal to", min.depth, "in OTU table", i, "\n"))
      X[[i]] = X[[i]][-remove.subject, ,drop=FALSE]
      OTU[[i]] = OTU[[i]][-remove.subject, ,drop=FALSE]
      if(!is.null(class.id[[i]])){
        class.id[[i]] = factor(class.id[[i]][-remove.subject])
      }
    }
    remove.col = which(colSums(OTU[[i]] == 0)/nrow(OTU[[i]]) > max.zero.prop)
    if(length(remove.col)>0){
      cat(paste0("Remove ",length(remove.col), " taxa with more than ", max.zero.prop*100,"% of zero counts in OTU table ", i, "\n"))
      OTU[[i]] = OTU[[i]][ , -remove.col, drop=FALSE]
    }
    if((length(remove.col) == ncol(OTU[[i]]))|(length(remove.subject) == nrow(OTU[[i]]))){
      remove.study = append(remove.study, i)
    }
    # drop.col = union(drop.col,which(colSums(OTU[[i]])==0))
  }
  if(length(remove.study) == n.OTU){
    stop("Please provide proper data")
  }
  if(missing(X.index)){
    X.index = 1:ncol(X[[1]])
  }
  if(!is.null(remove.study)){
    X = X[-remove.study]
    OTU = OTU[-remove.study]
    class.id = class.id[-remove.study]
    n.OTU = length(OTU)
  }
  # join all the OTU table
  OTU.comb = as.data.frame(OTU[[1]])
  if(n.OTU>1){
    for (i in 2:n.OTU){
      OTU.comb = bind_rows(OTU.comb,as.data.frame(OTU[[i]]))
    }
  }
  OTU.comb[is.na(OTU.comb)] <- 0 # substitute NA to zero (do not matter because those columns will be deleted during calculation)
  OTU.comb <- as.matrix(OTU.comb)
  taxa.names = colnames(OTU.comb)
  if(!is.null(Tax)){
    # preserve those columns which have taxonomy information
    col.save = intersect(rownames(Tax),taxa.names)
    OTU.comb = OTU.comb[,taxa.names %in% col.save]
    Tax = Tax[rownames(Tax) %in% col.save,]
    Tax = Tax[match(taxa.names, rownames(Tax)),]
  }
  # divide the combined OTU table according to the number of observations for each study
  batch.cnt <- unlist(lapply(OTU, function(x) nrow(x)))
  batch.cnt <- append(1,batch.cnt)
  batch.cnt <- cumsum(batch.cnt)
  count = list()
  for (i in 1:n.OTU) {
    count[[i]] <- OTU.comb[batch.cnt[i]:(batch.cnt[i+1]-1),]
  }
  X = lapply(1:n.OTU,function(j) cbind(1, X[[j]])) # add the intercept term
  X.index = X.index + 1
  rarefy = non.rarefy = F
  if(meta.method %in% c("FE-MetaQCAT*",  "RE-MetaQCAT*")){
    rarefy = T
    non.rarefy = F
    count.rare = lapply(count,function(X){Rarefy(X)$otu.tab.rff})
  }else if(meta.method %in% c("FE-MetaQCAT", "RE-MetaQCAT")){
    rarefy = F
    non.rarefy = T
  }else{
    rarefy = non.rarefy = T
    count.rare = lapply(count,function(X){Rarefy(X)$otu.tab.rff})
  }
  if(is.null(Tax)){
    zero.ave.lst <- lapply(count, function(x){colMeans(x==0)})
    zero.cnt = Reduce('+', zero.ave.lst)
    if(length(zero.ave.lst) > 1){
      zero.cnt.mat <- Reduce(rbind, zero.ave.lst)
    }else{
      zero.cnt.mat <- as.matrix(zero.ave.lst[[1]], nrow = 1)
    }
    non.base.id <- apply(zero.cnt.mat, 2, function(x){any(x == 1)})
    if(sum(non.base.id)>0){
      zero.cnt[non.base.id] <- Inf
    }
    base.column = rev(which(zero.cnt == min(zero.cnt)))[1]
    count =  lapply(count, function(x){x[, c(c(1:ncol(x))[-base.column], base.column)]})# perform one test using all OTUs

    if(is.null(n.perm)){
      if(meta.method == "FE-MetaQCAT"){
        if(sum(non.base.id) == ncol(zero.cnt.mat)){
          stop('can not find proper base')
        }
        tmp = try(Score.test.meta(count, X, X.index, tax.names = NULL, class.id, permute.strata, meta.method = meta.method))
        if(!("try-error" %in% class(tmp))){
          pval = c(tmp$score.pvalue)
          names(pval) = paste0("Asymptotic-",meta.method)
          est.par.index = tmp$save.index.pos
          est.par.index[est.par.index>=base.column] = est.par.index[est.par.index>=base.column] + 1
          est.par = tmp$beta.hat
          est.par.var = tmp$var.beta.hat
          names(est.par) = names(est.par.var) =  taxa.names[est.par.index]
          pos.coef <- list(est.par = est.par, est.par.var = est.par.var, reference.taxon = taxa.names[base.column])
        }else{
          pval = NA
          pos.coef = NA
          names(pval) = paste0("Asymptotic-",meta.method)
        }
        res <- list(pval = pval, pos.coef = pos.coef)
      }else if(meta.method == "FE-MetaQCAT*"){
        tmp = try(Score.test.zero.meta(count.rare, X, X.index, tax.names = NULL, class.id, permute.strata, meta.method = meta.method))
        if(!("try-error" %in% class(tmp))){
          pval = c( tmp$score.pvalue )
          names(pval) = paste0("Asymptotic-",meta.method)
          est.par = tmp$alpha.hat
          est.par.var = tmp$var.alpha.hat
          names(est.par) = names(est.par.var) =  taxa.names[tmp$save.index.zero]
          zero.coef <- list(est.par = est.par, est.par.var = est.par.var)
        }else{
          pval = NA
          zero.coef = NA
          names(pval) = paste0("Asymptotic-",meta.method)
        }
        res <- list(pval = pval, zero.coef = zero.coef)
      }else if(meta.method == "FE-MetaQCAT-O"){
        Pval_c = c()
        for(i in c("FE-MetaQCAT","FE-MetaQCAT*")){
          if(i == "FE-MetaQCAT"){
            if(sum(non.base.id) == ncol(zero.cnt.mat)){
              next
            }
            tmp = try(Score.test.meta(count, X, X.index, tax.names = NULL, class.id, permute.strata, meta.method = i))
            pos.coef <- NA
            if(!("try-error" %in% class(tmp))){
              Pval_c = append(Pval_c, tmp$score.pvalue)
              est.par.index = tmp$save.index.pos
              est.par.index[est.par.index>=base.column] = est.par.index[est.par.index>=base.column] + 1
              est.par = tmp$beta.hat
              est.par.var = tmp$var.beta.hat
              names(est.par) = names(est.par.var) =  taxa.names[est.par.index]
              pos.coef <- list(est.par = est.par, est.par.var = est.par.var, reference.taxon = taxa.names[base.column])
            }
          }else{
            tmp = try(Score.test.zero.meta(count.rare, X, X.index, tax.names = NULL, class.id, permute.strata, meta.method = i))
            zero.coef <- NA
            if(!("try-error" %in% class(tmp))){
              Pval_c = append(Pval_c, tmp$score.pvalue)
              est.par = tmp$alpha.hat
              est.par.var = tmp$var.alpha.hat
              names(est.par) = names(est.par.var) =  taxa.names[tmp$save.index.zero]
              zero.coef <- list(est.par = est.par, est.par.var = est.par.var)
            }
          }
        }
        Pval_c = na.omit(Pval_c)
        if(length(Pval_c) == 0){
          pval = NA
          names(pval) = paste0("Asymptotic-",meta.method)
        }else{
          pval = ACAT(Pval_c)
          names(pval) = paste0("Asymptotic-",meta.method)
        }
        res <- list(pval = pval, pos.coef = pos.coef, zero.coef = zero.coef)
      }
    }else{ # resampling test + asymptotic test
      if(meta.method %in% c("FE-MetaQCAT", "RE-MetaQCAT")){
        if(sum(non.base.id) == ncol(zero.cnt.mat)){
          stop('can not find proper base')
        }
        tmp = try(Score.test.meta(count, X, X.index, tax.names = NULL, class.id, permute.strata, resample=TRUE, n.replicates=n.perm, meta.method = meta.method))
        if(!("try-error" %in% class(tmp))){
          est.par.index = tmp$save.index.pos
          est.par.index[est.par.index>=base.column] = est.par.index[est.par.index>=base.column] + 1
          est.par = tmp$beta.hat
          est.par.var = tmp$var.beta.hat
          names(est.par) = names(est.par.var) =  taxa.names[est.par.index]
          pos.coef <- list(est.par = est.par, est.par.var = est.par.var, reference.taxon = taxa.names[base.column])
          if(meta.method == "RE-MetaQCAT"){
            pval = c(tmp$score.Rpvalue)
            names(pval) = paste0("Resampling-",meta.method)
          }else{
            pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
            names(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
          }
        }else{
          pos.coef <- NA
          if(meta.method == "RE-MetaQCAT"){
            pval = NA
            names(pval) = paste0("Resampling-",meta.method)
          }else{
            pval = c(NA, NA)
            names(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
          }
        }
        res <- list(pval=pval, pos.coef = pos.coef)
      }else if(meta.method %in% c("FE-MetaQCAT*", "RE-MetaQCAT*")){
        tmp = try(Score.test.zero.meta(count.rare, X, X.index, tax.names = NULL, class.id, permute.strata, resample=TRUE, n.replicates=n.perm, meta.method = meta.method))
        if(!("try-error" %in% class(tmp))){
          est.par = tmp$alpha.hat
          est.par.var = tmp$var.alpha.hat
          names(est.par) = names(est.par.var) =  taxa.names[tmp$save.index.zero]
          zero.coef <- list(est.par = est.par, est.par.var = est.par.var)
          if(meta.method == "RE-MetaQCAT*"){
            pval = c(tmp$score.Rpvalue)
            names(pval) = paste0("Resampling-",meta.method)
          }else{
            pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
            names(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
          }
        }else{
          zero.coef <- NA
          if(meta.method == "RE-MetaQCAT*"){
            pval = NA
            names(pval) = paste0("Resampling-",meta.method)
          }else{
            pval = c(NA, NA)
            names(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
          }
        }
        res <- list(pval=pval, zero.coef = zero.coef)
      }else if(meta.method == 'FE-MetaQCAT-O'){
        Pval_asy_c = c()
        Pval_res_c = c()
        for(i in c("FE-MetaQCAT", "FE-MetaQCAT*")){
          if(i == "FE-MetaQCAT"){
            if(sum(non.base.id) == ncol(zero.cnt.mat)){
              next
            }
            tmp = try(Score.test.meta(count, X, X.index, tax.names = NULL, class.id, permute.strata, resample=TRUE, n.replicates=n.perm, meta.method = i))
            pos.coef <- NA
            if(!("try-error" %in% class(tmp))){
              est.par.index = tmp$save.index.pos
              est.par.index[est.par.index>=base.column] = est.par.index[est.par.index>=base.column] + 1
              est.par = tmp$beta.hat
              est.par.var = tmp$var.beta.hat
              names(est.par) = names(est.par.var) =  taxa.names[est.par.index]
              pos.coef <- list(est.par = est.par, est.par.var = est.par.var, reference.taxon = taxa.names[base.column])
              Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
              Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
            }
          }else{
            tmp = try(Score.test.zero.meta(count.rare, X, X.index, tax.names = NULL, class.id, permute.strata, resample=TRUE, n.replicates=n.perm, meta.method = i))
            zero.coef <- NA
            if(!("try-error" %in% class(tmp))){
              est.par = tmp$alpha.hat
              est.par.var = tmp$var.alpha.hat
              names(est.par) = names(est.par.var) =  taxa.names[tmp$save.index.zero]
              zero.coef <- list(est.par = est.par, est.par.var = est.par.var)
              Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
              Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
            }
          }
        }
        Pval_asy_c = na.omit(Pval_asy_c)
        Pval_res_c = na.omit(Pval_res_c)
        if((length(Pval_asy_c) == 0)|(length(Pval_res_c) == 0)){
          pval = c(NA, NA)
          names(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
        }else{
          pval = c(ACAT(Pval_asy_c),ACAT(Pval_res_c))
          names(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
        }
        res <- list(pval=pval, pos.coef = pos.coef, zero.coef = zero.coef)
      }else if(meta.method == "RE-MetaQCAT-O"){
        Pval_res_c = c()
        for(i in c("RE-MetaQCAT", "RE-MetaQCAT*")){
          if(i == "RE-MetaQCAT"){
            if(sum(non.base.id) == ncol(zero.cnt.mat)){
              next
            }
            tmp = try(Score.test.meta(count, X, X.index, tax.names = NULL, class.id, permute.strata, resample=TRUE, n.replicates=n.perm, meta.method = i))
            pos.coef <- NA
            if(!("try-error" %in% class(tmp))){
              est.par.index = tmp$save.index.pos
              est.par.index[est.par.index>=base.column] = est.par.index[est.par.index>=base.column] + 1
              est.par = tmp$beta.hat
              est.par.var = tmp$var.beta.hat
              names(est.par) = names(est.par.var) =  taxa.names[est.par.index]
              pos.coef <- list(est.par = est.par, est.par.var = est.par.var, reference.taxon = taxa.names[base.column])
              Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
            }
          }else{
            tmp = try(Score.test.zero.meta(count.rare, X, X.index, tax.names = NULL, class.id, permute.strata, resample=TRUE, n.replicates=n.perm, meta.method = i))
            zero.coef <- NA
            if(!("try-error" %in% class(tmp))){
              est.par = tmp$alpha.hat
              est.par.var = tmp$var.alpha.hat
              names(est.par) = names(est.par.var) =  taxa.names[tmp$save.index.zero]
              zero.coef <- list(est.par = est.par, est.par.var = est.par.var)
              Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
            }
          }
        }
        Pval_res_c = na.omit(Pval_res_c)
        if(length(Pval_res_c) == 0){
          pval = c(NA)
          names(pval) = c(paste0("Resampling-",meta.method))
        }else{
          pval = ACAT(Pval_res_c)
          names(pval) = c(paste0("Resampling-",meta.method))
        }
        res <- list(pval=pval, pos.coef = pos.coef, zero.coef = zero.coef)
      }
    }
    return(res)
  }else{ # perform tests for lineages

    if(!is.matrix(Tax)){
      Tax = as.matrix(Tax)
    }

    for(i in 1:n.OTU)
    {
      if( sum(!(colnames(count[[i]]) %in% rownames(Tax)))>0 ){
        stop(paste0("Error: OTU IDs in OTU table ",i," are not consistent with OTU IDs in Tax table"))
      }
    }

    n.rank = ncol(Tax)
    # merge the tax and count together for later partition and merge OTU table according to taxonomy information
    W.data.list = NULL
    W.data.rare.list = NULL
    if(non.rarefy){
      W.data.list = lapply(1:n.OTU,function(j) data.table(data.frame(Tax, t(count[[j]]))))
    }
    if(rarefy){
      W.data.rare.list = lapply(1:n.OTU,function(j) data.table(data.frame(Tax, t(count.rare[[j]]))))
    }

    if(non.rarefy){
      otucols = lapply(1:n.OTU,function(j) names(W.data.list[[j]])[-(1:n.rank)])
    }else{
      otucols = lapply(1:n.OTU,function(j) names(W.data.rare.list[[j]])[-(1:n.rank)])
    }

    n.level = n.rank-1

    subtree = NULL
    pval = NULL
    lineage.coef = list()
    tax.tab = NULL
    for(k in 1:n.level){
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")

      tmp = table(Tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)
      tmp = cbind(Rank.low, Rank.high,level.uni)
      tax.tab = rbind(tax.tab, tmp)
    }
    tax.info.lst <- lapply(seq_len(nrow(tax.tab)), function(i) tax.tab[i,])
    n.cores = availableCores()
    if(n.threads > n.cores){
      n.threads = n.cores
      warning("n.threads is larger than the number of available cores")
    }
    if(n.threads > 1){
      res_lst <- mclapply(tax.info.lst, lineage.score.test.meta, W.data.list, W.data.rare.list, otucols, X, X.index, class.id, permute.strata, n.perm, meta.method = meta.method, mc.cores = n.threads)
      k = 1
      for(j in 1:length(tax.info.lst)){
        if(any(is.na(res_lst[[j]]))|is.null(res_lst[[j]])){
          next
        }else{
          subtree = c(subtree, unname(tax.info.lst[[j]]["level.uni"]))
          pval = cbind(pval, res_lst[[j]]$pval)
          if(is.null(res_lst[[j]]$pos.coef)){
            lineage.coef[[k]] = list(zero.coef = res_lst[[j]]$zero.coef)
          }else if(is.null(res_lst[[j]]$zero.coef)){
            lineage.coef[[k]] = list(pos.coef = res_lst[[j]]$pos.coef)
          }else{
            lineage.coef[[k]] = list(pos.coef = res_lst[[j]]$pos.coef, zero.coef = res_lst[[j]]$zero.coef)
          }
          k = k + 1
        }
      }
    }else{
      k = 1
      for(j in 1:length(tax.info.lst)){
        # print(tax.info.lst[[j]])
        res_tmp = lineage.score.test.meta(tax.info.lst[[j]], W.data.list, W.data.rare.list, otucols, X, X.index, class.id, permute.strata = permute.strata,  n.replicates=n.perm, meta.method = meta.method)
        if(any(is.na(res_tmp))|(is.null(res_tmp))){
          next
        }else{
          subtree = c(subtree, unname(tax.info.lst[[j]]["level.uni"]))
          pval = cbind(pval, res_tmp$pval)
          if(is.null(res_tmp$pos.coef)){
            lineage.coef[[k]] = list(zero.coef = res_tmp$zero.coef)
          }else if(is.null(res_tmp$zero.coef)){
            lineage.coef[[k]] = list(pos.coef = res_tmp$pos.coef)
          }else{
            lineage.coef[[k]] = list(pos.coef = res_tmp$pos.coef, zero.coef = res_tmp$zero.coef)
          }
          k = k + 1
        }
      }
    }

    colnames(pval) = subtree
    if(is.null(n.perm)){
      rownames(pval) = paste0("Asymptotic-",meta.method)
      score.tmp = pval[1,]
    }else{
      if(meta.method %in% c("RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")){
        rownames(pval) = paste0("Resampling-",meta.method)
        score.tmp = pval[1,]
      }else if(meta.method %in% c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O")){
        rownames(pval) = c(paste0("Asymptotic-",meta.method),paste0("Resampling-",meta.method))
        score.tmp = pval[2,]
      }
      #print(pval)
    }

    # identify significant lineages
    subtree.tmp = subtree
    index.na = which(is.na(score.tmp))
    if(length(index.na)>0){
      # drop those lineages which have NA values
      score.tmp = score.tmp[-index.na]
      subtree.tmp = subtree.tmp[-index.na]
    }

    #score.tmp[score.tmp==0] = 1e-4
    m.test = length(score.tmp)

    # Benjamini-Hochberg FDR control
    index.p = order(score.tmp)
    p.sort = sort(score.tmp)
    #fdr.alpha = 0.05

    # change 2022/12/01
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1
    }

    sig.lineage = subtree.tmp[reject==1]

    # return all the p-values as well as significant lineages
    return( list(lineage.pval=pval, lineage.coef = lineage.coef, sig.lineage=sig.lineage) )
  }

}

