
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

arma::uvec within_class_perm_index_map(const arma::uvec &labels){
  const arma::uword n = labels.n_elem;
  std::unordered_map<arma::uword, std::vector<arma::uword>> groups;
  for (arma::uword i = 0; i < n; ++i) {
    groups[labels[i]].push_back(i);
  }
  arma::uvec perm_idx = arma::linspace<arma::uvec>(0, n - 1, n);
  for (auto &kv : groups) {
    std::vector<arma::uword> &rows_vec = kv.second;
    if (rows_vec.size() <= 1) continue;
    arma::uvec rows(rows_vec);

    // Fisher–Yates shuffle
    for (int k = (int)rows.n_elem - 1; k > 0; --k) {
      int j = (int)std::floor(arma::randu() * (k + 1));
      std::swap(rows[k], rows[j]);
    }
    for (std::size_t k = 0; k < rows_vec.size(); ++k) {
      perm_idx[rows_vec[k]] = rows[k];
    }
  }
  return perm_idx;
}

arma::uvec between_class_perm_index_map(const arma::uvec &labels){
  arma::uword n = labels.n_elem;
  std::unordered_map<arma::uword, std::vector<arma::uword>> groups;
  for (arma::uword i = 0; i < n; ++i) {
    groups[labels[i]].push_back(i);
  }
  arma::uvec perm_idx = arma::linspace<arma::uvec>(0, n - 1, n);
  arma::uword ncluster = groups.size();
  arma::uvec cluster_perm = arma::randperm(ncluster);
  for(std::size_t i = 0; i < ncluster; ++i){
    std::vector<arma::uword> &rows_vec = groups[i];
    arma::uword idx = groups[cluster_perm[i]][0];
    for (std::size_t k = 0; k < rows_vec.size(); ++k){
      perm_idx[rows_vec[k]] = idx;
    }
  }
  return perm_idx;
}

arma::mat fun_score_i_beta(const arma::mat &Y_R, const arma::mat &X){
  unsigned int n = Y_R.n_rows;
  unsigned int m = Y_R.n_cols;
  unsigned int p = X.n_cols;
  unsigned int n_beta = m * p;
  arma::mat Score_beta_i(n, n_beta, arma::fill::zeros);
  for(unsigned int i = 0; i < n; ++i){
    Score_beta_i.row(i) = arma::kron(Y_R.row(i), X.row(i));
  }
  return(Score_beta_i);
}

arma::mat fun_hessian_beta(const arma::mat &Y_I, const arma::mat &X){
  unsigned int m = Y_I.n_cols;
  unsigned int p = X.n_cols;
  unsigned int n_beta = m * p;
  arma::mat Hessian_beta(p, n_beta, arma::fill::zeros);
  for(unsigned int i = 0; i < m; ++i){
    arma::mat XYI = X.each_col() % Y_I.col(i);
    Hessian_beta.cols(i*p, (i+1)*p - 1) = XYI.t() * X;
  }
  return(Hessian_beta);
}

// [[Rcpp::export]]
Rcpp::List score_test_stat_meta(const Rcpp::List &X_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, const Rcpp::List& col_index_list,
                                const Rcpp::List &Y_R_list, const Rcpp::List &Y_I_list, const arma::uvec &X_par_index, const unsigned int n_par_interest_beta,
                                const std::string &meta_method){
  // return a vector of score statistics with each element corresponding to one subject.
  // the total study number
  unsigned int total_num = X_list.length();
  double U_theta = 0; double V_theta = 0;
  // initialize for later use
  arma::vec score_stat_beta = arma::vec(n_par_interest_beta, arma::fill::zeros);
  // arma::vec score_stat_beta_vec = arma::vec(total_num);
  arma::mat est_cov_meta(n_par_interest_beta, n_par_interest_beta, arma::fill::zeros);
  Rcpp::List cid;
  bool flag = cluster_list.isNotNull();
  if(flag){cid = cluster_list.get();}
  Rcpp::List cov_lst(total_num); Rcpp::List est_beta_lst(total_num);
  // create null list to save the results for each study
  // Rcpp::List score_beta(total_num); Rcpp::List est_cov(total_num);
  for(unsigned int j = 0; j < total_num; j++){
    arma::mat Y_I = Y_I_list[j];
    arma::mat Y_R = Y_R_list[j];
    arma::mat X = X_list[j];
    arma::uvec cluster_index;
    unsigned int n = Y_I.n_rows;
    unsigned int m = Y_I.n_cols;
    unsigned int p = X.n_cols;
    unsigned int n_par_int = X_par_index.n_elem;
    unsigned int ncate = 0; bool has_cluster = FALSE;
    if(flag){has_cluster = !Rf_isNull(cid[j]);}
    if(has_cluster){
      cluster_index = Rcpp::as<arma::uvec>(cid[j]);
      ncate = arma::max(cluster_index); cluster_index -= 1;
    }else{ has_cluster = 0u; ncate = n; }
    arma::mat Score_reduce_beta = fun_score_i_beta(Y_R, X);
    arma::mat Hessian_beta = fun_hessian_beta(Y_I, X);
    arma::rowvec colsums_Score_beta = arma::sum(Score_reduce_beta, 0);
    arma::uvec X_par_index_u = X_par_index - 1; // change to 0-based index
    arma::mat est_beta(n_par_int, m, arma::fill::zeros);
    arma::uvec mask = arma::ones<arma::uvec>(p);
    mask.elem(X_par_index_u).zeros();
    arma::uvec X_dis_par_index_u = arma::find(mask == 1);
    arma::mat U_mat_combine(n_par_int * m, ncate, arma::fill::zeros);
    arma::mat InvI1_combine(n_par_int * m, n_par_int * m, arma::fill::zeros);
    // check the dimension of beta
    for(unsigned int i = 0; i < m; ++i){
      arma::uvec idx_i = X_par_index_u + (i * p);
      arma::mat Hessian_beta_i = Hessian_beta.cols(i*p, (i+1)*p - 1);
      arma::mat Score_reduce_beta_i = Score_reduce_beta.cols(i*p, (i+1)*p - 1);
      arma::mat I1 = Hessian_beta_i(X_par_index_u, X_par_index_u) - Hessian_beta_i(X_par_index_u, X_dis_par_index_u) *
                     arma::pinv(Hessian_beta_i(X_dis_par_index_u, X_dis_par_index_u)) * Hessian_beta_i(X_dis_par_index_u, X_par_index_u);
      arma::mat InvI1 = arma::pinv(I1);
      est_beta.col(i) = InvI1 * colsums_Score_beta.cols(idx_i).t();
      arma::mat U_mat(n_par_int, ncate, arma::fill::zeros);
      for(unsigned int k = 0; k < n; ++k){
        arma::uvec rowidx = arma::uvec{k};
        if(has_cluster){
          unsigned int cls_id = cluster_index(k);
          U_mat.col(cls_id) += Score_reduce_beta_i(rowidx, X_par_index_u).t() -  Hessian_beta_i(X_par_index_u, X_dis_par_index_u) *
            arma::pinv(Hessian_beta_i(X_dis_par_index_u, X_dis_par_index_u)) * Score_reduce_beta_i(rowidx, X_dis_par_index_u).t();
        }else{
          U_mat.col(k) = Score_reduce_beta_i(rowidx, X_par_index_u).t() -  Hessian_beta_i(X_par_index_u, X_dis_par_index_u) *
            arma::pinv(Hessian_beta_i(X_dis_par_index_u, X_dis_par_index_u)) * Score_reduce_beta_i(rowidx, X_dis_par_index_u).t();
        }
      }
      InvI1_combine.submat(i*n_par_int, i*n_par_int, (i+1)*n_par_int - 1, (i+1)*n_par_int - 1) = InvI1;
      U_mat_combine.rows(i*n_par_int, (i+1)*n_par_int - 1) = U_mat;
    }
    arma::mat cov_mat = InvI1_combine * (U_mat_combine * U_mat_combine.t()) * InvI1_combine;
    cov_lst[j] = cov_mat; est_beta_lst[j] = est_beta;
    // score_beta[j] = est_beta;
    // est_cov[j] = cov_mat;
    if(meta_method == "RE-MetaQCAT"){
      arma::mat est_inv = arma::pinv(cov_mat);
      arma::vec tmp = est_inv * arma::vectorise(est_beta);
      U_theta = U_theta + 0.5 * arma::as_scalar(tmp.t() * tmp - arma::trace(est_inv));
      V_theta = V_theta + 0.5 * arma::as_scalar(arma::trace(est_inv.t() * est_inv));
    }
    arma::uvec col_index_u  = col_index_list[j]; col_index_u -= 1;
    est_cov_meta.submat(col_index_u,col_index_u)  += arma::pinv(cov_mat);
    score_stat_beta(col_index_u) = score_stat_beta(col_index_u) + arma::pinv(cov_mat) * arma::vectorise(est_beta);
  }
  arma::mat cov_beta_hat = arma::pinv(est_cov_meta);
  arma::vec beta_hat = cov_beta_hat * score_stat_beta;
  arma::vec var_beta_hat = cov_beta_hat.diag();
  double score_stat_meta = arma::as_scalar(score_stat_beta.t() * cov_beta_hat * score_stat_beta);
  if(meta_method == "RE-MetaQCAT"){
    if (!std::isfinite(V_theta) || V_theta <= 0) Rcpp::stop("Degenerate random-effect covariance.");
    score_stat_meta = score_stat_meta + U_theta * U_theta / V_theta;
    return Rcpp::List::create(Rcpp::Named("score.stat.meta") = score_stat_meta,
                              Rcpp::Named("var.beta.hat") = var_beta_hat,
                              Rcpp::Named("beta.hat") = beta_hat,
                              Rcpp::Named("est.beta.lst") = est_beta_lst,
                              Rcpp::Named("cov.lst") = cov_lst);
  }else{
    double score_pvalue = R::pchisq(score_stat_meta, n_par_interest_beta, 0, 0);
    return Rcpp::List::create(Rcpp::Named("score.stat.meta") = score_stat_meta,
                              Rcpp::Named("score.pvalue") = score_pvalue,
                              Rcpp::Named("var.beta.hat") = var_beta_hat,
                              Rcpp::Named("beta.hat") = beta_hat,
                              Rcpp::Named("est.beta.lst") = est_beta_lst,
                              Rcpp::Named("cov.lst") = cov_lst);
  }
}


double score_test_stat_meta_resampling(const Rcpp::List &X_perm_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, const Rcpp::List& col_index_list,
                                       const Rcpp::List &Y_R_list, const Rcpp::List& Y_I_list, const arma::uvec &X_par_index, const unsigned int n_par_interest_beta,
                                       const std::string &meta_method){
  // the total study number
  int total_num = X_perm_list.length();
  double U_theta = 0; double V_theta = 0;
  // initialize for later use
  arma::vec score_stat_beta = arma::vec(n_par_interest_beta, arma::fill::zeros);
  // arma::vec score_stat_beta_vec = arma::vec(total_num);
  arma::mat est_cov_meta(n_par_interest_beta,n_par_interest_beta, arma::fill::zeros);
  Rcpp::List cid; bool flag = cluster_list.isNotNull();
  if(flag){cid = cluster_list.get();}
  for(int j = 0; j < total_num; j++){
    arma::mat Y_I = Y_I_list[j];
    arma::mat Y_R = Y_R_list[j];
    arma::mat X_perm = X_perm_list[j];
    arma::uvec cluster_index;
    unsigned int n = Y_I.n_rows;
    unsigned int m = Y_I.n_cols;
    unsigned int p = X_perm.n_cols;
    unsigned int n_par_int = X_par_index.n_elem;
    unsigned int ncate = 0; bool has_cluster = FALSE;
    if(flag){has_cluster = !Rf_isNull(cid[j]);}
    if(has_cluster){
      cluster_index = Rcpp::as<arma::uvec>(cid[j]);
      ncate = arma::max(cluster_index); cluster_index -= 1;
    }else{has_cluster = 0u; ncate = n; }
    arma::mat Score_reduce_beta = fun_score_i_beta(Y_R, X_perm);
    arma::mat Hessian_beta = fun_hessian_beta(Y_I, X_perm);
    arma::rowvec colsums_Score_beta = arma::sum(Score_reduce_beta, 0);
    arma::uvec X_par_index_u = X_par_index - 1; // change to 0-based index
    arma::mat est_beta(n_par_int, m, arma::fill::zeros);
    arma::uvec mask = arma::ones<arma::uvec>(p);
    mask.elem(X_par_index_u).zeros();
    arma::uvec X_dis_par_index_u = arma::find(mask == 1);
    arma::mat U_mat_combine(n_par_int * m, ncate, arma::fill::zeros);
    arma::mat InvI1_combine(n_par_int * m, n_par_int * m, arma::fill::zeros);
    // check the dimension of beta
    for(unsigned int i = 0; i < m; ++i){
      arma::uvec idx_i = X_par_index_u + (i * p);
      arma::mat Hessian_beta_i = Hessian_beta.cols(i*p, (i+1)*p - 1);
      arma::mat Score_reduce_beta_i = Score_reduce_beta.cols(i*p, (i+1)*p - 1);
      arma::mat I1 = Hessian_beta_i(X_par_index_u, X_par_index_u) - Hessian_beta_i(X_par_index_u, X_dis_par_index_u) *
                     arma::pinv(Hessian_beta_i(X_dis_par_index_u, X_dis_par_index_u)) * Hessian_beta_i(X_dis_par_index_u, X_par_index_u);
      arma::mat InvI1 = arma::pinv(I1);
      est_beta.col(i) = InvI1 * colsums_Score_beta.cols(idx_i).t();
      arma::mat U_mat(n_par_int, ncate, arma::fill::zeros);
      for(unsigned int k = 0; k < n; ++k){
        arma::uvec rowidx = arma::uvec{k};
        if(has_cluster){
          unsigned int cls_id = cluster_index(k);
          U_mat.col(cls_id) += Score_reduce_beta_i(rowidx, X_par_index_u).t() -  Hessian_beta_i(X_par_index_u, X_dis_par_index_u) *
            arma::pinv(Hessian_beta_i(X_dis_par_index_u, X_dis_par_index_u)) * Score_reduce_beta_i(rowidx, X_dis_par_index_u).t();
        }else{
          U_mat.col(k) = Score_reduce_beta_i(rowidx, X_par_index_u).t() -  Hessian_beta_i(X_par_index_u, X_dis_par_index_u) *
            arma::pinv(Hessian_beta_i(X_dis_par_index_u, X_dis_par_index_u)) * Score_reduce_beta_i(rowidx, X_dis_par_index_u).t();
        }
      }
      InvI1_combine.submat(i*n_par_int, i*n_par_int, (i+1)*n_par_int - 1, (i+1)*n_par_int - 1) = InvI1;
      U_mat_combine.rows(i*n_par_int, (i+1)*n_par_int - 1) = U_mat;
    }
    arma::mat cov_mat = InvI1_combine * (U_mat_combine * U_mat_combine.t()) * InvI1_combine;
    if(meta_method == "RE-MetaQCAT"){
      arma::mat est_inv = arma::pinv(cov_mat);
      arma::vec tmp = est_inv * arma::vectorise(est_beta);
      U_theta = U_theta + 0.5 * arma::as_scalar(tmp.t() * tmp - arma::trace(est_inv));
      V_theta = V_theta + 0.5 * arma::as_scalar(arma::trace(est_inv.t() * est_inv));
    }
    arma::uvec col_index_u = col_index_list[j]; col_index_u -= 1;
    est_cov_meta.submat(col_index_u,col_index_u) += arma::pinv(cov_mat);
    score_stat_beta(col_index_u) = score_stat_beta(col_index_u) + arma::pinv(cov_mat) * arma::vectorise(est_beta);
  }
  double score_stat_meta_perm = arma::as_scalar(score_stat_beta.t() * arma::pinv(est_cov_meta) * score_stat_beta);
  if(meta_method == "RE-MetaQCAT"){
    if (!std::isfinite(V_theta) || V_theta <= 0) Rcpp::stop("Degenerate random-effect covariance.");
    score_stat_meta_perm = score_stat_meta_perm + U_theta * U_theta / V_theta;
  }    // for fixed-effect MetaQCAT, the score statistics and estimated covariance matrix are calculated as above}
  return (score_stat_meta_perm);
}

// [[Rcpp::export]]
double resample_pvalue(const Rcpp::List &X_list, const Rcpp::Nullable<Rcpp::List>& cluster_id_list, const Rcpp::List& col_index_list, const Rcpp::List &Y_R_list,
                       const Rcpp::List& Y_I_list, const arma::uvec &X_par_index, const double score_stat_meta, const int n_par_interest_beta, const unsigned int &n_replicates,
                       const std::string &meta_method, const Rcpp::Nullable<Rcpp::String> &permute_strata){
  // arma::arma_rng::set_seed(seed);
  std::string strata = permute_strata.isNotNull() ? Rcpp::as<std::string>(permute_strata) : "";
  Rcpp::List cid; bool class_flag = cluster_id_list.isNotNull();
  if(class_flag){cid = cluster_id_list.get();}
  unsigned int total_num = X_list.length();
  arma::uvec X_par_index_u = X_par_index - 1;
  unsigned int n_one = 0; unsigned int one_acc = 0;
  unsigned int start_nperm = 1;
  unsigned int end_nperm = (n_replicates < 100)? n_replicates : 100;
  bool flag = true;
  while(flag && (start_nperm <= end_nperm) && (end_nperm <= n_replicates)){
    for(unsigned int k = start_nperm; k <= end_nperm; ++k){
      Rcpp::List X_perm_list(total_num);
      if((!class_flag)||((strata != "within")&(strata != "between"))){
        for(unsigned int i = 0; i < total_num; ++i){
          arma::mat X_perm = X_list[i];
          unsigned int n = X_perm.n_rows;
          arma::uvec perm_idx = arma::randperm(n);
          // std::vector<unsigned int> v = arma::conv_to< std::vector<unsigned int> >::from(perm_idx);
          X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
          X_perm_list[i] = X_perm;
        }
      }else{
        if(strata == "within"){
          for(unsigned int i = 0; i < total_num; ++i){
            bool cluster_flag = !Rf_isNull(cid[i]);
            arma::mat X_perm = X_list[i];
            unsigned int n = X_perm.n_rows;
            if(cluster_flag){
              arma::uvec label_vec = Rcpp::as<arma::uvec>(cid[i]); label_vec -= 1;
              arma::uvec perm_idx = within_class_perm_index_map(label_vec);
              // std::vector<unsigned int> v = arma::conv_to< std::vector<unsigned int> >::from(perm_idx);
              // std::vector<unsigned int> v1 = arma::conv_to< std::vector<unsigned int> >::from(label_vec);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }else{
              arma::uvec perm_idx = arma::randperm(n);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }
            X_perm_list[i] = X_perm;
          }
        }else if (strata == "between"){
          for(unsigned int i = 0; i < total_num; ++i){
            bool cluster_flag = !Rf_isNull(cid[i]);
            arma::mat X_perm = X_list[i];
            unsigned int n = X_perm.n_rows;
            if(cluster_flag){
              arma::uvec label_vec = Rcpp::as<arma::uvec>(cid[i]); label_vec -= 1;
              arma::uvec perm_idx = between_class_perm_index_map(label_vec);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }else{
              arma::uvec perm_idx = arma::randperm(n);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }
            X_perm_list[i] = X_perm;
          }
        }
      }
      double score_stat_meta_perm = score_test_stat_meta_resampling(X_perm_list, cluster_id_list, col_index_list, Y_R_list, Y_I_list,
                                                                    X_par_index, n_par_interest_beta, meta_method);
      if (!std::isfinite(score_stat_meta_perm)) Rcpp::stop("Non-finite permutation statistic; p-value not computed.");
      n_one += 1;
      if(score_stat_meta_perm >= score_stat_meta){one_acc += 1;}
    }
    if(one_acc < 1){
      start_nperm = end_nperm + 1;
      end_nperm = std::min(n_replicates, (end_nperm + 1) * 100 - 1);
      flag = true;
    }else if(one_acc < 10){
      start_nperm = end_nperm + 1;
      end_nperm = std::min(n_replicates, (end_nperm + 1) * 10 - 1);
      flag = true;
    }else{flag = false;}
  }
  if (n_one == 0) Rcpp::stop("No successful permutations; p-value not computed.");
  double score_Rpvalue = (one_acc + 1.0)/(n_one + 1.0);
  return(score_Rpvalue);
}



arma::vec Ei_beta(const arma::mat &X_i, const arma::vec &beta, const unsigned int m, const unsigned int p){
  arma::vec Ei_out(m, arma::fill::zeros);
  for(unsigned int j = 0; j < m-1; ++j){
    arma::vec beta_j = beta.subvec(j * p, (j + 1) * p - 1);
    double eta = arma::dot(beta_j, X_i);
    Ei_out(j) = std::exp(eta);
  }
  Ei_out(m-1) = 1.0;
  return(Ei_out);
}

void score_summary_beta(const arma::vec &coef_vec, const arma::mat &Y, const arma::mat &X,
                         arma::mat &Score_reduce_beta, arma::mat &Hessian_beta){
  unsigned int n = Y.n_rows;
  unsigned int m = Y.n_cols;
  unsigned int p = X.n_cols;
  arma::vec nY = arma::sum(Y, 1);
  unsigned int n_beta = (m-1) * p;
  arma::vec Score_beta = arma::vec(n_beta, arma::fill::zeros);
  for(unsigned int i = 0; i < n; ++i){
    arma::vec Ei_vec = Ei_beta(X.row(i).t(), coef_vec, m, p); // m x 1
    double sum_Ei = arma::sum(Ei_vec);
    arma::vec Pi_vec = Ei_vec / sum_Ei; 
    arma::uvec idx = arma::regspace<arma::uvec>(0, m-1);
    idx.shed_row(m-1); 
    arma::vec Pi_tmp = Pi_vec.elem(idx);
    arma::vec a = Y.row(i).t();  // m x 1
    arma::vec s = a.elem(idx) - nY(i) * Pi_tmp;  // (m-1) x 1
    Score_reduce_beta.row(i) = arma::kron(s, X.row(i).t()).t();
    arma::mat tmp_beta = Pi_tmp * Pi_tmp.t() - arma::diagmat(Pi_tmp);
    tmp_beta = tmp_beta * nY(i);
    Hessian_beta -= arma::kron(tmp_beta, X.row(i).t() * X.row(i));
  }
}

// [[Rcpp::export]]
Rcpp::List score_test_stat_QCAT_meta(const Rcpp::List &X_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, const Rcpp::List& col_index_list, 
                                     const Rcpp::List &Y_list, const Rcpp::List &coef_List, const arma::vec &X_par_index, const unsigned int n_par_interest_beta,
                                     const std::string &meta_method){
  // return a vector of score statistics with each element corresponding to one subject.
  // the total study number
  unsigned int total_num = X_list.length();
  double U_theta = 0; double V_theta = 0;
  // initialize for later use
  arma::vec score_stat_beta = arma::vec(n_par_interest_beta, arma::fill::zeros);
  // arma::vec score_stat_beta_vec = arma::vec(total_num);
  arma::mat est_cov_meta(n_par_interest_beta, n_par_interest_beta, arma::fill::zeros);
  Rcpp::List cid;
  bool flag = cluster_list.isNotNull();
  if(flag){cid = cluster_list.get();}
  Rcpp::List cov_lst(total_num); Rcpp::List est_beta_lst(total_num);
  // create null list to save the results for each study
  // Rcpp::List score_beta(total_num); Rcpp::List est_cov(total_num);
  for(unsigned int j = 0; j < total_num; ++j){
    arma::mat Y = Y_list[j];
    arma::mat X = X_list[j];
    unsigned int p = X.n_cols; 
    unsigned int n = X.n_rows;
    unsigned int m = Y.n_cols;
    unsigned int n_beta = (m-1) * p;
    // arma::mat vA_mat(m, n, arma::fill::zeros);
    // arma::cube Vinv_cube(m, m, n, arma::fill::zeros);
    // arma::mat VY_mat(m, n, arma::fill::zeros);
    arma::mat Score_reduce_beta(n, n_beta, arma::fill::zeros);
    arma::mat Hessian_beta(n_beta, n_beta, arma::fill::zeros);
    arma::uvec cluster_index;
    unsigned int ncate = 0; bool has_cluster = FALSE;
    if(flag){has_cluster = !Rf_isNull(cid[j]);}
    if(has_cluster){
      cluster_index = Rcpp::as<arma::uvec>(cid[j]);
      ncate = arma::max(cluster_index); cluster_index -= 1;
    }else{ has_cluster = 0u; ncate = n; }
    // arma::uvec Z_par_index_u = Z_par_index - 1;
    arma::mat idx = arma::kron(arma::linspace(0,(m-2) * p, (m-1)), arma::ones(X_par_index.n_elem)) + arma::kron(arma::ones(m-1), X_par_index);
    arma::uvec uidx = arma::conv_to< arma::uvec>::from(arma::vectorise(idx)-1);
    arma::uvec mask = arma::ones<arma::uvec>(n_beta);
    mask.elem(uidx).zeros();
    arma::uvec compl_uidx = arma::find(mask==1); 

    arma::vec coef_vec = coef_List[j]; 
    arma::vec coef_vec_full = arma::zeros(n_beta);
    coef_vec_full(compl_uidx) = coef_vec;
    // score_summary_alpha(coef_vec_full, Y_bin, Z, vA_mat, Vinv_cube,
    //                     VY_mat, Score_reduce_alpha, Hessian_alpha);
    score_summary_beta(coef_vec_full, Y, X, Score_reduce_beta, Hessian_beta); 
    // vA_lst[j] = vA_mat; Vinv_lst[j] = Vinv_cube; VY_lst[j] = VY_mat;
    arma::mat Score_reduce_beta_reorg = arma::join_rows(Score_reduce_beta.cols(uidx), Score_reduce_beta.cols(compl_uidx));
    arma::mat Hess_reduce_beta_reorg = arma::join_cols(arma::join_rows(Hessian_beta.submat(uidx, uidx), 
                                                        Hessian_beta.submat(uidx, compl_uidx)),
                                                        arma::join_rows(Hessian_beta.submat(compl_uidx, uidx), 
                                                        Hessian_beta.submat(compl_uidx, compl_uidx)));

    unsigned int n_par_interest = uidx.n_elem;
    arma::span par_idx(0,(n_par_interest-1)); 
    arma::span compl_par_idx(n_par_interest,(n_beta-1));
    arma::rowvec Score_reduce_rowsums = arma::sum(Score_reduce_beta_reorg, 0);
    arma::vec A = Score_reduce_rowsums.cols(par_idx).t();
    arma::mat B1 = arma::pinv(Hess_reduce_beta_reorg.submat(par_idx, par_idx) - Hess_reduce_beta_reorg.submat(par_idx, compl_par_idx) *
      arma::pinv(Hess_reduce_beta_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_beta_reorg.submat(compl_par_idx, par_idx));
    arma::vec beta_hat = B1 * A;
    arma::mat U = Score_reduce_beta_reorg.cols(par_idx) - Score_reduce_beta_reorg.cols(compl_par_idx) *
      arma::pinv(Hess_reduce_beta_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_beta_reorg.submat(par_idx, compl_par_idx).t();
    arma::mat B2(n_par_interest, n_par_interest, arma::fill::zeros);
    if(has_cluster){
      arma::mat within_class_sum(ncate, n_par_interest, arma::fill::zeros);
      for(int j = 0; j < n; j++){
        unsigned int obs_cls_id = cluster_index(j);
        arma::rowvec tmp_vec2 = arma::conv_to< arma::rowvec >::from(U.row(j));
        within_class_sum.row(obs_cls_id) += tmp_vec2;
      }
      B2 = within_class_sum.t() * within_class_sum;
    }else{ B2 = U.t() * U ;}
    arma::mat cov_beta = B1 * B2 * B1;
    est_beta_lst[j] = beta_hat; // restore the score statistics and estimated covariance matrix in lists which are of necessity for some meta-analysis methods
    cov_lst[j] = cov_beta;
    arma::uvec col_index_u  = col_index_list[j]; col_index_u -= 1;
    est_cov_meta.submat(col_index_u,col_index_u) += arma::pinv(cov_beta);
    score_stat_beta(col_index_u) = score_stat_beta(col_index_u) + arma::pinv(cov_beta) * beta_hat;
    if(meta_method == "RE-MetaQCAT"){
      arma::mat est_inv = arma::pinv(cov_beta);
      arma::vec tmp = est_inv * arma::vectorise(beta_hat);
      U_theta = U_theta + 0.5 * arma::as_scalar(tmp.t() * tmp - arma::trace(est_inv));
      V_theta = V_theta + 0.5 * arma::as_scalar(arma::trace(est_inv.t() * est_inv));
    }
  }
  arma::mat cov_beta_hat = arma::pinv(est_cov_meta);
  arma::vec beta_hat = cov_beta_hat * score_stat_beta;
  arma::vec var_beta_hat = cov_beta_hat.diag();
  double score_stat_meta = arma::as_scalar(score_stat_beta.t() * cov_beta_hat * score_stat_beta);
  if(meta_method == "RE-MetaQCAT"){
    if (!std::isfinite(V_theta) || V_theta <= 0) Rcpp::stop("Degenerate random-effect covariance.");
    score_stat_meta = score_stat_meta + U_theta * U_theta / V_theta;
    return Rcpp::List::create(Rcpp::Named("score.stat.meta") = score_stat_meta,
                              Rcpp::Named("var.beta.hat") = var_beta_hat,
                              Rcpp::Named("beta.hat") = beta_hat,
                              Rcpp::Named("est.beta.lst") = est_beta_lst,
                              Rcpp::Named("cov.lst") = cov_lst);
  }else{
    double score_pvalue = R::pchisq(score_stat_meta, n_par_interest_beta, 0, 0);
    return Rcpp::List::create(Rcpp::Named("score.stat.meta") = score_stat_meta,
                              Rcpp::Named("score.pvalue") = score_pvalue,
                              Rcpp::Named("var.beta.hat") = var_beta_hat,
                              Rcpp::Named("beta.hat") = beta_hat,
                              Rcpp::Named("est.beta.lst") = est_beta_lst,
                              Rcpp::Named("cov.lst") = cov_lst);
  }
}


double score_test_stat_QCAT_meta_resampling(const Rcpp::List& X_perm_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, 
                                            const Rcpp::List& col_index_list, const Rcpp::List &Y_list, const Rcpp::List &coef_List,
                                            const arma::vec& X_par_index, const unsigned int n_par_interest_beta, const std::string &meta_method){
  
  //alpha.meta.results$vA.lst, alpha.meta.results$Vinv.lst, alpha.meta.results$VY.lst,
  // the total study number
  int total_num = X_perm_list.length();
  // int total_num_tmp = coef_List.length();
  // int total_num_ttt = Z_perm_list.length();
  double U_theta = 0; double V_theta = 0;
  arma::vec score_stat_beta = arma::vec(n_par_interest_beta,arma::fill::zeros);
  arma::mat est_cov_meta(n_par_interest_beta,n_par_interest_beta);
  est_cov_meta.zeros();
  Rcpp::List cid; bool flag = cluster_list.isNotNull();
  if(flag){cid = cluster_list.get();}   
  // create null list to save the results for each study
  Rcpp::List score_alpha(total_num);
  Rcpp::List est_cov(total_num);
  for(unsigned int j = 0; j < total_num; ++j){
    arma::mat Y = Y_list[j];
    arma::mat X_perm = X_perm_list[j];
    unsigned int p = X_perm.n_cols; 
    unsigned int n = X_perm.n_rows;
    unsigned int m = Y.n_cols;
    unsigned int n_beta = (m-1) * p;
    // arma::mat vA_mat(m, n, arma::fill::zeros);
    // arma::cube Vinv_cube(m, m, n, arma::fill::zeros);
    // arma::mat VY_mat(m, n, arma::fill::zeros);
    arma::mat Score_reduce_beta(n, n_beta, arma::fill::zeros);
    arma::mat Hessian_beta(n_beta, n_beta, arma::fill::zeros);
    arma::uvec cluster_index;
    unsigned int ncate = 0; bool has_cluster = FALSE;
    if(flag){has_cluster = !Rf_isNull(cid[j]);}
    if(has_cluster){
      cluster_index = Rcpp::as<arma::uvec>(cid[j]);
      ncate = arma::max(cluster_index); cluster_index -= 1;
    }else{ has_cluster = 0u; ncate = n; }
    // arma::uvec Z_par_index_u = Z_par_index - 1;
    arma::mat idx = arma::kron(arma::linspace(0,(m-2) * p, (m-1)), arma::ones(X_par_index.n_elem)) + arma::kron(arma::ones(m-1), X_par_index);
    arma::uvec uidx = arma::conv_to< arma::uvec>::from(arma::vectorise(idx)-1);
    arma::uvec mask = arma::ones<arma::uvec>(n_beta);
    mask.elem(uidx).zeros();
    arma::uvec compl_uidx = arma::find(mask==1); 
    arma::vec coef_vec = coef_List[j]; 
    arma::vec coef_vec_full = arma::zeros(n_beta);
    coef_vec_full(compl_uidx) = coef_vec;
    // score_summary_alpha(coef_vec_full, Y_bin, Z, vA_mat, Vinv_cube,
    //                     VY_mat, Score_reduce_alpha, Hessian_alpha);
    score_summary_beta(coef_vec_full, Y, X_perm, Score_reduce_beta, Hessian_beta); 
    // vA_lst[j] = vA_mat; Vinv_lst[j] = Vinv_cube; VY_lst[j] = VY_mat;
    arma::mat Score_reduce_beta_reorg = arma::join_rows(Score_reduce_beta.cols(uidx), Score_reduce_beta.cols(compl_uidx));
    arma::mat Hess_reduce_beta_reorg = arma::join_cols(arma::join_rows(Hessian_beta.submat(uidx, uidx), 
                                                        Hessian_beta.submat(uidx, compl_uidx)),
                                                        arma::join_rows(Hessian_beta.submat(compl_uidx, uidx), 
                                                        Hessian_beta.submat(compl_uidx, compl_uidx)));

    unsigned int n_par_interest = uidx.n_elem;
    arma::span par_idx(0,(n_par_interest-1)); 
    arma::span compl_par_idx(n_par_interest,(n_beta-1));
    arma::rowvec Score_reduce_rowsums = arma::sum(Score_reduce_beta_reorg, 0);
    arma::vec A = Score_reduce_rowsums.cols(par_idx).t();
    arma::mat B1 = arma::pinv(Hess_reduce_beta_reorg.submat(par_idx, par_idx) - Hess_reduce_beta_reorg.submat(par_idx, compl_par_idx) *
      arma::pinv(Hess_reduce_beta_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_beta_reorg.submat(compl_par_idx, par_idx));
    arma::vec beta_hat = B1 * A;
    arma::mat U = Score_reduce_beta_reorg.cols(par_idx) - Score_reduce_beta_reorg.cols(compl_par_idx) *
      arma::pinv(Hess_reduce_beta_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_beta_reorg.submat(par_idx, compl_par_idx).t();
    arma::mat B2(n_par_interest, n_par_interest, arma::fill::zeros);
    if(has_cluster){
      arma::mat within_class_sum(ncate, n_par_interest, arma::fill::zeros);
      for(int j = 0; j < n; j++){
        unsigned int obs_cls_id = cluster_index(j);
        arma::rowvec tmp_vec2 = arma::conv_to< arma::rowvec >::from(U.row(j));
        within_class_sum.row(obs_cls_id) += tmp_vec2;
      }
      B2 = within_class_sum.t() * within_class_sum;
    }else{ B2 = U.t() * U ;}
    arma::mat cov_beta = B1 * B2 * B1;
    arma::uvec col_index_u  = col_index_list[j]; col_index_u -= 1;
    est_cov_meta.submat(col_index_u,col_index_u) += arma::pinv(cov_beta);
    score_stat_beta(col_index_u) = score_stat_beta(col_index_u) + arma::pinv(cov_beta) * beta_hat;
    if(meta_method == "RE-MetaQCAT"){
      arma::mat est_inv = arma::pinv(cov_beta);
      arma::vec tmp = est_inv * arma::vectorise(beta_hat);
      U_theta = U_theta + 0.5 * arma::as_scalar(tmp.t() * tmp - arma::trace(est_inv));
      V_theta = V_theta + 0.5 * arma::as_scalar(arma::trace(est_inv.t() * est_inv));
    }
  }
  double score_stat_meta = arma::as_scalar(score_stat_beta.t() * arma::pinv(est_cov_meta) * score_stat_beta);
  if(meta_method == "RE-MetaQCAT"){
    if (!std::isfinite(V_theta) || V_theta <= 0) Rcpp::stop("Degenerate random-effect covariance.");
    score_stat_meta = score_stat_meta + U_theta * U_theta / V_theta;
  }
  return(score_stat_meta);
}

// [[Rcpp::export]]
double resample_QCAT_pvalue(const Rcpp::List& X_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, const Rcpp::List& col_index_list, 
                            const Rcpp::List &Y_list, const Rcpp::List& coef_List, const arma::vec& X_par_index, const double score_stat_meta, 
                            const unsigned int n_par_interest_beta, const unsigned int n_perm, const std::string &meta_method, 
                            const Rcpp::Nullable<Rcpp::String> &permute_strata){
  std::string strata = permute_strata.isNotNull() ? Rcpp::as<std::string>(permute_strata) : "";
  Rcpp::List cid; bool class_flag = cluster_list.isNotNull();
  if(class_flag){cid = cluster_list.get();}
  unsigned int total_num = X_list.length();
  arma::uvec X_par_index_u = arma::conv_to< arma::uvec >::from(X_par_index - 1);
  unsigned int n_one = 0; unsigned int one_acc = 0;
  unsigned int start_nperm = 1;
  unsigned int end_nperm = (n_perm < 100)? n_perm : 100;
  bool flag = true;
  while(flag && (start_nperm <= end_nperm) && (end_nperm <= n_perm)){
    for(unsigned int k = start_nperm; k <= end_nperm; ++k){
      Rcpp::List X_perm_list(total_num);
      if((!class_flag)||((strata != "within")&(strata != "between"))){
        for(unsigned int i = 0; i < total_num; ++i){
          arma::mat X_perm = X_list[i];
          unsigned int n = X_perm.n_rows;
          arma::uvec perm_idx = arma::randperm(n);
          // std::vector<unsigned int> v = arma::conv_to< std::vector<unsigned int> >::from(perm_idx);
          X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
          X_perm_list[i] = X_perm;
        }
      }else{
        if(strata == "within"){
          for(unsigned int i = 0; i < total_num; ++i){
            bool has_cluster = !Rf_isNull(cid[i]);
            arma::mat X_perm = X_list[i];
            unsigned int n = X_perm.n_rows;
            if(has_cluster){
              arma::uvec label_vec = Rcpp::as<arma::uvec>(cid[i]); label_vec -= 1;
              arma::uvec perm_idx = within_class_perm_index_map(label_vec);
              // std::vector<unsigned int> v = arma::conv_to< std::vector<unsigned int> >::from(perm_idx);
              // std::vector<unsigned int> v1 = arma::conv_to< std::vector<unsigned int> >::from(label_vec);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }else{
              arma::uvec perm_idx = arma::randperm(n);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }
            X_perm_list[i] = X_perm;
          }
        }else if (strata == "between"){
          for(unsigned int i = 0; i < total_num; ++i){
            bool has_cluster = !Rf_isNull(cid[i]);
            arma::mat X_perm = X_list[i];
            unsigned int n = X_perm.n_rows;
            if(has_cluster){
              arma::uvec label_vec = Rcpp::as<arma::uvec>(cid[i]); label_vec -= 1;
              arma::uvec perm_idx = between_class_perm_index_map(label_vec);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }else{
              arma::uvec perm_idx = arma::randperm(n);
              X_perm.cols(X_par_index_u) = X_perm.submat(perm_idx, X_par_index_u);
            }
            X_perm_list[i] = X_perm;
          }
        }
      }
      double score_stat_meta_perm = score_test_stat_QCAT_meta_resampling(X_perm_list, cluster_list, col_index_list, Y_list, 
                                                                         coef_List, X_par_index, n_par_interest_beta, meta_method);
      if (!std::isfinite(score_stat_meta_perm)) Rcpp::stop("Non-finite permutation statistic; p-value not computed.");
      n_one += 1;
      if(score_stat_meta_perm >= score_stat_meta){one_acc += 1;}
    }
    if(one_acc < 1){
      start_nperm = end_nperm + 1;
      end_nperm = std::min(n_perm, (end_nperm + 1) * 100 - 1);
      flag = true;
    }else if(one_acc < 10){
      start_nperm = end_nperm + 1;
      end_nperm = std::min(n_perm, (end_nperm + 1) * 10 - 1);
      flag = true;
    }else{flag = false;}
  }
  if (n_one == 0) Rcpp::stop("No successful permutations; p-value not computed.");
  double score_Rpvalue = (one_acc + 1.0)/(n_one + 1.0);
  return(score_Rpvalue);
}


arma::vec Pi_alpha(const arma::mat &X_i, const arma::vec &alpha, const unsigned int m, const unsigned int p){
  arma::vec Pi_out(m, arma::fill::zeros);
  for(unsigned int j = 0; j < m; ++j){
    arma::vec alpha_j = alpha.subvec(j * p, (j + 1) * p - 1);
    double eta = arma::dot(alpha_j, X_i);
    double tmp = std::exp(eta);
    if (!std::isfinite(tmp)) {
      Pi_out(j) = 1.0;
    } else {
      Pi_out(j) = tmp / (1.0 + tmp);
    }
  }
  return(Pi_out);
}

// void score_summary_alpha(const arma::vec &coef_vec, const arma::mat &Y_bin, const arma::mat &Z,
//                          arma::mat &vA_mat, arma::cube &Vinv_cube, arma::mat &VY_mat,
//                          arma::mat &Score_reduce_alpha, arma::mat &Hessian_alpha){
//   unsigned int n = Y_bin.n_rows;
//   unsigned int m = Y_bin.n_cols;
//   unsigned int p = Z.n_cols;
//   for(unsigned int i = 0; i < n; ++i){
//     arma::vec Pi_vec = Pi_alpha(Z.row(i).t(), coef_vec, m, p); // m x 1
//     arma::vec vA_tmp = Pi_vec % (1 - Pi_vec);
//     arma::mat A_i = arma::diagmat(vA_tmp); // m x m
//     arma::mat t_D_i = arma::kron(A_i, Z.row(i).t()); // n_alpha x m
//     arma::mat A_i_inv = arma::pinv(A_i);
//     arma::vec VY = A_i_inv * (Y_bin.row(i).t() - Pi_vec); // m x 1
//     Score_reduce_alpha.row(i) = (t_D_i * VY).t();
//     Hessian_alpha += t_D_i * A_i_inv * t_D_i.t();
//     vA_mat.col(i) = vA_tmp;
//     Vinv_cube.slice(i) = A_i_inv;
//     VY_mat.col(i) = VY;
//   }
// }


void score_summary_alpha(const arma::vec &coef_vec, const arma::mat &Y_bin, const arma::mat &Z,
                         arma::mat &Score_reduce_alpha, arma::mat &Hessian_alpha){
  unsigned int n = Y_bin.n_rows;
  unsigned int m = Y_bin.n_cols;
  unsigned int p = Z.n_cols;
  for(unsigned int i = 0; i < n; ++i){
    arma::vec Pi_vec = Pi_alpha(Z.row(i).t(), coef_vec, m, p); // m x 1
    arma::vec vA_tmp = Pi_vec % (1 - Pi_vec);
    arma::mat A_i = arma::diagmat(vA_tmp); // m x m
    arma::mat t_D_i = arma::kron(A_i, Z.row(i).t()); // n_alpha x m
    arma::mat A_i_inv = arma::pinv(A_i);
    arma::vec VY = A_i_inv * (Y_bin.row(i).t() - Pi_vec); // m x 1
    Score_reduce_alpha.row(i) = (t_D_i * VY).t();
    Hessian_alpha += t_D_i * A_i_inv * t_D_i.t();
  }
}


// [[Rcpp::export]]
Rcpp::List score_test_stat_zero_meta(const Rcpp::List &Z_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, 
                                     const Rcpp::List& col_index_list, const Rcpp::List &Y_bin_list, const Rcpp::List &coef_List, 
                                     const arma::vec &Z_par_index, const unsigned int n_par_interest_alpha,
                                     const std::string &meta_method){
  // return a vector of score statistics with each element corresponding to one subject.
  // the total study number
  unsigned int total_num = Z_list.length();
  double U_theta = 0; double V_theta = 0;
  arma::vec score_stat_alpha = arma::vec(n_par_interest_alpha, arma::fill::zeros);
  arma::mat est_cov_meta(n_par_interest_alpha, n_par_interest_alpha, arma::fill::zeros);
  Rcpp::List cid;
  bool flag = cluster_list.isNotNull();
  if(flag){cid = cluster_list.get();}
  Rcpp::List cov_lst(total_num); 
  Rcpp::List est_alpha_lst(total_num);   
  // Rcpp::List vA_lst(total_num); 
  // Rcpp::List Vinv_lst(total_num); 
  // Rcpp::List VY_lst(total_num); 
  for(unsigned int j = 0; j < total_num; ++j){
    arma::mat Y_bin = Y_bin_list[j];
    arma::mat Z = Z_list[j];
    unsigned int p = Z.n_cols; 
    unsigned int n = Z.n_rows;
    unsigned int m = Y_bin.n_cols;
    unsigned int n_alpha = m * p;
    // arma::mat vA_mat(m, n, arma::fill::zeros);
    // arma::cube Vinv_cube(m, m, n, arma::fill::zeros);
    // arma::mat VY_mat(m, n, arma::fill::zeros);
    arma::mat Score_reduce_alpha(n, n_alpha, arma::fill::zeros);
    arma::mat Hessian_alpha(n_alpha, n_alpha, arma::fill::zeros);
    arma::uvec cluster_index;
    unsigned int ncate = 0; bool has_cluster = FALSE;
    if(flag){has_cluster = !Rf_isNull(cid[j]);}
    if(has_cluster){
      cluster_index = Rcpp::as<arma::uvec>(cid[j]);
      ncate = arma::max(cluster_index); cluster_index -= 1;
    }else{ has_cluster = 0u; ncate = n; }
    // arma::uvec Z_par_index_u = Z_par_index - 1;
    arma::mat idx = arma::kron(arma::linspace(0,(m-1) * p, m),arma::ones(Z_par_index.n_elem)) + arma::kron(arma::ones(m), Z_par_index);
    arma::uvec uidx = arma::conv_to< arma::uvec>::from(arma::vectorise(idx)-1);
    arma::uvec mask = arma::ones<arma::uvec>(n_alpha);
    mask.elem(uidx).zeros();
    arma::uvec compl_uidx = arma::find(mask==1); 

    arma::vec coef_vec = coef_List[j]; 
    arma::vec coef_vec_full = arma::zeros(n_alpha);
    coef_vec_full(compl_uidx) = coef_vec;
    // score_summary_alpha(coef_vec_full, Y_bin, Z, vA_mat, Vinv_cube,
    //                     VY_mat, Score_reduce_alpha, Hessian_alpha);
    score_summary_alpha(coef_vec_full, Y_bin, Z, Score_reduce_alpha, Hessian_alpha); 
    // vA_lst[j] = vA_mat; Vinv_lst[j] = Vinv_cube; VY_lst[j] = VY_mat;
    arma::mat Score_reduce_alpha_reorg = arma::join_rows(Score_reduce_alpha.cols(uidx), Score_reduce_alpha.cols(compl_uidx));
    arma::mat Hess_reduce_alpha_reorg = arma::join_cols(arma::join_rows(Hessian_alpha.submat(uidx, uidx), 
                                                        Hessian_alpha.submat(uidx, compl_uidx)),
                                                        arma::join_rows(Hessian_alpha.submat(compl_uidx, uidx), 
                                                        Hessian_alpha.submat(compl_uidx, compl_uidx)));

    unsigned int n_par_interest = uidx.n_elem;
    arma::span par_idx(0,(n_par_interest-1)); 
    arma::span compl_par_idx(n_par_interest,(n_alpha-1));
    arma::rowvec Score_reduce_rowsums = arma::sum(Score_reduce_alpha_reorg, 0);
    arma::vec A = Score_reduce_rowsums.cols(par_idx).t();
    arma::mat B1 = arma::pinv(Hess_reduce_alpha_reorg.submat(par_idx, par_idx) - Hess_reduce_alpha_reorg.submat(par_idx, compl_par_idx) *
      arma::pinv(Hess_reduce_alpha_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_alpha_reorg.submat(compl_par_idx, par_idx));
    arma::vec alpha_hat = B1 * A;
    arma::mat U = Score_reduce_alpha_reorg.cols(par_idx) - Score_reduce_alpha_reorg.cols(compl_par_idx) *
      arma::pinv(Hess_reduce_alpha_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_alpha_reorg.submat(par_idx, compl_par_idx).t();
    arma::mat B2(n_par_interest, n_par_interest, arma::fill::zeros);
    if(has_cluster){
      arma::mat within_class_sum(ncate, n_par_interest, arma::fill::zeros);
      for(int j = 0; j < n; j++){
        unsigned int obs_cls_id = cluster_index(j);
        arma::rowvec tmp_vec2 = arma::conv_to< arma::rowvec >::from(U.row(j));
        within_class_sum.row(obs_cls_id) += tmp_vec2;
      }
      B2 = within_class_sum.t() * within_class_sum;
    }else{ B2 = U.t() * U ;}
    arma::mat cov_alpha = B1 * B2 * B1;
    est_alpha_lst[j] = alpha_hat; // restore the score statistics and estimated covariance matrix in lists which are of necessity for some meta-analysis methods
    cov_lst[j] = cov_alpha;
    arma::uvec col_index_u  = col_index_list[j]; col_index_u -= 1;
    est_cov_meta.submat(col_index_u,col_index_u) += arma::pinv(cov_alpha);
    score_stat_alpha(col_index_u) = score_stat_alpha(col_index_u) + arma::pinv(cov_alpha) * alpha_hat;
    if(meta_method == "RE-MetaQCAT*"){
      arma::mat est_inv = arma::pinv(cov_alpha);
      arma::vec tmp = est_inv * arma::vectorise(alpha_hat);
      U_theta = U_theta + 0.5 * arma::as_scalar(tmp.t() * tmp - arma::trace(est_inv));
      V_theta = V_theta + 0.5 * arma::as_scalar(arma::trace(est_inv.t() * est_inv));
    }
  }
  arma::mat cov_alpha_hat = arma::pinv(est_cov_meta);
  arma::vec alpha_hat = cov_alpha_hat * score_stat_alpha;
  arma::vec var_alpha_hat = cov_alpha_hat.diag();
  double score_stat_meta = arma::as_scalar(score_stat_alpha.t() * cov_alpha_hat * score_stat_alpha);
  if(meta_method == "RE-MetaQCAT*"){
    if (!std::isfinite(V_theta) || V_theta <= 0) Rcpp::stop("Degenerate random-effect covariance.");
    score_stat_meta = score_stat_meta + U_theta * U_theta / V_theta;
    return Rcpp::List::create(Rcpp::Named("score.stat.meta") = score_stat_meta,
                              Rcpp::Named("var.alpha.hat") = var_alpha_hat,
                              Rcpp::Named("alpha.hat") = alpha_hat,
                              Rcpp::Named("est.alpha.lst") = est_alpha_lst,
                              Rcpp::Named("cov.lst") = cov_lst);
                              // ,
                              // Rcpp::Named("vA.lst") = vA_lst,
                              // Rcpp::Named("Vinv.lst") = Vinv_lst,
                              // Rcpp::Named("VY.lst") = VY_lst);
  }else{
    double score_pvalue = R::pchisq(score_stat_meta, n_par_interest_alpha, 0, 0);
    return Rcpp::List::create(Rcpp::Named("score.stat.meta") = score_stat_meta,
                              Rcpp::Named("score.pvalue") = score_pvalue,
                              Rcpp::Named("var.alpha.hat") = var_alpha_hat,
                              Rcpp::Named("alpha.hat") = alpha_hat,
                              Rcpp::Named("est.alpha.lst") = est_alpha_lst,
                              Rcpp::Named("cov.lst") = cov_lst);
                              // ,
                              // Rcpp::Named("vA.lst") = vA_lst,
                              // Rcpp::Named("Vinv.lst") = Vinv_lst,
                              // Rcpp::Named("VY.lst") = VY_lst);
  }
}

double score_test_stat_zero_meta_resampling(const Rcpp::List& Z_perm_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, 
                                            const Rcpp::List& col_index_list, const Rcpp::List &Y_bin_list, const Rcpp::List &coef_List,
                                            const arma::vec& Z_par_index, const unsigned int n_par_interest_alpha, const std::string &meta_method){
  
  //alpha.meta.results$vA.lst, alpha.meta.results$Vinv.lst, alpha.meta.results$VY.lst,
  // the total study number
  int total_num = Z_perm_list.length();
  // int total_num_tmp = coef_List.length();
  // int total_num_ttt = Z_perm_list.length();
  double U_theta = 0; double V_theta = 0;
  arma::vec score_stat_alpha = arma::vec(n_par_interest_alpha,arma::fill::zeros);
  arma::mat est_cov_meta(n_par_interest_alpha,n_par_interest_alpha);
  est_cov_meta.zeros();
  Rcpp::List cid; bool flag = cluster_list.isNotNull();
  if(flag){cid = cluster_list.get();}   
  // create null list to save the results for each study
  Rcpp::List score_alpha(total_num);
  Rcpp::List est_cov(total_num);
  for(unsigned int j = 0; j < total_num; ++j){
    arma::mat Y_bin = Y_bin_list[j];
    arma::mat Z = Z_perm_list[j];
    // arma::mat vA_mat = vA_lst[j];
    // arma::cube Vinv_cube = Vinv_lst[j];
    // arma::mat VY_mat = VY_lst[j];
    unsigned int p = Z.n_cols; 
    unsigned int n = Z.n_rows;
    unsigned int m = Y_bin.n_cols;
    unsigned int n_alpha = m * p;
    arma::mat Score_reduce_alpha(n, n_alpha, arma::fill::zeros);
    arma::mat Hessian_alpha(n_alpha, n_alpha, arma::fill::zeros);
    arma::uvec cluster_index;
    unsigned int ncate = 0; bool has_cluster = FALSE;
    if(flag){has_cluster = !Rf_isNull(cid[j]);}
    if(has_cluster){
      cluster_index = Rcpp::as<arma::uvec>(cid[j]);
      ncate = arma::max(cluster_index); cluster_index -= 1;
    }else{ has_cluster = 0u; ncate = n; }

    arma::mat idx = arma::kron(arma::linspace(0,(m-1) * p, m),arma::ones(Z_par_index.n_elem)) + arma::kron(arma::ones(m), Z_par_index);
    arma::uvec uidx = arma::conv_to< arma::uvec>::from(arma::vectorise(idx)-1);
    arma::uvec mask = arma::ones<arma::uvec>(n_alpha);
    mask.elem(uidx).zeros();
    arma::uvec compl_uidx = arma::find(mask==1); 

    arma::vec coef_vec = coef_List[j]; 
    arma::vec coef_vec_full = arma::zeros(n_alpha);
    coef_vec_full(compl_uidx) = coef_vec;
    score_summary_alpha(coef_vec_full, Y_bin, Z, Score_reduce_alpha, Hessian_alpha);
    // for(unsigned int i = 0; i < n; ++i){
    //   arma::vec vA_tmp = vA_mat.col(i);
    //   arma::mat A_i = arma::diagmat(vA_tmp); // m x m
    //   arma::mat t_D_i = arma::kron(A_i, Z.row(i).t()); // n_alpha x m
    //   arma::vec VY = VY_mat.col(i); // m x 1
    //   Score_reduce_alpha.row(i) = (t_D_i * VY).t();
    //   Hessian_alpha += t_D_i * Vinv_cube.slice(i) * t_D_i.t();
    // }
    arma::mat Score_reduce_alpha_reorg = arma::join_rows(Score_reduce_alpha.cols(uidx), Score_reduce_alpha.cols(compl_uidx));
    arma::mat Hess_reduce_alpha_reorg = arma::join_cols(arma::join_rows(Hessian_alpha.submat(uidx, uidx), 
                                                        Hessian_alpha.submat(uidx, compl_uidx)),
                                                        arma::join_rows(Hessian_alpha.submat(compl_uidx, uidx), 
                                                        Hessian_alpha.submat(compl_uidx, compl_uidx)));

    unsigned int n_par_interest = uidx.n_elem;
    arma::span par_idx(0,(n_par_interest-1)); 
    arma::span compl_par_idx(n_par_interest,(n_alpha-1));
    arma::rowvec Score_reduce_rowsums = arma::sum(Score_reduce_alpha_reorg, 0);
    arma::vec A = Score_reduce_rowsums.cols(par_idx).t();
    arma::mat B1 = arma::pinv(Hess_reduce_alpha_reorg.submat(par_idx, par_idx) - Hess_reduce_alpha_reorg.submat(par_idx, compl_par_idx) *
      arma::pinv(Hess_reduce_alpha_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_alpha_reorg.submat(compl_par_idx, par_idx));
    arma::vec alpha_hat = B1 * A;
    arma::mat U = Score_reduce_alpha_reorg.cols(par_idx) - Score_reduce_alpha_reorg.cols(compl_par_idx) *
      arma::pinv(Hess_reduce_alpha_reorg.submat(compl_par_idx, compl_par_idx)) * Hess_reduce_alpha_reorg.submat(par_idx, compl_par_idx).t();
    arma::mat B2(n_par_interest, n_par_interest, arma::fill::zeros);
    if(has_cluster){
      arma::mat within_class_sum(ncate, n_par_interest, arma::fill::zeros);
      for(int j = 0; j < n; j++){
        unsigned int obs_cls_id = cluster_index(j);
        arma::rowvec tmp_vec2 = arma::conv_to< arma::rowvec>::from(U.row(j));
        within_class_sum.row(obs_cls_id) += tmp_vec2;
      }
      B2 = within_class_sum.t() * within_class_sum;
    }else{ B2 = U.t() * U ;}
    arma::mat cov_alpha = B1 * B2 * B1;
    arma::uvec col_index_u  = col_index_list[j]; col_index_u -= 1;
    est_cov_meta.submat(col_index_u,col_index_u)  += arma::pinv(cov_alpha);
    score_stat_alpha(col_index_u) = score_stat_alpha(col_index_u) + arma::pinv(cov_alpha) * alpha_hat;
    if(meta_method == "RE-MetaQCAT*"){
      arma::mat est_inv = arma::pinv(cov_alpha);
      arma::vec tmp = est_inv * arma::vectorise(alpha_hat);
      U_theta = U_theta + 0.5 * arma::as_scalar(tmp.t() * tmp - arma::trace(est_inv));
      V_theta = V_theta + 0.5 * arma::as_scalar(arma::trace(est_inv.t() * est_inv));
    }
  }
  double score_stat_meta = arma::as_scalar(score_stat_alpha.t() * arma::pinv(est_cov_meta) * score_stat_alpha);
  if(meta_method == "RE-MetaQCAT*"){
    if (!std::isfinite(V_theta) || V_theta <= 0) Rcpp::stop("Degenerate random-effect covariance.");
    score_stat_meta = score_stat_meta + U_theta * U_theta / V_theta;
  }
  return (score_stat_meta);
}

// [[Rcpp::export]]
double resample_zero_pvalue(const Rcpp::List& Z_list, const Rcpp::Nullable<Rcpp::List>& cluster_list, const Rcpp::List& col_index_list, 
                            const Rcpp::List &Y_bin_list, const Rcpp::List& coef_List,
                            const arma::vec& Z_par_index, const double score_stat_meta, 
                            const unsigned int n_par_interest_alpha, const unsigned int n_perm, const std::string &meta_method, 
                            const Rcpp::Nullable<Rcpp::String> &permute_strata){
  std::string strata = permute_strata.isNotNull() ? Rcpp::as<std::string>(permute_strata) : "";
  Rcpp::List cid; bool class_flag = cluster_list.isNotNull();
  if(class_flag){cid = cluster_list.get();}
  unsigned int total_num = Z_list.length();
  arma::uvec Z_par_index_u = arma::conv_to< arma::uvec >::from(Z_par_index - 1);
  unsigned int n_one = 0; unsigned int one_acc = 0;
  unsigned int start_nperm = 1;
  unsigned int end_nperm = (n_perm < 100)? n_perm : 100;
  bool flag = true;
  while(flag && (start_nperm <= end_nperm) && (end_nperm <= n_perm)){
    for(unsigned int k = start_nperm; k <= end_nperm; ++k){
      Rcpp::List Z_perm_list(total_num);
      if((!class_flag)||((strata != "within")&(strata != "between"))){
        for(unsigned int i = 0; i < total_num; ++i){
          arma::mat Z_perm = Z_list[i];
          unsigned int n = Z_perm.n_rows;
          arma::uvec perm_idx = arma::randperm(n);
          std::vector<unsigned int> v = arma::conv_to< std::vector<unsigned int> >::from(perm_idx);
          Z_perm.cols(Z_par_index_u) = Z_perm.submat(perm_idx, Z_par_index_u);
          Z_perm_list[i] = Z_perm;
        }
      }else{
        if(strata == "within"){
          for(unsigned int i = 0; i < total_num; ++i){
            bool has_cluster = !Rf_isNull(cid[i]);
            arma::mat Z_perm = Z_list[i];
            unsigned int n = Z_perm.n_rows;
            if(has_cluster){
              arma::uvec label_vec = Rcpp::as<arma::uvec>(cid[i]); label_vec -= 1;
              arma::uvec perm_idx = within_class_perm_index_map(label_vec);
              // std::vector<unsigned int> v = arma::conv_to< std::vector<unsigned int> >::from(perm_idx);
              // std::vector<unsigned int> v1 = arma::conv_to< std::vector<unsigned int> >::from(label_vec);
              Z_perm.cols(Z_par_index_u) = Z_perm.submat(perm_idx, Z_par_index_u);
            }else{
              arma::uvec perm_idx = arma::randperm(n);
              Z_perm.cols(Z_par_index_u) = Z_perm.submat(perm_idx, Z_par_index_u);
            }
            Z_perm_list[i] = Z_perm;
          }
        }else if (strata == "between"){
          for(unsigned int i = 0; i < total_num; ++i){
            bool has_cluster = !Rf_isNull(cid[i]);
            arma::mat Z_perm = Z_list[i];
            unsigned int n = Z_perm.n_rows;
            if(has_cluster){
              arma::uvec label_vec = Rcpp::as<arma::uvec>(cid[i]); label_vec -= 1;
              arma::uvec perm_idx = between_class_perm_index_map(label_vec);
              Z_perm.cols(Z_par_index_u) = Z_perm.submat(perm_idx, Z_par_index_u);
            }else{
              arma::uvec perm_idx = arma::randperm(n);
              Z_perm.cols(Z_par_index_u) = Z_perm.submat(perm_idx, Z_par_index_u);
            }
            Z_perm_list[i] = Z_perm;
          }
        }
      }
      double score_stat_meta_perm = score_test_stat_zero_meta_resampling(Z_perm_list, cluster_list, col_index_list, Y_bin_list, 
                                                                         coef_List, Z_par_index, n_par_interest_alpha, meta_method);
      if (!std::isfinite(score_stat_meta_perm)) Rcpp::stop("Non-finite permutation statistic; p-value not computed.");
      n_one += 1;
      if(score_stat_meta_perm >= score_stat_meta){one_acc += 1;}
    }
    if(one_acc < 1){
      start_nperm = end_nperm + 1;
      end_nperm = std::min(n_perm, (end_nperm + 1) * 100 - 1);
      flag = true;
    }else if(one_acc < 10){
      start_nperm = end_nperm + 1;
      end_nperm = std::min(n_perm, (end_nperm + 1) * 10 - 1);
      flag = true;
    }else{flag = false;}
  }
  if (n_one == 0) Rcpp::stop("No successful permutations; p-value not computed.");
  double score_Rpvalue = (one_acc + 1.0)/(n_one + 1.0);
  return(score_Rpvalue);
}
