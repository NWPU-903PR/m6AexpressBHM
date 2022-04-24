#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//  Function that computes multivariate normal pdf for a single point
double dmvnrm_lf(const arma::vec &x,
                 const arma::vec &mean,
                 const arma::mat &sigma){
  int xdim = x.n_elem;
  arma::mat rooti;
  arma::vec z;
  double out;
  double tmp;
  double rootisum, constants;
  const double log2pi = std::log(2.0 * M_PI);
  constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  arma::mat sigma_copy = sigma;
  // Do the Cholesky decomposition
  arma::mat csigma;
  bool success = false;
  while(success == false)
  {
    success = chol(csigma, sigma_copy, "upper");
    
    if(success == false)
    {
      sigma_copy += arma::eye(sigma.n_rows,sigma.n_rows) * 1e-4;
      //Rcpp::Rcout << "Inv error dmvnorm" << std::endl;
    }
  }
  rooti = arma::trans(arma::inv(trimatu(csigma)));
  rootisum = arma::sum(log(rooti.diag()));
  z = rooti * arma::trans(x.t() - mean.t());
  tmp = constants - 0.5 * arma::sum(z%z) + rootisum;
  out = exp(tmp);
  return(out);
}
// Function that generates random multivariate normal points
arma::mat rmvnormal(const int &n,
                    const arma::vec &mu,
                    const arma::mat &sigma) {
  const int d = mu.size();
  arma::vec mu_tmp = mu;
  
  // Create matrix of standard normal random values
  arma::mat ans(n, d, arma::fill::randn);
  arma::mat sigma_copy = sigma;
  
  // Do the Cholesky decomposition
  arma::mat csigma;
  bool success = false;
  
  while(success == false)
  {
    success = chol(csigma, sigma_copy);
    
    if(success == false)
    {
      sigma_copy += arma::eye(sigma.n_rows,sigma.n_rows) * 1e-4;
      //Rcpp::Rcout << "Inv error rmvnorm" << std::endl;
    }
  }
  
  // Do the transformation
  ans = ans * csigma;
  ans.each_row() += mu.t(); // Add mu to each row in transformed ans
  
  return ans;
}
// Update beta and b0 together
arma::vec update_beta(const arma::rowvec &beta_cur,
                      const arma::rowvec &counts,
                      double &alpha_cur,
                      const arma::rowvec &b0_cur,
                      const arma::vec &log_offset,
                      const arma::mat kX,
                      const arma::mat kZ,
                      const arma::vec beta_mean,
                      const arma::vec beta_sd,
                      const double &sigma2_cur,
                      const int &n_beta,
                      const int &n_b0,
                      const double &sd0,
                      const int &n_sample){
  arma::vec all_beta_curtmp = join_cols(beta_cur.t(),b0_cur.t());
  arma::vec all_beta_prop(n_beta+n_b0), mean_cur(n_sample), mean_prop(n_sample);
  arma::vec mean_wls_cur(n_beta+n_b0),mean_wls_prop(n_beta+n_b0);
  arma::mat W_mat(n_sample, n_sample), R_mat(n_b0 + n_beta, n_b0 + n_beta), R_mat_i(n_b0 + n_beta, n_b0 + n_beta);
  arma::mat W_mat_prop(n_sample, n_sample);
  arma::mat cov_mat_cur(n_b0 + n_beta, n_b0 + n_beta), cov_mat_prop(n_b0 + n_beta, n_b0 + n_beta);
  arma::vec y_tilde(n_sample),y_tilde_prop(n_sample),eta(n_sample), eta_prop(n_sample);
  arma::mat design_mat=join_horiz(kX, kZ);
  arma::vec prior_mean_b0(n_b0);
  prior_mean_b0.zeros();
  arma::vec prior_mean_betas(n_b0+n_beta), R_mat_diag(n_b0+n_beta);
  arma::vec prior_var_beta(n_beta);
  for(int k=0; k < n_beta; k++){
    prior_var_beta(k) = pow(beta_sd(k),2);
  }
  arma::vec sigma2_cur_rw(n_b0);
  sigma2_cur_rw.fill(sigma2_cur);
  R_mat_diag=join_cols(prior_var_beta,sigma2_cur_rw);
  prior_mean_betas = join_cols(beta_mean,prior_mean_b0);
  
  double ll_prop, ll_cur, mh_prop, mh_cur, rho_cur = 1.0 / alpha_cur;
  eta = design_mat*all_beta_curtmp+log_offset;
  mean_cur= arma::exp(eta);
  y_tilde = eta + arma::trans(counts - mean_cur.t()) % arma::exp(-eta) - log_offset;
  arma::mat W_mat_i = arma::diagmat(1.0 / (alpha_cur + arma::exp(-eta)));
  R_mat_i.zeros();
  R_mat_i.diag() = 1.0 / R_mat_diag;
  R_mat.zeros();
  R_mat.diag() = R_mat_diag;
  arma::mat ds_W_mat_i = design_mat.t() * W_mat_i;
  arma::mat csigma = R_mat_i + ds_W_mat_i * design_mat;
  bool success = false;
  
  success = arma::inv_sympd(cov_mat_cur, csigma);
  
  if(success == false)
  {
    return(all_beta_curtmp);
  }
  
  mean_wls_cur = cov_mat_cur * (ds_W_mat_i * y_tilde);
  
  all_beta_prop = arma::trans(rmvnormal(1, mean_wls_cur, sd0*cov_mat_cur));
  eta_prop = design_mat * all_beta_prop + log_offset;
  mean_prop = arma::exp(eta_prop);
  
  y_tilde_prop = eta_prop + arma::trans(counts - mean_prop.t()) % arma::exp(-eta_prop) - log_offset;
  
  arma::mat W_mat_prop_i = arma::diagmat(1.0 / (alpha_cur + arma::exp(-eta_prop)));
  
  arma::mat ds_W_mat_prop_i = design_mat.t() * W_mat_prop_i;
  arma::mat csigma2 = R_mat_i + ds_W_mat_prop_i * design_mat;
  
  success = arma::inv_sympd(cov_mat_prop, csigma2);
  
  if(success == false)
  {
    return(all_beta_curtmp);
  }
  
  mean_wls_prop = cov_mat_prop * (ds_W_mat_prop_i * y_tilde_prop);
  
  ll_cur = arma::sum(counts * arma::log(mean_cur) - (counts + rho_cur) * arma::log(1.0 + mean_cur * alpha_cur));
  ll_prop = arma::sum(counts * arma::log(mean_prop) - (counts + rho_cur) * arma::log(1.0 + mean_prop * alpha_cur));
  
  mh_cur = ll_cur +
    log(dmvnrm_lf(all_beta_curtmp, prior_mean_betas, R_mat)) -
    log(dmvnrm_lf(all_beta_curtmp, mean_wls_prop, cov_mat_prop));
  
  mh_prop = ll_prop +
    log(dmvnrm_lf(all_beta_prop, prior_mean_betas, R_mat)) -
    log(dmvnrm_lf(all_beta_prop, mean_wls_cur, cov_mat_cur));
  
  if((R::runif(0, 1) < exp(mh_prop - mh_cur))){
    all_beta_curtmp = all_beta_prop;
  }
  return all_beta_curtmp;
  
}
// Update alpha
double update_rho(const arma::rowvec &beta_cur,
                  const arma::rowvec &counts,
                  const arma::rowvec &b0_cur,
                  const double &alpha_cur,
                  const double &log_dispFitted,
                  const double &varPrior,
                  const arma::vec &log_offset,
                  const arma::mat kX,
                  const arma::mat kZ,
                  const double sd1,
                  const int &n_beta,
                  const int &n_sample){
  arma::vec mean_cur(n_sample);
  double ll_prop, ll_cur, mh_prop, mh_cur, rho_prop, rho_cur=1.0/alpha_cur, alpha_prop;
  mean_cur = arma::exp(kX*beta_cur.t()+kZ*b0_cur.t()+log_offset);
  alpha_prop = (exp(log(alpha_cur)+ R::rnorm(0,sd1)));
  rho_prop = 1.0/alpha_prop;
  ll_cur = arma::sum(arma::lgamma(rho_cur + counts)) - arma::sum((counts + rho_cur) * arma::log(1.0 + mean_cur * alpha_cur)) - n_sample * lgamma(rho_cur) + arma::sum(counts * log(alpha_cur));
  ll_prop = arma::sum(arma::lgamma(rho_prop + counts)) - arma::sum((counts + rho_prop) * arma::log(1.0 + mean_cur * alpha_prop)) - n_sample * lgamma(rho_prop) + arma::sum(counts * log(alpha_prop));
  mh_cur = ll_cur + R::dnorm4(log(alpha_cur), log_dispFitted, varPrior, 1);
  mh_prop = ll_prop + R::dnorm4(log(alpha_prop), log_dispFitted, varPrior, 1);
  if(R::runif(0,1) < exp(mh_prop-mh_cur)){
    return alpha_prop;
    
  }
  else{
    return alpha_cur;
  }
}
// Update sigma2
double sigm2sugar(const double & est_alph, 
                  const double & est_beta, 
                  const arma::rowvec & b0_cur) {   
  double sigma2_cur;
  const int d = b0_cur.size();
  double a_rand_int_post = est_alph + d / 2.0;
  double b_rand_int_post = est_beta+arma::dot(b0_cur.t(), b0_cur.t()) / 2.0;
  sigma2_cur = 1.0 / (R::rgamma(a_rand_int_post, 1.0 / b_rand_int_post));
  return sigma2_cur;
}

// [[Rcpp::export]]
Rcpp::List MH_MC_iter(const arma::rowvec counts,
                      const arma::mat kX,
                      const arma::mat kZ,
                      const arma::rowvec betas,
                      const arma::rowvec b0s,
                      const arma::vec beta_mean,
                      const arma::vec beta_sd,
                      const arma::vec log_offset,
                      const double alpha,
                      const double sigma2,
                      const double disp_variance,
                      const double log_dispersion_fit,
                      const double est_alpha,
                      const double est_beta,
                      double sd0,
                      double sd1,
                      const int n_b0,
                      const int n_beta,
                      const int n_sample,
                      const int it_num){
  
  arma::mat betas_b0s_sample(it_num,n_beta+n_b0);
  arma::vec disp_sample(it_num);
  arma::vec sigma2_sample(it_num);
  arma::rowvec betas_last(n_beta), b0s_last(n_b0);
  arma::rowvec beta_b0s_cur(n_beta+n_b0), beta_b0s_last(n_beta+n_b0);
  betas_b0s_sample.row(0) = arma::trans(join_cols(betas.t(),b0s.t()));
  disp_sample.zeros();
  disp_sample(0) = alpha;
  sigma2_sample(0) = sigma2;
  betas_last = betas;
  b0s_last = b0s;
  for(int i=1; i<it_num; i++){
    beta_b0s_cur = arma::trans(update_beta(betas_last,
                                        counts,
                                        disp_sample(i-1),
                                        b0s_last,
                                        log_offset,
                                        kX, 
                                        kZ,
                                        beta_mean,
                                        beta_sd,
                                        sigma2_sample(i-1),
                                        n_beta,
                                        n_b0,
                                        sd0,
                                        n_sample));
    
    beta_b0s_last = beta_b0s_cur;
    betas_b0s_sample.row(i) = beta_b0s_cur;
    betas_last = beta_b0s_last.cols(0,n_beta-1);
    b0s_last = beta_b0s_last.cols(n_beta,((n_b0 + n_beta)-1));
    disp_sample(i) = update_rho(betas_last,
                                    counts,
                                  b0s_last,
                          disp_sample(i-1),
                        log_dispersion_fit,
                             disp_variance,
                                log_offset,
                                        kX, 
                                        kZ,
                                        sd1,
                                     n_beta,
                                  n_sample);
    
    sigma2_sample(i) = sigm2sugar(est_alpha, est_beta, b0s_last);
  }
  
  return Rcpp::List::create(Rcpp::Named("betas_b0s_est")= betas_b0s_sample,
                            Rcpp::Named("alphas")=disp_sample,
                            Rcpp::Named("sigma2")=sigma2_sample);
}
