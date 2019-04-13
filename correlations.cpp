#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat remove_col(arma::mat x, int which){
  x.shed_col(which);
  return(x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec select_col(arma::mat x, int which){
  arma::vec z = x.col(which);
  return(z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat select_row(arma::mat x, int which){
  arma::mat z = x.row(which);
  return(z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat select_row_remove_col(arma::mat x, int which_row, int which_col){
  arma::mat z = x.row(which_row);
  z.shed_col(which_col);
  return(z);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec offdiag_extract(arma::mat& A, int k) {
  return A.diag(k);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double select_row_select_col(arma::mat x, int which_row, int which_col){
  arma::vec z = x.row(which_row);
  double y = z[which_col];
  return(y);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat slice_extract(arma::cube X, int which){
  arma::mat temp = X.slice(which);
  return(temp);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat remove_row(arma::mat x, int which){
  x.shed_row(which);
  return(x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mv_ordinal(arma::mat X, arma::mat model_mat, int nu, int samples) {
  
  // number of rows
  float n = X.n_rows;
  
  // number of columns
  int p = X.n_cols;
  
  // number of predictors
  int k = model_mat.n_cols;
  
  mat vb(p, p);
  
  vec n_levels = unique(X.col(0));
  
  Range thresh_sampled =  seq(min(X.col(0)) + 1, max(X.col(0)) - 1);
  
  // for R in MH step
  colvec mh_step(p);
  mat R_top(n,p);
  mat  R_bottom(n,p);
  
  // intialize beta
  cube  b(p, k, 1, arma::fill::zeros);
  
  // intitialize variances
  mat ss(p, 1);
  
  // intialize thresholds
  cube  thresh_mat(p, n_levels.size() + 1, 1, arma::fill::zeros);
  
  // current thresholds
  cube c_thresh_mat(p, n_levels.size() + 1, 1, arma::fill::zeros); 
  
  // latent data
  arma::mat zstar(n, p, arma::fill::zeros); 
  
  // intialize correlation matrix
  cube cor_mat(p,p, 1, arma::fill::zeros);
  cor_mat.slice(0) = arma::eye<arma::mat>(p, p);
  
 // store correlations    
 cube cor_save(p, p, samples);
 
 // save coefficients
 cube beta_save(p, k, samples);
  
  int size = thresh_sampled.size();
  
  int arr[] = {0, size + 2};
  
    for(int i = 0; i < p; i ++ ){
      for(int j = 0; j < 2; j ++){
        
        if(j == 0){
          
          c_thresh_mat(i,arr[j],0) = -arma::datum::inf;
          
          thresh_mat(i,arr[j], 0) =  -arma::datum::inf;
        } if(j == 1){
          
          c_thresh_mat(i,arr[j],0) = arma::datum::inf;
          
          thresh_mat(i, arr[j],0) = arma::datum::inf;
        }
      }
    }
    
  for(int k = 0; k < samples; k ++){
    
  
    mat m = model_mat * b.slice(0).t();
    
    for(int i = 0; i < max(X.col(0)) - 1; i ++){
      
      for(int j = 0; j <= p - 1; j ++){
        
        thresh_mat(j, i + 2, 0) = R::qnorm5(sum(X.col(j) <= thresh_sampled[i] ) / n,
                                             - R::qnorm5(sum(X.col(j) == 1) / n, 0, 1, TRUE, FALSE), 1, TRUE, FALSE);
        
      }
    }
    
    // define for scoping
    mat  mm(n,1);
    
    for(int i = 0; i < p ; i++){
      
      // row i and remove column i--i.e., mat[ i, -i]
      arma::mat mat_i_not_i = select_row_remove_col(cor_mat.slice(0), i, i);
      
      
      // remove row i and remove column i--i.e., mat[ - i, -i]
      arma::mat mat_not_i_not_i = inv(remove_row(remove_col(cor_mat.slice(0), i), i));
      
      // updated: subset of multivariate location
      mm.col(0) = trans(trans(m.col(i)) + mat_i_not_i * mat_not_i_not_i *  trans(remove_col(zstar, i) - remove_col(m, i)));
      
      // updated: subset of multivariate scale
      arma::mat ss =  select_row(cor_mat.slice(0).col(i) , i) -  mat_i_not_i  * mat_not_i_not_i * trans(mat_i_not_i);
      
      
      for(int j = 0; j < n ; j ++ ){
        zstar(j, i) = R::qnorm(R::runif(
          // minimum
          R::pnorm(thresh_mat.slice(0)(i , (X.col(i)[j] - 1)), mm.col(0)(j), sqrt(ss(0)), TRUE, FALSE),
          // maximum
          R::pnorm(thresh_mat.slice(0)(i , X.col(i)[j]), mm.col(0)(j), sqrt(ss(0)), TRUE, FALSE)),
          // location and scale
          mm.col(0)(j), sqrt(ss(0)), TRUE, FALSE) ;
        
        
      }
    }
    
    for(int i = 0; i < max(X.col(0)) - 1; i ++){
      // note 0 = 1; 1 - max gives largest threshold
      for(int j = 0; j < p ; j ++){
        c_thresh_mat(j, thresh_sampled[i], 0) = R::qnorm5(R::runif(R::pnorm(0, thresh_mat.slice(0)(j , thresh_sampled[i]), .05, TRUE, FALSE), 1),
                                                          thresh_mat.slice(0)(j , thresh_sampled[i]),   0.05, TRUE, FALSE);


      }
    }
    
    for(int i = 0; i < p ; i ++){
      for(int j = 0; j < n ; j ++ ){
        
        R_top(j, i) =  (R::pnorm5(c_thresh_mat(i , X.col(i)[j], 0) - mm(j), 0, 1, TRUE, FALSE) -
                          R::pnorm5(c_thresh_mat(i , X.col(i)[j] - 1, 0) - mm(j), 0, 1, TRUE, FALSE));
        
        R_bottom(j, i) = (R::pnorm5(thresh_mat(i , X.col(i)[j], 0) - mm(j), 0, 1, TRUE, FALSE) -
                            R::pnorm5(thresh_mat(i , X.col(i)[j] - 1, 0) - mm(j), 0, 1, TRUE, FALSE));
        
        
      }
    }
    
    for(int i = 0; i < p; i ++){
      mh_step(i) = accu(log(R_top.col(i)) - log(R_bottom.col(i))) ;
      
    }
    
    
    for(int j = 0; j < p; j ++){
      if( mh_step[j] > log(R::runif(0, 1) )){
        thresh_mat.slice(0).row(j) = c_thresh_mat.slice(0).row(j);
        
      }  
      
    }
    // sample MVN
    arma::mat vb(inv(inv(cor_mat.slice(0) / n)));
    arma::mat mn = arma::mat(vb * trans((model_mat.t() * zstar *  trans(inv(cor_mat.slice(0))))));
    // sampled beta (current)
    arma::mat c_beta = mn.t() + (b.slice(0).randn().t() * chol(vb));
  
    // residuals
    arma::mat e = zstar - (trans(c_beta.t() * model_mat.t()));
    
    // residual scatter matrix  
    arma::mat  vt =  trans(e) * e ;
    
    //  sample covariance matrix
    arma::mat  Sigma = iwishrnd(vt, n + nu - p - 1);
  
    // current correlation matrix
    arma::mat c_cor_mat = diagmat(1 / sqrt(Sigma.diag())) * Sigma * diagmat(1 / sqrt(Sigma.diag()));
    
     // store correlation matrices
    cor_save.slice(k) = c_cor_mat;
    
    // store betas 
    beta_save.slice(k) = c_beta.t();
    
     // return to top 
    cor_mat.slice(0) = c_cor_mat;
    b.slice(0) = c_beta.t();
    
  
  }
  
  
  List ret;
  ret["betas"] = beta_save;
  ret["correlations"] = cor_save;
  return ret;
  
}
